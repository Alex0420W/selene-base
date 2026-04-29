"""Reproject every available raw raster onto the common 240 m grid.

For each dataset we know how to load:

- check whether the raw bytes are present (skip with a clear message
  if not),
- load with the matching :mod:`selene_base.data.load` helper,
- reproject to the south-polar stereographic grid via
  :func:`selene_base.data.reproject.reproject_to_grid`,
- cache as a COG in ``data/processed/``.

Robbins is intentionally not in this loop — it is vector data and is
rasterised by the hazard criterion directly in week 3.

Filled in week 2.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import typer
import xarray as xr
import yaml

from selene_base.data.load import (
    load_crater_catalog,
    load_diviner,
    load_illumination,
    load_lend,
    load_lola_ldem,
)
from selene_base.data.rasterize import rasterize_crater_density
from selene_base.data.reproject import cache_processed, reproject_to_grid

DEFAULT_REGION_CONFIG = Path("config/region_southpole.yaml")
DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PROCESSED_DIR = Path("data/processed")


@dataclass(frozen=True)
class DatasetSpec:
    """How to detect, load, and resample one raw dataset."""

    name: str
    raw_check: Path
    loader: Callable[[], xr.DataArray]
    resampling: str
    note: str


def _lola() -> xr.DataArray:
    return load_lola_ldem()


def _illumination() -> xr.DataArray:
    return load_illumination()


def _diviner_tmax() -> xr.DataArray:
    return load_diviner().tbol_max


def _diviner_tmin() -> xr.DataArray:
    return load_diviner().tbol_min


def _lend() -> xr.DataArray:
    return load_lend()


# Resampling rationale documented per dataset; see also README §Methodology.
RASTER_DATASETS: tuple[DatasetSpec, ...] = (
    DatasetSpec(
        name="lola",
        raw_check=DEFAULT_RAW_DIR / "lola" / "ldem_80s_80m.img",
        loader=_lola,
        resampling="bilinear",
        note="elevation is a smooth continuous surface; bilinear is correct.",
    ),
    DatasetSpec(
        name="illumination",
        raw_check=DEFAULT_RAW_DIR / "illumination" / "avgvisib_65s_240m_201608.img",
        loader=_illumination,
        resampling="bilinear",
        note="continuous illumination percentage; bilinear preserves it.",
    ),
    DatasetSpec(
        name="diviner_tmax",
        raw_check=DEFAULT_RAW_DIR / "diviner" / "diviner_tbol_max_sp.tif",
        loader=_diviner_tmax,
        resampling="bilinear",
        note="continuous Tbol field; bilinear is correct.",
    ),
    DatasetSpec(
        name="diviner_tmin",
        raw_check=DEFAULT_RAW_DIR / "diviner" / "diviner_tbol_min_sp.tif",
        loader=_diviner_tmin,
        resampling="bilinear",
        note="continuous Tbol field; bilinear is correct.",
    ),
    DatasetSpec(
        name="lend",
        raw_check=DEFAULT_RAW_DIR / "lend" / "lend_csetn_sp.img",
        loader=_lend,
        resampling="bilinear",
        note="coarse neutron map; bilinear at 240 m smooths cleanly.",
    ),
)


@dataclass
class PreprocessResult:
    """One row of the ``selene preprocess`` summary table."""

    name: str
    status: str  # "cached", "skipped", "missing"
    output_path: Path | None
    bytes_written: int


def load_region_config(path: Path) -> dict[str, object]:
    """Read the south-polar grid definition (CRS, bounds, resolution)."""
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def run(
    *,
    region_config: Path = DEFAULT_REGION_CONFIG,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    overwrite: bool = False,
    datasets: tuple[DatasetSpec, ...] | None = None,
    echo: Callable[[str], None] = typer.echo,
) -> list[PreprocessResult]:
    """Reproject + cache every available raster dataset.

    Args:
        region_config: YAML defining ``crs``, ``bounds_m``, ``resolution_m``.
        processed_dir: Output directory for cached COGs.
        overwrite: Re-cache even when the COG already exists.
        datasets: Override the dataset list (used by tests).
        echo: Logging sink, defaults to ``typer.echo``.

    Returns:
        One :class:`PreprocessResult` per dataset, in input order.
    """
    cfg = load_region_config(region_config)
    target_crs = str(cfg["crs"])
    bounds_m = tuple(cfg["bounds_m"])  # type: ignore[arg-type]
    resolution_m = float(cfg["resolution_m"])

    if len(bounds_m) != 4:
        raise ValueError(f"bounds_m must have 4 entries, got {bounds_m!r}")

    if datasets is None:
        datasets = RASTER_DATASETS

    results: list[PreprocessResult] = []
    for spec in datasets:
        if not spec.raw_check.exists():
            echo(f"[skip] {spec.name}: raw file not present ({spec.raw_check})")
            results.append(PreprocessResult(spec.name, "missing", None, 0))
            continue

        echo(f"[load] {spec.name} from {spec.raw_check.parent}")
        src = spec.loader()
        echo(f"[warp] {spec.name} ({spec.resampling})")
        warped = reproject_to_grid(
            src,
            target_crs=target_crs,
            bounds_m=bounds_m,  # type: ignore[arg-type]
            resolution_m=resolution_m,
            resampling=spec.resampling,
        )
        out_path = cache_processed(warped, spec.name, processed_dir, overwrite=overwrite)
        size = out_path.stat().st_size
        echo(f"[done] {spec.name} -> {out_path} ({size:,} bytes)")
        results.append(PreprocessResult(spec.name, "cached", out_path, size))

    # ------- Robbins → crater-density raster (vector → grid) -------
    robbins_path = DEFAULT_RAW_DIR / "robbins" / "robbins_southpole.csv.gz"
    crater_cog = processed_dir / "crater_density_southpole_240m.tif"
    if not robbins_path.exists():
        echo(f"[skip] crater_density: Robbins catalog not present ({robbins_path})")
        results.append(PreprocessResult("crater_density", "missing", None, 0))
    elif crater_cog.exists() and not overwrite:
        size = crater_cog.stat().st_size
        echo(f"[skip] crater_density: {crater_cog.name} already cached")
        results.append(PreprocessResult("crater_density", "cached", crater_cog, size))
    else:
        # Need a target grid to rasterise onto — reuse the LOLA COG if present,
        # otherwise build a zeros grid in the right CRS.
        target_path = processed_dir / "lola_southpole_240m.tif"
        if target_path.exists():
            import rioxarray  # local import to keep top-level imports light

            target_grid = rioxarray.open_rasterio(target_path, masked=True).squeeze(
                "band", drop=True
            )
        else:
            import numpy as np
            import rioxarray  # noqa: F401  (registers .rio accessor)
            from rasterio.transform import from_origin

            xmin, ymin, xmax, ymax = bounds_m  # type: ignore[misc]
            width = int(round((xmax - xmin) / resolution_m))
            height = int(round((ymax - ymin) / resolution_m))
            target_grid = xr.DataArray(
                np.zeros((height, width), dtype=np.float32),
                dims=("y", "x"),
                coords={
                    "y": np.linspace(ymax - resolution_m / 2, ymin + resolution_m / 2, height),
                    "x": np.linspace(xmin + resolution_m / 2, xmax - resolution_m / 2, width),
                },
            ).rio.write_crs(target_crs, inplace=False)
            target_grid.rio.write_transform(
                from_origin(xmin, ymax, resolution_m, resolution_m), inplace=True
            )

        echo("[load] robbins crater catalog")
        craters = load_crater_catalog(robbins_path)
        echo(f"[rasterise] crater density (n={len(craters):,}, radius=3 km)")
        density = rasterize_crater_density(
            craters, target_grid, radius_km=3.0, diameter_col="diam_km"
        )
        out_path = cache_processed(density, "crater_density", processed_dir, overwrite=overwrite)
        size = out_path.stat().st_size
        echo(f"[done] crater_density -> {out_path} ({size:,} bytes)")
        results.append(PreprocessResult("crater_density", "cached", out_path, size))

    return results


def format_summary(results: list[PreprocessResult]) -> str:
    """Render a fixed-width summary of preprocess results."""
    header = f"{'dataset':<14} {'status':<9} {'size':>12}  path"
    rows = []
    for r in results:
        size_str = f"{r.bytes_written:,}" if r.bytes_written else "-"
        path_str = str(r.output_path) if r.output_path else "-"
        rows.append(f"{r.name:<14} {r.status:<9} {size_str:>12}  {path_str}")
    return "\n".join([header, *rows])
