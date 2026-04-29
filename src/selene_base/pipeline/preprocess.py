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
    """How to detect, load, and resample one raw dataset.

    ``raw_check`` is a tuple of candidate file paths; the first entry
    that exists is the one used. This lets the same DatasetSpec serve
    both the full-resolution PDS3 product and the downsampled GeoTIFF
    that ships in the sample-data tarball.
    """

    name: str
    raw_check: tuple[Path, ...]
    loader: Callable[[Path], xr.DataArray]
    resampling: str
    note: str

    def resolve(self) -> Path | None:
        for candidate in self.raw_check:
            if candidate.exists():
                return candidate
        return None


def _load_lola_path(path: Path) -> xr.DataArray:
    # Sample-data TIF carries elevation already in metres; the PDS3
    # IMG path applies the 0.5 m scale via load_lola_ldem.
    if path.suffix.lower() == ".tif":
        from selene_base.data.load import load_raster

        return load_raster(path).rename("elevation_m")
    return load_lola_ldem(path)


def _load_illumination_path(path: Path) -> xr.DataArray:
    if path.suffix.lower() == ".tif":
        # The sample-data TIF was produced AFTER load_illumination's
        # DN -> fraction scaling, so values are already in [0, 1].
        from selene_base.data.load import load_raster

        return load_raster(path).rename("illumination_fraction")
    return load_illumination(path)


def _load_lend_path(path: Path) -> xr.DataArray:
    return load_lend(path)


# Resampling rationale documented per dataset; see also README §Methodology.
RASTER_DATASETS: tuple[DatasetSpec, ...] = (
    DatasetSpec(
        name="lola",
        raw_check=(
            # PDS3 driver opens the .lbl (which references the .img); put it
            # first so spec.resolve() returns the file rasterio can actually
            # parse. Sample-data path is a plain GeoTIFF.
            DEFAULT_RAW_DIR / "lola" / "ldem_80s_80m.lbl",
            DEFAULT_RAW_DIR / "lola" / "sample_lola.tif",
        ),
        loader=_load_lola_path,
        resampling="bilinear",
        note="elevation is a smooth continuous surface; bilinear is correct.",
    ),
    DatasetSpec(
        name="illumination",
        raw_check=(
            DEFAULT_RAW_DIR / "illumination" / "avgvisib_65s_240m_201608.lbl",
            DEFAULT_RAW_DIR / "illumination" / "sample_illumination.tif",
        ),
        loader=_load_illumination_path,
        resampling="bilinear",
        note="continuous illumination percentage; bilinear preserves it.",
    ),
    # Diviner PRP is rasterised separately below — it's a triangular
    # mesh, not a regular raster, and goes through the
    # selene_base.data.diviner_prp loader.
    DatasetSpec(
        name="lend",
        raw_check=(DEFAULT_RAW_DIR / "lend" / "lend_csetn_sp.img",),
        loader=_load_lend_path,
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
        resolved = spec.resolve()
        if resolved is None:
            echo(
                f"[skip] {spec.name}: raw file not present "
                f"(checked {[str(p) for p in spec.raw_check]})"
            )
            results.append(PreprocessResult(spec.name, "missing", None, 0))
            continue

        echo(f"[load] {spec.name} from {resolved.parent}")
        src = spec.loader(resolved)
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

    # ------- Diviner PRP → temp_avg + temp_max + ice_depth rasters -------
    diviner_xml = DEFAULT_RAW_DIR / "diviner" / "dlre_prp_south.xml"
    diviner_tab = DEFAULT_RAW_DIR / "diviner" / "dlre_prp_south.tab"
    diviner_layers = ("temp_avg", "temp_max", "ice_depth")
    diviner_cog_paths = {
        layer: processed_dir / f"diviner_{layer}_southpole_240m.tif" for layer in diviner_layers
    }
    diviner_caches_present = all(p.exists() for p in diviner_cog_paths.values())

    if not diviner_caches_present and not (diviner_xml.exists() and diviner_tab.exists()):
        echo(f"[skip] diviner_prp: source not present ({diviner_xml.parent})")
        for layer in diviner_layers:
            results.append(PreprocessResult(f"diviner_{layer}", "missing", None, 0))
    elif diviner_caches_present and not overwrite:
        for layer, path in diviner_cog_paths.items():
            size = path.stat().st_size
            echo(f"[skip] diviner_{layer}: {path.name} already cached")
            results.append(PreprocessResult(f"diviner_{layer}", "cached", path, size))
    else:
        from selene_base.data.diviner_prp import load_diviner_prp

        target_path = processed_dir / "lola_southpole_240m.tif"
        if not target_path.exists():
            echo(
                "[skip] diviner_prp: needs the LOLA COG as a target grid; "
                "ensure LOLA preprocess ran first"
            )
            for layer in diviner_layers:
                results.append(PreprocessResult(f"diviner_{layer}", "missing", None, 0))
        else:
            import rioxarray  # local import keeps top-level imports light

            target_grid = rioxarray.open_rasterio(target_path, masked=True).squeeze(
                "band", drop=True
            )
            echo("[load] diviner PRP (PDS4 character table)")
            layers = load_diviner_prp(
                raw_dir=diviner_xml.parent,
                target_grid=target_grid,
                cache_dir=processed_dir,
                overwrite=overwrite,
                echo=echo,
            )
            for layer in diviner_layers:
                path = diviner_cog_paths[layer]
                size = path.stat().st_size if path.exists() else 0
                results.append(PreprocessResult(f"diviner_{layer}", "cached", path, size))
            del layers  # release file handles

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
