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
    load_diviner,
    load_illumination,
    load_lend,
    load_lola_ldem,
)
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

    echo("[skip] robbins: vector data, rasterised inside the hazard criterion (week 3)")
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
