"""Warp rasters onto the common south-polar analysis grid.

The downstream pipeline assumes every score map shares one grid: the
stereographic projection defined in ``config/region_southpole.yaml`` at
240 m / pixel over a ±304 km half-extent. This module is the only place
that touches CRS transforms, plus the COG cache writer used by
``selene preprocess``.

Filled in week 2.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr
from rasterio.enums import Resampling
from rasterio.transform import from_origin

Bounds = tuple[float, float, float, float]

_RESAMPLING_NAMES: dict[str, Resampling] = {
    "nearest": Resampling.nearest,
    "bilinear": Resampling.bilinear,
    "cubic": Resampling.cubic,
    "min": Resampling.min,
    "max": Resampling.max,
    "mean": Resampling.average,
    "median": Resampling.med,
}


def reproject_to_grid(
    src: xr.DataArray,
    *,
    target_crs: str,
    bounds_m: Bounds,
    resolution_m: float,
    resampling: str = "bilinear",
) -> xr.DataArray:
    """Reproject ``src`` onto the common south-polar grid.

    The output grid is rectangular in ``target_crs`` with origin at the
    top-left of ``bounds_m`` and square pixels of ``resolution_m``.
    Nodata is preserved as NaN in the output so downstream criterion
    code can rely on ``np.isnan`` masks.

    Args:
        src: Source DataArray; must have a CRS available via
            ``src.rio.crs``.
        target_crs: Target CRS as a PROJ string or EPSG code.
        bounds_m: ``(xmin, ymin, xmax, ymax)`` in target-CRS units.
        resolution_m: Pixel size in target-CRS units (square pixels).
        resampling: One of ``"nearest"``, ``"bilinear"``, ``"cubic"``,
            ``"min"``, ``"max"``, ``"mean"``, ``"median"``. Choose by
            data character — bilinear for smooth continuous fields,
            nearest for categorical, min/max for hazard envelopes.

    Returns:
        DataArray on the target grid with dims ``("y", "x")`` and CRS
        set to ``target_crs``.

    Raises:
        ValueError: If ``src`` has no CRS, if ``bounds_m`` is degenerate,
            if ``resolution_m`` is non-positive, or if ``resampling`` is
            unknown.
    """
    if src.rio.crs is None:
        raise ValueError(
            "reproject_to_grid: src has no CRS. Set one before reprojecting "
            "(e.g. via src.rio.write_crs(...))."
        )
    if resolution_m <= 0:
        raise ValueError(f"resolution_m must be positive, got {resolution_m!r}")
    xmin, ymin, xmax, ymax = bounds_m
    if not (xmax > xmin and ymax > ymin):
        raise ValueError(f"bounds_m must satisfy xmax>xmin and ymax>ymin, got {bounds_m!r}")
    if resampling not in _RESAMPLING_NAMES:
        raise ValueError(
            f"unknown resampling {resampling!r}; choose one of {sorted(_RESAMPLING_NAMES)}"
        )

    width = int(round((xmax - xmin) / resolution_m))
    height = int(round((ymax - ymin) / resolution_m))
    transform = from_origin(xmin, ymax, resolution_m, resolution_m)

    out = src.rio.reproject(
        dst_crs=target_crs,
        shape=(height, width),
        transform=transform,
        resampling=_RESAMPLING_NAMES[resampling],
        nodata=np.nan,
    )
    if "band" in out.dims and out.sizes["band"] == 1:
        out = out.squeeze("band", drop=True)
    return out


def cache_processed(
    da: xr.DataArray,
    name: str,
    processed_dir: Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write a reprojected raster as a Cloud-Optimized GeoTIFF.

    Output filename: ``{name}_southpole_240m.tif``. The COG is tiled
    512×512 with DEFLATE compression and internal overviews; idempotent
    unless ``overwrite=True``.

    Args:
        da: DataArray with a CRS set on its ``rio`` accessor.
        name: Logical dataset name (used in the filename).
        processed_dir: Output directory. Created if missing.
        overwrite: When ``False`` (default), skip if the output exists.

    Returns:
        Path to the written COG.

    Raises:
        ValueError: If ``da`` has no CRS.
    """
    if da.rio.crs is None:
        raise ValueError("cache_processed: input DataArray has no CRS")
    processed_dir = Path(processed_dir)
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path = processed_dir / f"{name}_southpole_240m.tif"
    if out_path.exists() and not overwrite:
        return out_path

    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)

    da.rio.to_raster(
        out_path,
        driver="COG",
        compress="DEFLATE",
        BLOCKSIZE=512,
        OVERVIEWS="AUTO",
        RESAMPLING="BILINEAR",
        LEVEL=6,
        BIGTIFF="IF_SAFER",
    )
    return out_path


def is_cog(path: Path) -> bool:
    """Cheap structural check that ``path`` looks like a COG.

    True iff the file opens as GeoTIFF/COG, is internally tiled with
    block sizes ≥ 256, and (for files larger than one block) carries at
    least one overview level.
    """
    with rasterio.open(path) as ds:
        if ds.driver not in {"COG", "GTiff"}:
            return False
        block_h, block_w = ds.block_shapes[0]
        if block_h < 256 or block_w < 256:
            return False
        if ds.width > block_w or ds.height > block_h:
            return bool(ds.overviews(1))
    return True
