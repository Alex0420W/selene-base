"""Rasterise vector inputs onto the common 240 m grid.

Currently exports one function — :func:`rasterize_crater_density` — used
by the hazard criterion to convert the Robbins crater catalog into a
per-pixel "how many craters live within R km of this pixel?" density
map. The implementation uses a KDTree on projected crater centres and
a chunked pixel query so 31k × 6.4M lookups stay tractable.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import xarray as xr
from pyproj import Transformer
from scipy.spatial import cKDTree

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS


def rasterize_crater_density(
    craters: gpd.GeoDataFrame,
    target_grid: xr.DataArray,
    *,
    radius_km: float = 3.0,
    diameter_col: str | None = None,
    chunk_rows: int = 256,
) -> xr.DataArray:
    """Count craters within ``radius_km`` of each target-grid pixel.

    Crater centres are reprojected from their native lon/lat CRS into
    the target grid's projected CRS, indexed in a 2-D KDTree, then
    queried in row-chunks against the pixel coordinate grid. The output
    has the same shape and georeferencing as ``target_grid``.

    Args:
        craters: GeoDataFrame with point geometries in
            :data:`selene_base.data.load.LUNAR_GEOGRAPHIC_CRS` (or any
            CRS pyproj recognises) and at minimum a geometry column.
        target_grid: DataArray on the common projected grid; values are
            unused (only ``rio.transform()`` and shape are read).
        radius_km: Search radius in kilometres around each pixel
            centre. Must be strictly positive.
        diameter_col: Optional column to use for filtering craters
            below 1 km diameter — pass ``"DIAM_CIRC_IMG_KM"`` /
            ``"diam_km"`` etc. When ``None`` (default), every crater
            in the input contributes one count.
        chunk_rows: Number of grid rows to query at a time. Lower
            values use less memory; the default fits comfortably on a
            laptop for the 2533² south polar grid.

    Returns:
        DataArray of integer counts (stored as float so NaN can mark
        out-of-region pixels later) with the same dims, coords, and
        CRS as ``target_grid``.

    Raises:
        ValueError: If the target grid has no CRS or transform, if
            ``radius_km`` is non-positive, or if ``craters`` has no
            CRS set.
    """
    if radius_km <= 0:
        raise ValueError(f"radius_km must be positive, got {radius_km!r}")
    if target_grid.rio.crs is None:
        raise ValueError("target_grid has no CRS; set one before rasterising")
    if craters.crs is None:
        raise ValueError("craters has no CRS; set one before rasterising")
    if "y" not in target_grid.dims or "x" not in target_grid.dims:
        raise ValueError(f"target_grid must have ('y', 'x') dims, got {target_grid.dims!r}")

    radius_m = radius_km * 1000.0

    if diameter_col is not None and diameter_col in craters.columns:
        craters = craters[craters[diameter_col].fillna(0) >= 1.0]

    if len(craters) == 0:
        zeros = np.zeros((target_grid.sizes["y"], target_grid.sizes["x"]), dtype=np.float32)
        out = xr.DataArray(
            zeros,
            coords=target_grid.transpose("y", "x").coords,
            dims=("y", "x"),
            name="crater_density",
        )
        return out.rio.write_crs(target_grid.rio.crs, inplace=False)

    # Project crater centres into the target grid's CRS.
    transformer = Transformer.from_crs(craters.crs, target_grid.rio.crs, always_xy=True)
    crater_x, crater_y = transformer.transform(
        craters.geometry.x.to_numpy(), craters.geometry.y.to_numpy()
    )
    tree = cKDTree(np.column_stack([crater_x, crater_y]))

    # Pixel-centre coordinates from the target grid's transform.
    transform = target_grid.rio.transform()
    height, width = target_grid.sizes["y"], target_grid.sizes["x"]
    cols = np.arange(width)
    pixel_x_row = transform.c + transform.a * (cols + 0.5)

    counts = np.zeros((height, width), dtype=np.float32)
    for r0 in range(0, height, chunk_rows):
        r1 = min(r0 + chunk_rows, height)
        rows = np.arange(r0, r1)
        pixel_y = transform.f + transform.e * (rows + 0.5)
        # Cartesian product of (x_row, y_chunk) → (chunk_h * width, 2)
        xs = np.tile(pixel_x_row, len(rows))
        ys = np.repeat(pixel_y, width)
        pts = np.column_stack([xs, ys])
        # query_ball_point returns a list-of-lists of indices per point.
        nbrs = tree.query_ball_point(pts, r=radius_m, return_length=True)
        counts[r0:r1, :] = np.asarray(nbrs, dtype=np.float32).reshape(r1 - r0, width)

    out = xr.DataArray(
        counts,
        coords=target_grid.transpose("y", "x").coords,
        dims=("y", "x"),
        name="crater_density",
    )
    return out.rio.write_crs(target_grid.rio.crs, inplace=False)


__all__ = ["LUNAR_GEOGRAPHIC_CRS", "rasterize_crater_density"]
