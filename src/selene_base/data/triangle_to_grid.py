"""Interpolate scalar values from an irregular triangular mesh to a regular grid.

The Diviner Polar Resource Product reports scalar quantities (Tavg,
Tmax, ice depth) at the centres of an icosahedral triangular mesh.
This module reprojects the irregular ``(lon, lat, value)`` cloud into
the project's south-polar stereographic CRS and runs
:func:`scipy.interpolate.griddata` to produce a regular raster on the
analysis grid.

Two methods are exposed:

- ``"linear"`` — Delaunay-triangulated linear interpolation; right
  default for continuous fields like temperature.
- ``"nearest"`` — Voronoi nearest-neighbour. Right for fields that
  are effectively discontinuous (ice-depth crossings, categorical
  masks).
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr
from pyproj import Transformer
from scipy.interpolate import griddata

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

ALLOWED_METHODS = ("linear", "nearest", "cubic")


def triangles_to_raster(
    points_lat: np.ndarray,
    points_lon: np.ndarray,
    values: np.ndarray,
    target_grid: xr.DataArray,
    *,
    method: str = "linear",
    fill_value: float = np.nan,
) -> xr.DataArray:
    """Rasterise a scalar field defined on triangle centres onto a regular grid.

    Args:
        points_lat: Latitudes of the triangle centres (degrees,
            planetocentric).
        points_lon: Longitudes of the triangle centres (degrees,
            positive east).
        values: Scalar value at each centre. NaN entries are dropped
            before interpolation.
        target_grid: DataArray on the projected analysis grid; only its
            CRS, transform, and shape are read.
        method: ``"linear"``, ``"nearest"``, or ``"cubic"``.
        fill_value: Value used outside the input mesh's convex hull.

    Returns:
        DataArray on ``target_grid``'s grid with the same dims and CRS.

    Raises:
        ValueError: On an unknown ``method`` or on shape mismatches
            between ``points_lat``, ``points_lon``, and ``values``.
    """
    if method not in ALLOWED_METHODS:
        raise ValueError(f"unknown method {method!r}; choose from {ALLOWED_METHODS}")
    points_lat = np.asarray(points_lat, dtype=np.float64)
    points_lon = np.asarray(points_lon, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if not (points_lat.shape == points_lon.shape == values.shape):
        raise ValueError(
            "points_lat, points_lon, values must share shape; "
            f"got {points_lat.shape}, {points_lon.shape}, {values.shape}"
        )
    if target_grid.rio.crs is None:
        raise ValueError("target_grid has no CRS")

    finite = np.isfinite(points_lat) & np.isfinite(points_lon) & np.isfinite(values)
    points_lat = points_lat[finite]
    points_lon = points_lon[finite]
    values = values[finite]

    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, target_grid.rio.crs, always_xy=True)
    src_xs, src_ys = transformer.transform(points_lon, points_lat)
    src_xy = np.column_stack([src_xs, src_ys])

    transform = target_grid.rio.transform()
    h, w = target_grid.sizes["y"], target_grid.sizes["x"]
    cols = np.arange(w)
    rows = np.arange(h)
    grid_x = transform.c + transform.a * (cols + 0.5)
    grid_y = transform.f + transform.e * (rows + 0.5)
    mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)

    interpolated = griddata(
        src_xy,
        values,
        (mesh_x, mesh_y),
        method=method,
        fill_value=fill_value,
    )

    out = xr.DataArray(
        interpolated.astype(np.float32),
        dims=("y", "x"),
        coords={"y": grid_y, "x": grid_x},
    )
    return out.rio.write_crs(target_grid.rio.crs, inplace=False)
