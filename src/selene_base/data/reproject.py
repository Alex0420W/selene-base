"""Warp rasters onto the common south-polar analysis grid.

The downstream pipeline assumes every score map shares one grid: the
stereographic projection defined in ``config/region_southpole.yaml`` at
240 m / pixel over a 304 km half-extent. This module is the only place
that touches CRS transforms.

Filled in week 2.
"""

from __future__ import annotations

import xarray as xr

Bounds = tuple[float, float, float, float]


def reproject_to_grid(
    src: xr.DataArray,
    target_crs: str,
    bounds: Bounds,
    res: float,
) -> xr.DataArray:
    """Reproject ``src`` onto a regular grid in ``target_crs``.

    Args:
        src: Source DataArray with CRS metadata available on its ``rio``
            accessor.
        target_crs: Target CRS as a PROJ string or EPSG code understood by
            pyproj (e.g. the south-polar stereographic in the region config).
        bounds: ``(xmin, ymin, xmax, ymax)`` extent in target CRS units.
        res: Pixel size in target CRS units (square pixels).

    Returns:
        DataArray on the target grid with ``y`` decreasing and ``x``
        increasing, NaN-filled outside the source footprint.

    Raises:
        NotImplementedError: Implementation is filled in week 2.
    """
    raise NotImplementedError("filled in week 2")
