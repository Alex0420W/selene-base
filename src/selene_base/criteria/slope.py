"""Slope criterion — penalises steep terrain.

Two functions:

- :func:`derive_slope_degrees` reads an elevation grid and returns
  per-cell slope in degrees via second-order central differences in
  metric units (Zevenbergen & Thorne 1987 convention).
- :func:`compute` maps slope degrees to a [0, 1] safety score using
  :func:`selene_base.scoring.normalize.inverse_threshold`.

Filled in week 2 (first criterion to land alongside reproject).
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from selene_base.scoring.normalize import inverse_threshold


def derive_slope_degrees(
    elevation: xr.DataArray,
    pixel_size_m: float,
) -> xr.DataArray:
    """Compute slope in degrees from a square-pixel elevation grid.

    Uses :func:`numpy.gradient` with explicit metric spacing to obtain
    ``dz/dy`` and ``dz/dx``; magnitude of the gradient gives the
    steepest local slope, which is then ``arctan``'d into degrees.
    This is the Zevenbergen & Thorne (1987) convention — the most
    common choice in planetary GIS — and gives values within ~5% of
    the Sobel-weighted Horn (1981) kernel on smooth surfaces.

    Edge pixels (first/last row, first/last column) are explicitly
    set to NaN rather than left as one-sided differences, because
    the rest of the pipeline treats NaN as "no signal here".

    Args:
        elevation: 2-D DataArray of elevation in metres on a square-
            pixel grid; must have dims that include ``("y", "x")``.
        pixel_size_m: Edge length of one pixel in metres. Both axes
            are assumed to have the same spacing (true on the south-
            polar stereographic grid).

    Returns:
        DataArray of slope in degrees, same coords as ``elevation``,
        with edges and propagated-NaN cells set to NaN.

    Raises:
        ValueError: If ``pixel_size_m`` is non-positive.
    """
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")
    if "y" not in elevation.dims or "x" not in elevation.dims:
        raise ValueError(f"elevation must have ('y', 'x') dims, got {elevation.dims!r}")

    z = elevation.transpose("y", "x").to_numpy().astype(np.float64)
    dzdy, dzdx = np.gradient(z, pixel_size_m, pixel_size_m)
    slope_rad = np.arctan(np.hypot(dzdx, dzdy))
    slope_deg = np.degrees(slope_rad)

    # Force edges to NaN — np.gradient uses one-sided differences there.
    slope_deg[0, :] = np.nan
    slope_deg[-1, :] = np.nan
    slope_deg[:, 0] = np.nan
    slope_deg[:, -1] = np.nan
    # Propagate input NaNs: a finite-diff stencil that touched any NaN
    # should yield NaN, but np.gradient doesn't enforce that, so be
    # explicit.
    slope_deg = np.where(np.isnan(z), np.nan, slope_deg)

    out = xr.DataArray(
        slope_deg,
        coords=elevation.transpose("y", "x").coords,
        dims=("y", "x"),
        name="slope_deg",
    )
    if elevation.rio.crs is not None:
        out = out.rio.write_crs(elevation.rio.crs, inplace=False)
    return out


def compute(
    slope_deg: xr.DataArray,
    *,
    max_slope_deg: float = 15.0,
) -> xr.DataArray:
    """Map slope degrees to a [0, 1] safety score.

    Linear ramp: 0° → 1.0, ``max_slope_deg`` → 0.0, anything beyond
    ``max_slope_deg`` → 0.0. Implemented via
    :func:`selene_base.scoring.normalize.inverse_threshold`.

    The default 15° threshold is conservative for Artemis-class
    landers (NASA design guidance is ≤10° at the touchdown footprint,
    ≤15° within a few-hundred-metre operational radius).

    Args:
        slope_deg: 2-D DataArray of slope in degrees.
        max_slope_deg: Slope (degrees) at and above which the score is
            zero; must be strictly positive.

    Returns:
        DataArray of [0, 1] scores aligned with ``slope_deg``; NaN where
        ``slope_deg`` is NaN.

    Raises:
        ValueError: If ``max_slope_deg`` is non-positive.
    """
    if max_slope_deg <= 0:
        raise ValueError(f"max_slope_deg must be positive, got {max_slope_deg!r}")

    score_arr = inverse_threshold(slope_deg.to_numpy(), threshold=max_slope_deg)
    out = xr.DataArray(
        score_arr,
        coords=slope_deg.coords,
        dims=slope_deg.dims,
        name="slope_score",
    )
    if slope_deg.rio.crs is not None:
        out = out.rio.write_crs(slope_deg.rio.crs, inplace=False)
    return out
