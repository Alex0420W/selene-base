"""Spatial-coupling criterion — scores cells by joint proximity.

NASA's Artemis III candidate regions cluster along the lunar south
polar rim system, where each candidate is simultaneously a few
kilometres from a permanently-shadowed water-ice deposit *and* a few
kilometres from a sustained-illumination ridge. A linear-sum aggregator
treats those two conditions independently and lets a cell that wins
either one (a deep PSR floor; a far-side flat plain) outrank a cell
that satisfies *both* (the rim band).

This criterion models NASA's actual selection logic — the conjunction
"near a PSR AND near a sunlit ridge" — by **multiplying** two distance
falloffs rather than summing them. Failing either condition drives the
score to zero; both must be satisfied to score well.

Filled in week 7.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr
from scipy.ndimage import distance_transform_edt


def derive_distance_to_psr(
    illumination: xr.DataArray,
    *,
    psr_threshold: float = 0.001,
    pixel_size_m: float = 240.0,
) -> xr.DataArray:
    """Distance (metres) from each cell to the nearest PSR pixel.

    A PSR is defined as a cell where ``illumination < psr_threshold``
    (the small positive cutoff lets numerical noise pass without
    producing false positives). Distance is computed via
    :func:`scipy.ndimage.distance_transform_edt` with explicit pixel
    sampling so the result is in real metres rather than pixels.

    Args:
        illumination: DataArray of illumination fraction in [0, 1].
        psr_threshold: Fraction below which a pixel counts as a PSR;
            must be in ``(0, 1)``.
        pixel_size_m: Edge length of one pixel in metres; must be
            strictly positive.

    Returns:
        DataArray of metres, same shape as ``illumination``. PSR
        pixels themselves return 0; if no PSR pixels exist the
        returned grid is all-NaN.

    Raises:
        ValueError: On out-of-range parameters.
    """
    if not (0 < psr_threshold < 1):
        raise ValueError(f"psr_threshold must be in (0, 1), got {psr_threshold!r}")
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")

    arr = illumination.to_numpy()
    psr_mask = (arr < psr_threshold) & np.isfinite(arr)
    if not psr_mask.any():
        out = np.full(arr.shape, np.nan, dtype=np.float64)
    else:
        # distance_transform_edt computes distance to the *nearest 0*; we
        # want distance to the nearest True (PSR) pixel, so invert.
        distance_pixels = distance_transform_edt(~psr_mask)
        out = distance_pixels.astype(np.float64) * pixel_size_m

    da = xr.DataArray(
        out,
        coords=illumination.coords,
        dims=illumination.dims,
        name="distance_to_psr_m",
    )
    if illumination.rio.crs is not None:
        da = da.rio.write_crs(illumination.rio.crs, inplace=False)
    return da


def derive_distance_to_sunlit_ridge(
    illumination: xr.DataArray,
    slope_deg: xr.DataArray,
    *,
    illumination_threshold: float = 0.70,
    slope_min_deg: float = 5.0,
    slope_max_deg: float = 25.0,
    pixel_size_m: float = 240.0,
) -> xr.DataArray:
    """Distance (metres) from each cell to the nearest "sunlit ridge".

    A sunlit ridge pixel satisfies all three conditions:

    - ``illumination >= illumination_threshold`` (consistent sun)
    - ``slope_deg >= slope_min_deg`` (it's a real ridge, not a plain)
    - ``slope_deg <= slope_max_deg`` (it's buildable, not a cliff)

    The slope band is the geometry of polar crater rims: gentler than
    cliffs, steeper than inter-crater plains. Plains are excluded so
    that a flat well-illuminated cell isn't mistaken for the ridge band
    NASA's candidates sit on.

    Args:
        illumination: DataArray of illumination fraction in [0, 1].
        slope_deg: DataArray of slope in degrees on the same grid.
        illumination_threshold: Lower bound on illumination fraction.
        slope_min_deg: Lower bound on slope (degrees).
        slope_max_deg: Upper bound on slope (degrees).
        pixel_size_m: Pixel edge length in metres.

    Returns:
        DataArray of metres on the same grid. Ridge pixels return 0;
        if no ridge pixels exist the returned grid is all-NaN.

    Raises:
        ValueError: On out-of-range parameters or a shape mismatch
            between ``illumination`` and ``slope_deg``.
    """
    if not (0 <= illumination_threshold <= 1):
        raise ValueError(
            f"illumination_threshold must be in [0, 1], got {illumination_threshold!r}"
        )
    if slope_min_deg < 0 or slope_max_deg <= slope_min_deg:
        raise ValueError(
            f"slope band must satisfy 0 <= min < max; got "
            f"min={slope_min_deg!r} max={slope_max_deg!r}"
        )
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")

    illum_arr = illumination.to_numpy()
    slope_arr = slope_deg.to_numpy()
    if illum_arr.shape != slope_arr.shape:
        raise ValueError(
            f"illumination shape {illum_arr.shape!r} != slope_deg shape {slope_arr.shape!r}"
        )

    ridge = (
        np.isfinite(illum_arr)
        & np.isfinite(slope_arr)
        & (illum_arr >= illumination_threshold)
        & (slope_arr >= slope_min_deg)
        & (slope_arr <= slope_max_deg)
    )
    if not ridge.any():
        out = np.full(illum_arr.shape, np.nan, dtype=np.float64)
    else:
        distance_pixels = distance_transform_edt(~ridge)
        out = distance_pixels.astype(np.float64) * pixel_size_m

    da = xr.DataArray(
        out,
        coords=illumination.coords,
        dims=illumination.dims,
        name="distance_to_sunlit_ridge_m",
    )
    if illumination.rio.crs is not None:
        da = da.rio.write_crs(illumination.rio.crs, inplace=False)
    return da


def compute(
    distance_to_psr_m: xr.DataArray,
    distance_to_ridge_m: xr.DataArray,
    *,
    coupling_distance_km: float = 5.0,
) -> xr.DataArray:
    """Score the spatial-coupling criterion.

    Per cell:

    .. code-block:: text

        d_km          = coupling_distance_km
        score_psr     = max(0, 1 - distance_to_psr_m   / (d_km * 1000))
        score_ridge   = max(0, 1 - distance_to_ridge_m / (d_km * 1000))
        score         = score_psr * score_ridge

    The **product** is the conjunction. The two falloffs are linear
    distance ramps to zero at ``coupling_distance_km`` and they're
    multiplied together so that failing either drives the score to
    zero — exactly the structural property a linear weighted sum
    cannot encode.

    NaN distances (e.g. when one of the two feature classes is empty
    in the input) propagate to NaN scores.

    Args:
        distance_to_psr_m: Distance grid from
            :func:`derive_distance_to_psr`.
        distance_to_ridge_m: Distance grid from
            :func:`derive_distance_to_sunlit_ridge`.
        coupling_distance_km: Distance at which each falloff hits zero.
            Both inputs use the same threshold; defaults to 5 km
            (rough geometry of polar rim crater corridors).

    Returns:
        DataArray of [0, 1] scores aligned with the inputs.

    Raises:
        ValueError: If ``coupling_distance_km`` is non-positive or the
            two input grids have mismatched shapes.
    """
    if coupling_distance_km <= 0:
        raise ValueError(f"coupling_distance_km must be positive, got {coupling_distance_km!r}")

    psr_arr = distance_to_psr_m.to_numpy().astype(np.float64)
    ridge_arr = distance_to_ridge_m.to_numpy().astype(np.float64)
    if psr_arr.shape != ridge_arr.shape:
        raise ValueError(f"shape mismatch: PSR {psr_arr.shape!r} vs ridge {ridge_arr.shape!r}")

    cap_m = coupling_distance_km * 1000.0
    score_psr = np.clip(1.0 - psr_arr / cap_m, 0.0, 1.0)
    score_ridge = np.clip(1.0 - ridge_arr / cap_m, 0.0, 1.0)
    score = score_psr * score_ridge
    nan_mask = np.isnan(psr_arr) | np.isnan(ridge_arr)
    score = np.where(nan_mask, np.nan, score)

    out = xr.DataArray(
        score,
        coords=distance_to_psr_m.coords,
        dims=distance_to_psr_m.dims,
        name="coupling_score",
    )
    if distance_to_psr_m.rio.crs is not None:
        out = out.rio.write_crs(distance_to_psr_m.rio.crs, inplace=False)
    return out
