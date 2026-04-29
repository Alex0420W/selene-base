"""Thermal criterion — rewards moderate, stable surface temperatures.

Sites that swing between cryogenic shadow and full sunlight stress
hardware; the ideal is a stable, mild regime around the
``target_temp_k`` mean (default 180 K) with a small diurnal range.
Score factors a Gaussian on mean temperature with a linear penalty on
the diurnal swing:

    score = exp(-(T_mean - target)² / (2σ²))
            * max(0, 1 - (T_max - T_min) / max_range)

Filled in week 3.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr


def compute(
    t_max: xr.DataArray,
    t_min: xr.DataArray,
    *,
    target_temp_k: float = 180.0,
    sigma_k: float = 50.0,
    max_range_k: float = 200.0,
) -> xr.DataArray:
    """Map Diviner Tmax/Tmin to a [0, 1] thermal-stability score.

    Args:
        t_max: DataArray of annual maximum bolometric temperature (K).
        t_min: DataArray of annual minimum bolometric temperature (K).
            Must broadcast against ``t_max``.
        target_temp_k: Mean temperature receiving the maximum mean
            score. Must be strictly positive.
        sigma_k: Gaussian width on mean temperature (K). Must be
            strictly positive.
        max_range_k: Diurnal range (K) above which the swing penalty
            saturates the score at zero. Must be strictly positive.

    Returns:
        DataArray of [0, 1] scores; NaN where either input is NaN.

    Raises:
        ValueError: If any of ``target_temp_k``, ``sigma_k``,
            ``max_range_k`` is non-positive, or if shapes mismatch.
    """
    if target_temp_k <= 0:
        raise ValueError(f"target_temp_k must be positive, got {target_temp_k!r}")
    if sigma_k <= 0:
        raise ValueError(f"sigma_k must be positive, got {sigma_k!r}")
    if max_range_k <= 0:
        raise ValueError(f"max_range_k must be positive, got {max_range_k!r}")

    tmax_arr = t_max.to_numpy().astype(np.float64)
    tmin_arr = t_min.to_numpy().astype(np.float64)
    if tmax_arr.shape != tmin_arr.shape:
        raise ValueError(f"t_max shape {tmax_arr.shape!r} != t_min shape {tmin_arr.shape!r}")

    t_mean = 0.5 * (tmax_arr + tmin_arr)
    t_range = tmax_arr - tmin_arr

    mean_score = np.exp(-((t_mean - target_temp_k) ** 2) / (2.0 * sigma_k**2))
    range_score = np.clip(1.0 - (t_range / max_range_k), 0.0, 1.0)
    score = mean_score * range_score

    nan_mask = np.isnan(tmax_arr) | np.isnan(tmin_arr)
    score = np.where(nan_mask, np.nan, score)

    out = xr.DataArray(
        score,
        coords=t_max.coords,
        dims=t_max.dims,
        name="thermal_score",
    )
    if t_max.rio.crs is not None:
        out = out.rio.write_crs(t_max.rio.crs, inplace=False)
    return out
