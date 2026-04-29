"""Thermal criterion — rewards habitat-friendly mean temperatures.

Operates on the Diviner Polar Resource Product's annual-mean surface
temperature at 2 cm depth (``temp_avg``). The earlier two-input form
that took ``temp_max`` and ``temp_min`` is gone — the PRP gives Tavg
and Tmax, not Tmin, and what matters for habitat thermal control is
the long-term average energy budget the regolith presents to a base
rather than the day/night swing (which is enormous everywhere).

Score is a Gaussian on the mean:

    score = exp(-(Tavg - target)^2 / (2 sigma^2))

with default ``target_temp_k = 230`` (~-43 °C, near typical lunar polar
averages and inside the engineering range for thermal control), and
``sigma = 50 K``.

Filled in week 3; rewritten week 6 for the PRP single-input signature.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr


def compute(
    temp_avg: xr.DataArray,
    *,
    target_temp_k: float = 230.0,
    sigma_k: float = 50.0,
) -> xr.DataArray:
    """Map annual-mean surface temperature to a [0, 1] thermal score.

    Args:
        temp_avg: DataArray of annual-mean surface temperature (K).
        target_temp_k: Mean temperature receiving the maximum score.
            230 K (~-43 °C) is a defensible target — close to typical
            lunar polar annual means and inside the engineering range
            for thermal control. Was ``180 K`` in the original
            week-3 spec when the criterion took Tmin; that value made
            sense for "instantaneous comfortable temperature" but is
            wrong for an annual-average input.
        sigma_k: Gaussian width on the mean temperature, in kelvin.

    Returns:
        DataArray of [0, 1] scores aligned with ``temp_avg``;
        NaN where the input is NaN.

    Raises:
        ValueError: If ``target_temp_k`` or ``sigma_k`` is non-positive.
    """
    if target_temp_k <= 0:
        raise ValueError(f"target_temp_k must be positive, got {target_temp_k!r}")
    if sigma_k <= 0:
        raise ValueError(f"sigma_k must be positive, got {sigma_k!r}")

    arr = temp_avg.to_numpy().astype(np.float64)
    score = np.exp(-((arr - target_temp_k) ** 2) / (2.0 * sigma_k**2))
    score = np.where(np.isnan(arr), np.nan, score)

    out = xr.DataArray(
        score,
        coords=temp_avg.coords,
        dims=temp_avg.dims,
        name="thermal_score",
    )
    if temp_avg.rio.crs is not None:
        out = out.rio.write_crs(temp_avg.rio.crs, inplace=False)
    return out
