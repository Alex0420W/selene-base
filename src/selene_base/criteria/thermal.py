"""Thermal criterion — rewards habitat-friendly mean temperatures.

Operates on the Diviner Polar Resource Product's annual-mean surface
temperature at 2 cm depth (``temp_avg``). The earlier two-input form
that took ``temp_max`` and ``temp_min`` is gone — the PRP gives Tavg
and Tmax, not Tmin, and what matters for habitat thermal control is
the long-term average energy budget the regolith presents to a base
rather than the day/night swing (which is enormous everywhere).

Score is a Gaussian on the mean:

    score = exp(-(Tavg - target)^2 / (2 sigma^2))

with default ``target_temp_k = 140`` and ``sigma_k = 30``.

**Week 8 parameter correction.** Earlier defaults (``target_temp_k = 230``,
``sigma_k = 50``) put the Gaussian's peak at 230 K, *outside* the support
of the data — Diviner PRP ``temp_avg`` for the south pole peaks at 211 K
with a median near 131 K. Every cell scored in the Gaussian's tail and
the criterion contributed almost no discriminative signal (week 7
diagnostic: 0.325 ± 0.198 at our top-20 vs 0.113 ± 0.149 at NASA centroids,
both deep in the tail of a 230 K peak).

The corrected defaults put the peak inside the data distribution
(140 K, near the median) with a narrower width (30 K) so a 30 K offset
from a 140 K target now scores ~0.6 — the kind of contrast the
criterion is designed to provide. This is a parameter correction, not
a tuning preference: the previous values were factually outside the
support of the input.

Filled in week 3; rewritten week 6 for the PRP single-input signature;
defaults corrected week 8.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr


def compute(
    temp_avg: xr.DataArray,
    *,
    target_temp_k: float = 140.0,
    sigma_k: float = 30.0,
) -> xr.DataArray:
    """Map annual-mean surface temperature to a [0, 1] thermal score.

    Args:
        temp_avg: DataArray of annual-mean surface temperature (K).
        target_temp_k: Mean temperature receiving the maximum score.
            Default 140 K, between the data median (131 K) and the
            actual habitat-relevant temperatures (regolith mean energy
            budget around 150–200 K). The previous default of 230 K was
            *outside* the data support and put the Gaussian in its tail
            at every cell; see the module docstring for the full
            correction rationale.
        sigma_k: Gaussian width on the mean temperature, in kelvin.
            Default 30 K — tight enough to discriminate ±30 K offsets
            cleanly (a 30 K offset scores ~0.61).

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
