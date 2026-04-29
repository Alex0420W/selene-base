"""Illumination criterion — rewards near-permanent sunlight.

Higher average illumination gives more solar power and longer comms
windows. Built from the Mazarico average-illumination raster (fraction
of the 18.6-year cycle a pixel sees the Sun); diminishing returns past
``target_pct`` because power and thermal systems are sized for the
worst-case duty cycle anyway.

Filled in week 3.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr


def compute(
    illumination: xr.DataArray,
    *,
    target_pct: float = 0.70,
) -> xr.DataArray:
    """Map illumination fraction to a [0, 1] score.

    ``score = clip(illumination / target_pct, 0, 1)``: a site that is
    sunlit ``target_pct`` of the time scores 1.0 and any "extra"
    illumination beyond that does not help (system sizing is locked to
    the worst case).

    Args:
        illumination: DataArray of illumination fraction in [0, 1] on
            the common 240 m grid (e.g. from ``load_illumination``).
        target_pct: Fraction at and above which the score saturates at
            1.0. Must be strictly positive.

    Returns:
        DataArray of [0, 1] scores aligned with ``illumination``;
        NaN where the input is NaN.

    Raises:
        ValueError: If ``target_pct`` is non-positive.
    """
    if target_pct <= 0:
        raise ValueError(f"target_pct must be positive, got {target_pct!r}")

    arr = illumination.to_numpy().astype(np.float64)
    raw = arr / target_pct
    score = np.clip(raw, 0.0, 1.0)

    out = xr.DataArray(
        score,
        coords=illumination.coords,
        dims=illumination.dims,
        name="illumination_score",
    )
    if illumination.rio.crs is not None:
        out = out.rio.write_crs(illumination.rio.crs, inplace=False)
    return out
