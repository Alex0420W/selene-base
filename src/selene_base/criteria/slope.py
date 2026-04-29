"""Slope criterion — penalises steep terrain.

Computes terrain slope in degrees from the LOLA elevation grid and maps
it to a [0, 1] safety score: roughly 1 on flat ground, 0 above the
``max_slope_deg`` threshold (default 10°, matching Artemis lander
trafficability guidance).

Filled in week 2 (first criterion to land alongside reproject).
"""

from __future__ import annotations

import xarray as xr


def compute(grid: xr.DataArray, **kwargs: object) -> xr.DataArray:
    """Score every cell on slope-derived trafficability.

    Args:
        grid: DataArray on the common south-polar grid containing
            elevation in metres.
        **kwargs: Tuning knobs. Recognised keys:

            * ``max_slope_deg`` (float, default 10.0) — slope at or above
              which the score is zero.

    Returns:
        DataArray of [0, 1] scores aligned with ``grid``; NaN where
        elevation is missing.

    Raises:
        NotImplementedError: Implementation is filled in week 2.
    """
    raise NotImplementedError("filled in week 2")
