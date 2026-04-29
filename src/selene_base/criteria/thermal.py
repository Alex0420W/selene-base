"""Thermal criterion — rewards moderate, stable surface temperatures.

Sites that swing between cryogenic shadow and full sunlight stress
hardware; the ideal is a stable, mild regime around the ``target_K``
mean (default 230 K). Built on Diviner bolometric temperature
climatology.

Filled in week 3.
"""

from __future__ import annotations

import xarray as xr


def compute(grid: xr.DataArray, **kwargs: object) -> xr.DataArray:
    """Score every cell on thermal stability.

    Args:
        grid: DataArray on the common south-polar grid holding mean
            surface temperature in kelvin.
        **kwargs: Tuning knobs. Recognised keys:

            * ``target_K`` (float, default 230.0) — peak-score temperature.
            * ``sigma_K`` (float, default 30.0) — Gaussian width.

    Returns:
        DataArray of [0, 1] scores aligned with ``grid``; NaN where
        temperature data is missing.

    Raises:
        NotImplementedError: Implementation is filled in week 3.
    """
    raise NotImplementedError("filled in week 3")
