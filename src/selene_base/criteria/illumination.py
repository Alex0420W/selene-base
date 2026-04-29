"""Illumination criterion — rewards near-permanent sunlight.

Higher average illumination gives more solar power and longer comms
windows. Built from the Mazarico average-illumination raster (fraction
of a year a pixel sees the Sun); score is roughly the raw fraction,
clipped to [0, 1].

Filled in week 3.
"""

from __future__ import annotations

import xarray as xr


def compute(grid: xr.DataArray, **kwargs: object) -> xr.DataArray:
    """Score every cell on average solar illumination.

    Args:
        grid: DataArray on the common south-polar grid holding the
            average illumination fraction in [0, 1].
        **kwargs: Tuning knobs. Recognised keys:

            * ``min_fraction`` (float, default 0.5) — values below this
              are floored to zero rather than scaled linearly.

    Returns:
        DataArray of [0, 1] scores aligned with ``grid``; NaN where
        illumination data is missing.

    Raises:
        NotImplementedError: Implementation is filled in week 3.
    """
    raise NotImplementedError("filled in week 3")
