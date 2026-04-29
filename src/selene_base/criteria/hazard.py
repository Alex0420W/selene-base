"""Hazard criterion — penalises proximity to fresh craters.

Uses the Robbins crater catalog: a site loses score the closer it sits
to a crater whose diameter exceeds ``min_diameter_km``, modelling impact
ejecta and rim instability. Distances are computed in projected metres.

Filled in week 3.
"""

from __future__ import annotations

import xarray as xr


def compute(grid: xr.DataArray, **kwargs: object) -> xr.DataArray:
    """Score every cell on distance from significant impact craters.

    Args:
        grid: DataArray on the common south-polar grid; values are
            unused — only its coordinates are read.
        **kwargs: Tuning knobs. Recognised keys:

            * ``catalog`` (geopandas.GeoDataFrame) — crater catalog,
              required.
            * ``min_diameter_km`` (float, default 5.0) — craters smaller
              than this are ignored.
            * ``buffer_factor`` (float, default 2.0) — ejecta extent as a
              multiple of crater radius.

    Returns:
        DataArray of [0, 1] scores aligned with ``grid``.

    Raises:
        NotImplementedError: Implementation is filled in week 3.
    """
    raise NotImplementedError("filled in week 3")
