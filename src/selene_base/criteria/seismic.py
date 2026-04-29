"""Seismic criterion — penalises proximity to active lobate scarps.

Civilini et al. (2023) re-located several Apollo-era shallow moonquakes
to within tens of kilometres of young lobate scarps from the Watters
catalog, suggesting active thrust faulting. The score decreases as a
site approaches any scarp in the catalog.

Filled in week 3.
"""

from __future__ import annotations

import xarray as xr


def compute(grid: xr.DataArray, **kwargs: object) -> xr.DataArray:
    """Score every cell on distance from known active lobate scarps.

    Args:
        grid: DataArray on the common south-polar grid; values are
            unused — only its coordinates are read.
        **kwargs: Tuning knobs. Recognised keys:

            * ``catalog`` (geopandas.GeoDataFrame) — scarp catalog,
              required.
            * ``threshold_km`` (float, default 50.0) — distance at and
              beyond which the score is 1.

    Returns:
        DataArray of [0, 1] scores aligned with ``grid``.

    Raises:
        NotImplementedError: Implementation is filled in week 3.
    """
    raise NotImplementedError("filled in week 3")
