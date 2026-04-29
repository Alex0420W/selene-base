"""Extract geographically-distinct top sites from an aggregate score map.

Implements non-maximum suppression in projected metres: pick the highest
remaining cell, emit it, blank a circular ``min_distance_m`` neighbourhood,
repeat until ``n`` sites are collected or the map is exhausted.

Filled in week 3.
"""

from __future__ import annotations

import geopandas as gpd
import xarray as xr


def top_n_sites(
    score_map: xr.DataArray,
    n: int,
    min_distance_m: float,
) -> gpd.GeoDataFrame:
    """Return the top-N sites separated by at least ``min_distance_m``.

    Operates entirely in the projected coordinates of ``score_map``, so
    distance is straightforward Euclidean metres.

    Args:
        score_map: Aggregated [0, 1] suitability scores on the common
            south-polar grid.
        n: Maximum number of sites to return; must be positive.
        min_distance_m: Minimum pairwise separation between returned
            sites, in metres; must be positive.

    Returns:
        GeoDataFrame with one row per site, sorted by descending score,
        containing at least:

        * ``geometry`` — point geometry in the score map's CRS
        * ``score`` — aggregate [0, 1] suitability
        * ``rank`` — 1-based rank

    Raises:
        NotImplementedError: Implementation is filled in week 3.
    """
    raise NotImplementedError("filled in week 3")
