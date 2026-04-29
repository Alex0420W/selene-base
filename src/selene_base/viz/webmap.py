"""Folium-based interactive web map of ranked sites.

Renders the aggregate score raster as a translucent overlay with the top
sites as ranked markers and the NASA Artemis III candidate regions for
visual cross-reference.

Filled in week 4.
"""

from __future__ import annotations

from pathlib import Path

import folium
import geopandas as gpd
import xarray as xr


def build_webmap(
    score_map: xr.DataArray,
    sites: gpd.GeoDataFrame,
    out_path: Path,
) -> folium.Map:
    """Render an interactive folium map and write it to ``out_path``.

    Args:
        score_map: Aggregate [0, 1] suitability map on the common grid.
        sites: Ranked sites from
            :func:`selene_base.scoring.ranking.top_n_sites`.
        out_path: Destination HTML file.

    Returns:
        The folium ``Map`` instance, returned for further customisation
        in notebooks.

    Raises:
        NotImplementedError: Implementation is filled in week 4.
    """
    raise NotImplementedError("filled in week 4")
