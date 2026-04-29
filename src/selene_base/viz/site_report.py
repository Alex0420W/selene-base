"""Per-site HTML report renderer.

For each top-N site, emits a standalone HTML page summarising:

- Its rank, aggregate score, and projected coordinates.
- A per-criterion breakdown (score and the value behind it).
- Local context plots (slope, illumination, thermal, ice, hazard).
- Distance to the closest NASA Artemis III candidate region.

Filled in week 4.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import xarray as xr


def render_site_report(
    site: gpd.GeoSeries,
    score_map: xr.DataArray,
    criterion_maps: dict[str, xr.DataArray],
    out_path: Path,
) -> Path:
    """Render the HTML report for a single ranked site.

    Args:
        site: One row from the ranked sites GeoDataFrame.
        score_map: Aggregate [0, 1] suitability map on the common grid.
        criterion_maps: Mapping from criterion name to its [0, 1] score
            grid; used for the per-criterion section.
        out_path: Destination HTML file.

    Returns:
        The path written to.

    Raises:
        NotImplementedError: Implementation is filled in week 4.
    """
    raise NotImplementedError("filled in week 4")
