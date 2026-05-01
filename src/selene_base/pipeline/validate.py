"""Run proximity analysis between top sites and NASA candidate regions.

Loads ``data/outputs/top_sites.geojson``, computes the headline metrics
via :func:`selene_base.validation.comparison.proximity_analysis`, prints
a summary, and writes the full result dict to
``data/outputs/validation.json``.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import typer

from selene_base.validation.comparison import (
    ProximityResult,
    proximity_analysis,
    render_summary,
)
from selene_base.validation.nasa_regions import (
    regions_polygons_to_geodataframe,
    regions_to_geodataframe,
)

DEFAULT_OUTPUTS_DIR = Path("data/outputs")


def run(
    *,
    sites_path: Path | None = None,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    near_km: float = 25.0,
    echo: Callable[[str], None] = typer.echo,
) -> ProximityResult:
    """Validate ranked sites against the NASA Artemis IV (formerly Artemis III) candidates.

    Args:
        sites_path: Path to ``top_sites.geojson``; defaults to
            ``<outputs_dir>/top_sites.geojson``.
        outputs_dir: Where ``validation.json`` lands.
        near_km: Threshold for the "within X km" headline metric.
        echo: Logging sink.

    Returns:
        The full :class:`ProximityResult` dict.
    """
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    if sites_path is None:
        sites_path = outputs_dir / "top_sites.geojson"
    if not sites_path.exists():
        raise FileNotFoundError(f"top sites not found at {sites_path}; run `selene rank` first")

    sites = gpd.read_file(sites_path)
    nasa = regions_to_geodataframe()
    nasa_polygons = regions_polygons_to_geodataframe()
    echo(
        f"[validate] {len(sites)} top sites vs {len(nasa)} NASA centroid disks "
        f"and {len(nasa_polygons)} USGS polygons"
    )
    result = proximity_analysis(sites, nasa, near_km=near_km, nasa_regions_polygons=nasa_polygons)

    json_path = outputs_dir / "validation.json"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    echo(f"[done] validation.json -> {json_path}")
    echo("")
    echo(render_summary(result))
    return result
