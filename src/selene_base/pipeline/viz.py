"""Generate the standalone web map and per-site HTML reports.

Loads ``data/outputs/score_southpole.tif`` plus ``top_sites.geojson``,
calls :func:`selene_base.viz.webmap.build_map` and
:func:`selene_base.viz.site_report.generate_site_report` (and the
matching ``index`` page), and writes everything as self-contained HTML
under ``data/outputs/`` so the artefacts are shareable with no server
or external CDN dependency.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import typer

from selene_base.validation.nasa_regions import regions_to_geodataframe
from selene_base.viz.site_report import generate_site_index, generate_site_report
from selene_base.viz.webmap import build_map

DEFAULT_OUTPUTS_DIR = Path("data/outputs")


def run(
    *,
    score_map_path: Path | None = None,
    sites_path: Path | None = None,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    processed_dir: Path | None = None,
    skip_reports: bool = False,
    echo: Callable[[str], None] = typer.echo,
) -> dict[str, Path]:
    """Build the web map and per-site reports.

    Args:
        score_map_path: Aggregate score COG path; defaults to
            ``<outputs_dir>/score_southpole.tif``.
        sites_path: Top sites GeoJSON path; defaults to
            ``<outputs_dir>/top_sites.geojson``.
        outputs_dir: Where the artefacts land.
        processed_dir: Where per-criterion COGs live (used for richer
            web-map layers); defaults to ``data/processed``.
        skip_reports: When True, build only the web map.
        echo: Logging sink.

    Returns:
        Mapping with at least ``"webmap"`` (always) and
        ``"sites_index"`` (when reports were generated).
    """
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    if score_map_path is None:
        score_map_path = outputs_dir / "score_southpole.tif"
    if sites_path is None:
        sites_path = outputs_dir / "top_sites.geojson"
    if processed_dir is None:
        processed_dir = Path("data/processed")

    if not score_map_path.exists():
        raise FileNotFoundError(
            f"score map not found at {score_map_path}; run `selene score` first"
        )
    if not sites_path.exists():
        raise FileNotFoundError(f"top sites not found at {sites_path}; run `selene rank` first")

    sites = gpd.read_file(sites_path)
    nasa = regions_to_geodataframe()
    artefacts: dict[str, Path] = {}

    webmap_path = outputs_dir / "webmap.html"
    echo(f"[viz] building web map -> {webmap_path}")
    build_map(
        score_cog=score_map_path,
        top_sites=sites,
        nasa_regions=nasa,
        output_path=webmap_path,
        processed_dir=Path(processed_dir),
    )
    artefacts["webmap"] = webmap_path

    if not skip_reports:
        sites_dir = outputs_dir / "sites"
        sites_dir.mkdir(parents=True, exist_ok=True)
        echo(f"[viz] writing per-site reports -> {sites_dir}")
        report_paths: list[Path] = []
        for _, row in sites.iterrows():
            path = generate_site_report(
                row,
                score_cog=score_map_path,
                output_dir=sites_dir,
                nasa_regions=nasa,
            )
            report_paths.append(path)
        index_path = generate_site_index(sites, sites_dir)
        artefacts["sites_index"] = index_path
        echo(f"[done] {len(report_paths)} site report(s) + index -> {index_path}")

    echo(f"[done] web map -> {webmap_path}")
    return artefacts
