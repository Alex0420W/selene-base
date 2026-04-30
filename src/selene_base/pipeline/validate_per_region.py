"""CLI driver for ``selene validate-per-region`` (week 11)."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import typer

from selene_base.validation.comparison import (
    PerRegionComplianceResult,
    per_region_compliance_analysis,
    render_per_region_compliance_summary,
)
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_PER_REGION_SUBDIR = "per_region"


def run(
    *,
    sites_path: Path | None = None,
    summary_path: Path | None = None,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    per_region_subdir: str = DEFAULT_PER_REGION_SUBDIR,
    echo: Callable[[str], None] = typer.echo,
) -> PerRegionComplianceResult:
    """Compute and persist per-region compliance summary.

    Loads ``data/outputs/per_region/sites.geojson`` (or another path
    via ``sites_path``) and the cached per-region summary JSON written
    by ``selene rank-per-region``. Computes the
    :class:`~selene_base.validation.comparison.PerRegionComplianceResult`
    and writes it to ``per_region_validation.json`` alongside the
    legacy ``validation.json``.
    """
    outputs_dir = Path(outputs_dir)
    per_region_dir = outputs_dir / per_region_subdir
    if sites_path is None:
        sites_path = per_region_dir / "sites.geojson"
    if summary_path is None:
        summary_path = per_region_dir / "per_region_summary.json"
    if not sites_path.exists():
        raise FileNotFoundError(
            f"per-region sites not found at {sites_path}; run `selene rank-per-region` first."
        )

    sites = gpd.read_file(sites_path)
    polygons = regions_polygons_to_geodataframe()

    eligible_area_km2: dict[str, float] = {}
    polygon_area_km2: dict[str, float] = {}
    if summary_path.exists():
        summary_blob = json.loads(summary_path.read_text(encoding="utf-8"))
        for entry in summary_blob.get("per_region", []):
            name = str(entry["name"])
            eligible_area_km2[name] = float(entry.get("eligible_area_km2", 0.0))
            poly_cells = entry.get("polygon_cell_area_km2") or entry.get("polygon_area_km2")
            if poly_cells:
                polygon_area_km2[name] = float(poly_cells)

    echo(f"[validate-per-region] {len(sites)} sites across {len(polygons)} USGS regions")
    result = per_region_compliance_analysis(
        sites,
        polygons,
        eligible_area_km2=eligible_area_km2 or None,
        polygon_area_km2=polygon_area_km2 or None,
    )

    out_path = outputs_dir / "per_region_validation.json"
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    echo(f"[done] per_region_validation.json -> {out_path}")
    echo("")
    echo(render_per_region_compliance_summary(result))
    return result
