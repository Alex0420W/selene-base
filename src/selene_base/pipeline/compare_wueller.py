"""CLI driver for ``selene compare-wueller`` (v1.4.1).

Loads selene-base's per-region sites and the bundled Wueller 2026
catalog (130 real sites from doi:10.5281/zenodo.17084058; CC-BY 4.0),
runs the pairwise comparison restricted to NASA's October 2024
down-selected nine regions by default, prints a summary, and writes
the full result + per-site distance table.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import pandas as pd
import typer

from selene_base.validation.wueller_comparison import (
    DEFAULT_MATCH_THRESHOLD_KM,
    WuellerComparisonResult,
    compare_sites,
    evaluate_per_launch_year,
    load_wueller_sites,
    render_launch_year_summary,
    render_summary,
)

DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_PER_REGION_SITES = DEFAULT_OUTPUTS_DIR / "per_region" / "sites.geojson"
V2_OUTPUTS_SUBDIR = "v2"
LAUNCH_YEAR_JSON_NAME = "wueller_comparison_with_launch_year.json"


def run(
    *,
    selene_sites_path: Path | None = None,
    wueller_source: Path | None = None,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    match_threshold_km: float = DEFAULT_MATCH_THRESHOLD_KM,
    filter_to_usgs_scope: bool = True,
    echo: Callable[[str], None] = typer.echo,
) -> WuellerComparisonResult:
    """Run the Wueller-comparison harness and persist artefacts."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    if selene_sites_path is None:
        selene_sites_path = DEFAULT_PER_REGION_SITES
    if not selene_sites_path.exists():
        raise FileNotFoundError(
            f"selene per-region sites not found at {selene_sites_path}; "
            "run `selene rank-per-region` first."
        )

    selene_sites = gpd.read_file(selene_sites_path)
    wueller_sites = load_wueller_sites(wueller_source)

    echo(
        f"[compare-wueller] {len(selene_sites)} selene sites vs "
        f"{len(wueller_sites)} Wueller sites; threshold = {match_threshold_km:.1f} km; "
        f"scope filter = {'on' if filter_to_usgs_scope else 'off'}"
    )

    result = compare_sites(
        selene_sites,
        wueller_sites,
        match_threshold_km=match_threshold_km,
        filter_to_usgs_scope=filter_to_usgs_scope,
    )

    json_path = outputs_dir / "wueller_comparison.json"
    csv_path = outputs_dir / "wueller_comparison.csv"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    pd.DataFrame(result["per_selene_site"]).to_csv(csv_path, index=False)

    echo(f"[done] wueller_comparison.json -> {json_path}")
    echo(f"[done] wueller_comparison.csv -> {csv_path}")

    # v2.0: per-launch-year HLS cross-check using SunDays25-32 columns.
    try:
        launch_year_result = evaluate_per_launch_year(result, wueller_sites)
    except ValueError as exc:
        echo(f"[skip] launch-year cross-check: {exc}")
    else:
        v2_dir = outputs_dir / V2_OUTPUTS_SUBDIR
        v2_dir.mkdir(parents=True, exist_ok=True)
        ly_path = v2_dir / LAUNCH_YEAR_JSON_NAME
        ly_path.write_text(json.dumps(launch_year_result, indent=2), encoding="utf-8")
        echo(f"[done] {LAUNCH_YEAR_JSON_NAME} -> {ly_path}")
        echo("")
        echo(render_launch_year_summary(launch_year_result))

    echo("")
    echo(render_summary(result))
    return result
