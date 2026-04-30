"""CLI driver for ``selene compare-wueller`` (week 12 / v1.4.0).

Loads selene-base's per-region sites and the bundled Wueller 2026
catalog (currently a synthetic placeholder; see
:mod:`selene_base.validation.wueller_comparison` for the data
acquisition status), runs the pairwise comparison, prints a three-block
summary, and writes the full result + per-site distance table.
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
    WUELLER_SITES_CSV,
    WuellerComparisonResult,
    compare_sites,
    load_wueller_sites,
    render_summary,
)

DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_PER_REGION_SITES = DEFAULT_OUTPUTS_DIR / "per_region" / "sites.geojson"


def run(
    *,
    selene_sites_path: Path | None = None,
    wueller_csv: Path | None = None,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    match_threshold_km: float = DEFAULT_MATCH_THRESHOLD_KM,
    echo: Callable[[str], None] = typer.echo,
) -> WuellerComparisonResult:
    """Run the Wueller-comparison harness and persist artefacts."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    if selene_sites_path is None:
        selene_sites_path = DEFAULT_PER_REGION_SITES
    if wueller_csv is None:
        wueller_csv = WUELLER_SITES_CSV
    if not selene_sites_path.exists():
        raise FileNotFoundError(
            f"selene per-region sites not found at {selene_sites_path}; "
            "run `selene rank-per-region` first."
        )

    selene_sites = gpd.read_file(selene_sites_path)
    wueller_sites = load_wueller_sites(wueller_csv)

    echo(
        f"[compare-wueller] {len(selene_sites)} selene sites vs "
        f"{len(wueller_sites)} Wueller sites; threshold = {match_threshold_km:.1f} km"
    )

    result = compare_sites(
        selene_sites,
        wueller_sites,
        match_threshold_km=match_threshold_km,
    )

    json_path = outputs_dir / "wueller_comparison.json"
    csv_path = outputs_dir / "wueller_comparison.csv"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    pd.DataFrame(result["per_selene_site"]).to_csv(csv_path, index=False)

    echo(f"[done] wueller_comparison.json -> {json_path}")
    echo(f"[done] wueller_comparison.csv -> {csv_path}")
    echo("")
    echo(render_summary(result))
    return result
