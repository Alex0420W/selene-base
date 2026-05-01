"""CLI driver for ``selene score-wueller-sites`` (v1.5 catalog report).

Evaluates selene's six active 240 m criterion rasters and aggregate
score at every in-scope Wueller 2026 site. Pairs with
``selene compare-wueller``: that subcommand answers "are the two
catalogs picking the same cells?", this one answers "do they agree on
what makes a good cell?". Writes a CSV + JSON for downstream report
generation; the v1.5 catalog report at
``docs/v1.5_catalog_report.md`` is the primary consumer.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import typer

from selene_base.validation.wueller_comparison import load_wueller_sites
from selene_base.validation.wueller_scoring import (
    DEFAULT_OUTPUTS_DIR,
    DEFAULT_PROCESSED_DIR,
    score_wueller_sites,
)


def run(
    *,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    csv_filename: str = "v1.5_catalog_wueller_evaluation.csv",
    json_filename: str = "wueller_sites_scored_by_selene.json",
    in_scope_only: bool = True,
    echo: Callable[[str], None] = typer.echo,
) -> Path:
    """Score Wueller sites against selene's criterion rasters and persist."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    wueller_sites = load_wueller_sites()
    df = score_wueller_sites(
        wueller_sites=wueller_sites,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        in_scope_only=in_scope_only,
    )

    csv_path = outputs_dir / csv_filename
    json_path = outputs_dir / json_filename
    df.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(
            {
                "n_sites": int(len(df)),
                "in_scope_only": bool(in_scope_only),
                "hls_compliant_count": int(df["hls_compliant"].sum()),
                "median_aggregate_score": float(df["aggregate_score"].median()),
                "median_aggregate_score_by_region": {
                    str(region): float(group["aggregate_score"].median())
                    for region, group in df.groupby("region", sort=True)
                },
                "rows": df.to_dict(orient="records"),
            },
            indent=2,
            default=lambda o: bool(o) if hasattr(o, "item") else str(o),
        ),
        encoding="utf-8",
    )

    echo(
        f"[score-wueller] {len(df)} sites evaluated; "
        f"hls_compliant={int(df['hls_compliant'].sum())}/{len(df)}; "
        f"median aggregate={float(df['aggregate_score'].median()):.3f}"
    )
    echo(f"[done] csv  -> {csv_path}")
    echo(f"[done] json -> {json_path}")
    return csv_path
