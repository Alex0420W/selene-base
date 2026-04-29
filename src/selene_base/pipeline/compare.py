"""Drive the per-criterion diagnostic comparison from the CLI.

Loads ``data/outputs/top_sites.geojson`` and the cached per-criterion
score COGs, runs :func:`selene_base.validation.diagnostic.per_criterion_comparison`,
prints a clean stdout table, and writes
``data/outputs/comparison.json`` plus ``data/outputs/comparison.png``
(a horizontal bar chart of signed delta per criterion).
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import typer

from selene_base.scoring.ranking import load_sub_scores
from selene_base.validation.diagnostic import per_criterion_comparison, render_summary
from selene_base.validation.nasa_regions import regions_to_geodataframe

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")


def _save_delta_plot(df: pd.DataFrame, out_path: Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sorted_df = df.sort_values("delta")
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    colors = ["#06d6a0" if v >= 0 else "#e63946" for v in sorted_df["delta"]]
    ax.barh(sorted_df.index.astype(str), sorted_df["delta"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("delta = our top-20 mean - NASA-centroid mean")
    ax.set_title("Per-criterion score: where we differ from NASA")
    for i, (_crit, row) in enumerate(sorted_df.iterrows()):
        label = f"{row['delta']:+.2f}"
        offset = 0.01 if row["delta"] >= 0 else -0.01
        ax.text(
            row["delta"] + offset,
            i,
            label,
            va="center",
            ha="left" if row["delta"] >= 0 else "right",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return out_path


def run(
    *,
    sites_path: Path | None = None,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    echo: Callable[[str], None] = typer.echo,
) -> pd.DataFrame:
    """Run the per-criterion comparison and persist artefacts."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    processed_dir = Path(processed_dir)

    if sites_path is None:
        sites_path = outputs_dir / "top_sites.geojson"
    if not sites_path.exists():
        raise FileNotFoundError(f"top sites not found at {sites_path}; run `selene rank` first")

    scored_dir = processed_dir / "scored"
    score_maps = load_sub_scores(scored_dir) if scored_dir.exists() else {}
    if not score_maps:
        raise FileNotFoundError(
            f"no per-criterion score COGs in {scored_dir}; run `selene score` first"
        )

    sites = gpd.read_file(sites_path)
    nasa = regions_to_geodataframe()

    df = per_criterion_comparison(sites, nasa, score_maps)
    json_path = outputs_dir / "comparison.json"
    json_path.write_text(
        json.dumps(
            {
                "criteria": df.reset_index().to_dict(orient="records"),
                "n_top_sites": int(len(sites)),
                "n_nasa_regions": int(len(nasa)),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    echo(f"[done] comparison.json -> {json_path}")

    png_path = _save_delta_plot(df, outputs_dir / "comparison.png")
    echo(f"[done] comparison.png -> {png_path}")

    echo("")
    echo(render_summary(df))
    return df
