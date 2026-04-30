"""Drive the weight sensitivity sweep from the CLI.

Loads the cached per-criterion score COGs that ``selene score`` wrote
to ``data/processed/scored/``, draws ``n_samples`` weight vectors via
Latin hypercube on the simplex, runs proximity analysis for each, and
emits ``data/outputs/sensitivity_results.parquet`` plus a histogram
PNG that's safe to embed in the README.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import typer
import yaml

from selene_base.scoring.ranking import DEFAULT_CRITERIA, load_sub_scores
from selene_base.validation.nasa_regions import (
    regions_polygons_to_geodataframe,
    regions_to_geodataframe,
)
from selene_base.validation.sensitivity import (
    best_weights,
    latin_hypercube_weights,
    render_summary,
    save_results,
    sweep_weights,
)

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_WEIGHTS = Path("config/weights_default.yaml")


def _load_default_weights(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    return {str(k): float(v) for k, v in raw.items()}


def _save_distribution_plot(
    results: pd.DataFrame,
    out_path: Path,
    *,
    default_match_count: int | None,
    proximity_threshold_km: float,
) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    target_col = "n_regions_within_proximity_km"
    n_regions = 9  # fixed in NASA list, but recover from data when possible.

    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    counts = results[target_col].astype(int)
    bin_edges = list(range(0, int(counts.max()) + 2))
    if len(bin_edges) < 3:
        bin_edges = [0, 1, 2]
    ax.hist(counts, bins=bin_edges, align="left", color="#06d6a0", edgecolor="#04362a")
    ax.set_xlabel(
        f"NASA regions with a top site within {proximity_threshold_km:.0f} km (out of {n_regions})"
    )
    ax.set_ylabel("number of weight samples")
    ax.set_title(f"sensitivity over {len(results)} weight samples")

    if default_match_count is not None:
        ax.axvline(
            default_match_count,
            color="#e63946",
            linestyle="--",
            linewidth=2,
            label=f"default weights: {default_match_count}",
        )
        ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return out_path


def run(
    *,
    n_samples: int = 200,
    top_n: int = 20,
    min_distance_km: float = 25.0,
    proximity_threshold_km: float = 25.0,
    far_threshold_km: float = 100.0,
    seed: int = 42,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    weights_path: Path = DEFAULT_WEIGHTS,
    echo: Callable[[str], None] = typer.echo,
) -> pd.DataFrame:
    """Sweep weight vectors and persist + summarise the results."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    processed_dir = Path(processed_dir)

    scored_dir = processed_dir / "scored"
    score_maps = load_sub_scores(scored_dir) if scored_dir.exists() else {}
    if not score_maps:
        raise FileNotFoundError(
            f"no per-criterion score COGs in {scored_dir}; run `selene score` first"
        )

    criterion_names = [name for name in DEFAULT_CRITERIA if name in score_maps]
    echo(
        f"[sweep] {n_samples} weight samples over {len(criterion_names)} criteria: "
        f"{', '.join(criterion_names)}"
    )

    samples = latin_hypercube_weights(n_samples, criterion_names, seed=seed)
    nasa = regions_to_geodataframe()
    nasa_polygons = regions_polygons_to_geodataframe()
    results = sweep_weights(
        score_maps,
        samples,
        nasa,
        top_n=top_n,
        min_distance_km=min_distance_km,
        proximity_threshold_km=proximity_threshold_km,
        far_threshold_km=far_threshold_km,
        criterion_order=criterion_names,
        nasa_regions_polygons=nasa_polygons,
    )

    parquet_path = outputs_dir / "sensitivity_results.parquet"
    save_results(results, parquet_path)
    echo(f"[done] sensitivity_results.parquet -> {parquet_path}")

    default_weights_full = _load_default_weights(weights_path)
    default_match_count: int | None = None
    if default_weights_full:
        default_subset_total = sum(default_weights_full.get(name, 0.0) for name in criterion_names)
        if default_subset_total > 0:
            default_subset = {
                name: default_weights_full.get(name, 0.0) / default_subset_total
                for name in criterion_names
            }
            # Closest sample to the renormalised default weight vector.
            diff = sum((results[f"w_{n}"] - default_subset[n]).abs() for n in criterion_names)
            nearest_idx = diff.idxmin()
            default_match_count = int(results.loc[nearest_idx, "n_regions_within_proximity_km"])

    plot_path = outputs_dir / "sensitivity_distribution.png"
    _save_distribution_plot(
        results,
        plot_path,
        default_match_count=default_match_count,
        proximity_threshold_km=proximity_threshold_km,
    )
    echo(f"[done] sensitivity_distribution.png -> {plot_path}")
    echo("")
    echo(
        render_summary(
            results,
            proximity_threshold_km=proximity_threshold_km,
        )
    )

    best = best_weights(results, criterion_names)
    best_idx = int(results["n_regions_within_proximity_km"].idxmax())
    best_row = results.iloc[best_idx]
    echo("")
    echo(
        f"  best sample: regions matched within {proximity_threshold_km:.0f} km = "
        f"{int(best_row['n_regions_within_proximity_km'])}"
    )
    echo("    weight vector: " + ", ".join(f"{k}={v:.2f}" for k, v in best.items()))
    return results
