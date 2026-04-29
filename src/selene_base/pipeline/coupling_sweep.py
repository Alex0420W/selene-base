"""Sweep the spatial-coupling distance parameter and report alignment.

Tunes the only free knob the coupling criterion exposes — the
``coupling_distance_km`` cap on the per-cell PSR-distance and
ridge-distance falloffs — and reports how NASA-region alignment shifts
across a small grid of values. Output: a parquet of per-distance
metrics plus a PNG plot of "regions matched within 25 km" against
``coupling_distance_km``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import rioxarray  # noqa: F401
import typer
import yaml

from selene_base.criteria import coupling as coupling_criterion
from selene_base.scoring.ranking import load_sub_scores
from selene_base.validation.nasa_regions import regions_to_geodataframe
from selene_base.validation.sensitivity import sweep_coupling_distance

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_WEIGHTS = Path("config/weights_default.yaml")


def _load_weights(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    return {str(k): float(v) for k, v in raw.items()}


def _save_plot(df: pd.DataFrame, out_path: Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.plot(
        df["coupling_distance_km"],
        df["n_regions_within_proximity_km"],
        marker="o",
        color="#06d6a0",
        linewidth=2,
        label="regions matched within 25 km",
    )
    ax.plot(
        df["coupling_distance_km"],
        df["n_within_proximity_km"],
        marker="s",
        color="#1d4ed8",
        linewidth=1.5,
        alpha=0.7,
        label="top sites within 25 km",
    )
    ax.set_xlabel("coupling_distance_km")
    ax.set_ylabel("count")
    ax.set_title("Spatial-coupling distance sweep vs NASA Artemis III alignment")
    ax.set_xticks(df["coupling_distance_km"])
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return out_path


def run(
    *,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    weights_path: Path = DEFAULT_WEIGHTS,
    distances_km: list[float] | None = None,
    top_n: int = 20,
    min_distance_km: float = 25.0,
    proximity_threshold_km: float = 25.0,
    echo: Callable[[str], None] = typer.echo,
) -> pd.DataFrame:
    """Run the coupling-distance sweep and persist results."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    processed_dir = Path(processed_dir)

    illum_cog = processed_dir / "illumination_southpole_240m.tif"
    slope_deg_cog = processed_dir / "lola_slope_deg_southpole_240m.tif"
    if not (illum_cog.exists() and slope_deg_cog.exists()):
        raise FileNotFoundError(
            f"need {illum_cog.name} and {slope_deg_cog.name}; run `selene preprocess`"
        )

    illum = rioxarray.open_rasterio(illum_cog, masked=True).squeeze("band", drop=True)
    slope_deg = rioxarray.open_rasterio(slope_deg_cog, masked=True).squeeze("band", drop=True)

    echo("[derive] distance-to-PSR + distance-to-sunlit-ridge")
    distance_psr = coupling_criterion.derive_distance_to_psr(illum, pixel_size_m=240.0)
    distance_ridge = coupling_criterion.derive_distance_to_sunlit_ridge(
        illum, slope_deg, pixel_size_m=240.0
    )

    scored_dir = processed_dir / "scored"
    sub_scores = load_sub_scores(scored_dir) if scored_dir.exists() else {}
    score_maps_no_coupling = {name: arr for name, arr in sub_scores.items() if name != "coupling"}
    if not score_maps_no_coupling:
        raise FileNotFoundError(
            f"no per-criterion score COGs in {scored_dir}; run `selene score` first"
        )

    weights = _load_weights(weights_path)
    crit_present = ["coupling", *score_maps_no_coupling.keys()]
    weights = {k: weights.get(k, 0.0) for k in crit_present}
    if "coupling" not in weights or weights.get("coupling", 0) == 0:
        weights["coupling"] = 0.20

    nasa = regions_to_geodataframe()

    if distances_km is None:
        distances_km = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    echo(f"[sweep] coupling_distance_km in {distances_km}")

    df = sweep_coupling_distance(
        score_maps_no_coupling,
        distance_psr,
        distance_ridge,
        weights,
        nasa,
        distances_km=distances_km,
        top_n=top_n,
        min_distance_km=min_distance_km,
        proximity_threshold_km=proximity_threshold_km,
    )

    parquet_path = outputs_dir / "coupling_sweep.parquet"
    df.to_parquet(parquet_path)
    echo(f"[done] coupling_sweep.parquet -> {parquet_path}")

    plot_path = _save_plot(df, outputs_dir / "coupling_distance_sweep.png")
    echo(f"[done] coupling_distance_sweep.png -> {plot_path}")

    echo("")
    cols = [
        "coupling_distance_km",
        "n_within_proximity_km",
        "n_regions_within_proximity_km",
        "n_regions_with_top_site",
    ]
    echo(df[cols].to_string(index=False))
    echo("")
    best = df.loc[df["n_regions_within_proximity_km"].idxmax()]
    echo(
        f"  best coupling_distance_km = {best['coupling_distance_km']:.1f} "
        f"({int(best['n_regions_within_proximity_km'])} / 9 NASA regions matched)"
    )
    return df


__all__ = ["run"]
