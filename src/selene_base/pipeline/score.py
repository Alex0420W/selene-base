"""Compute every available criterion score map and aggregate.

Pipeline shape:

1. Load weights YAML.
2. For each criterion that has the inputs it needs:
   a. derive its [0, 1] score grid,
   b. cache it as ``data/processed/scored/<name>_score.tif``.
3. Aggregate the cached score grids via
   :func:`selene_base.scoring.aggregate.weighted_sum`, which warns and
   renormalises across the criteria actually present.
4. Write the final score COG to ``data/outputs/score_southpole.tif``.

Week 2 ships ``slope`` only; week 3 lights up the remaining five.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rioxarray  # noqa: F401
import typer
import xarray as xr
import yaml

from selene_base.criteria import slope as slope_criterion
from selene_base.data.reproject import cache_processed
from selene_base.pipeline.preprocess import load_region_config
from selene_base.scoring.aggregate import weighted_sum

DEFAULT_REGION_CONFIG = Path("config/region_southpole.yaml")
DEFAULT_WEIGHTS = Path("config/weights_default.yaml")
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
SCORED_SUBDIR = "scored"


@dataclass
class ScoreSummary:
    """Per-cell statistics of the final aggregate score."""

    n_finite: int
    minimum: float
    maximum: float
    mean: float
    pct_above_0_7: float
    pct_nodata: float
    output_path: Path

    def render(self) -> str:
        return (
            f"score map: {self.output_path}\n"
            f"  finite cells: {self.n_finite:,}\n"
            f"  min / mean / max: {self.minimum:.3f} / {self.mean:.3f} / {self.maximum:.3f}\n"
            f"  cells with score > 0.7: {self.pct_above_0_7:.1f}%\n"
            f"  cells with no data:     {self.pct_nodata:.1f}%"
        )


def _load_weights(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"weights file must be a mapping, got {type(raw).__name__}")
    return {str(k): float(v) for k, v in raw.items()}


def _compute_slope_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
) -> Path | None:
    """If the LOLA COG is cached, derive slope and write its score COG."""
    lola_cog = processed_dir / "lola_southpole_240m.tif"
    if not lola_cog.exists():
        echo(f"[skip] slope: LOLA COG not present at {lola_cog}; run `selene preprocess`")
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "slope_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] slope: {out_cog.name} already cached")
        return out_cog

    echo(f"[compute] slope from {lola_cog.name}")
    elevation = (
        rioxarray.open_rasterio(lola_cog, masked=True)
        .squeeze("band", drop=True)
        .rename("elevation_m")
    )
    slope_deg = slope_criterion.derive_slope_degrees(elevation, pixel_size_m=pixel_size_m)
    score = slope_criterion.compute(slope_deg)
    out_path = cache_processed(score, "slope_score", scored_dir, overwrite=overwrite)
    echo(f"[done] slope -> {out_path}")
    return out_path


CRITERION_FUNCS: dict[str, Callable[..., Path | None]] = {
    "slope": _compute_slope_score,
}


def run(
    *,
    weights_path: Path = DEFAULT_WEIGHTS,
    region_config: Path = DEFAULT_REGION_CONFIG,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    overwrite: bool = False,
    echo: Callable[[str], None] = typer.echo,
) -> ScoreSummary:
    """Compute all available criterion scores and aggregate.

    Args:
        weights_path: YAML mapping criterion name to non-negative weight.
        region_config: South-polar grid YAML (used for pixel size).
        processed_dir: Where ``preprocess`` wrote cached input COGs and
            where this function writes per-criterion score COGs (under
            ``scored/``).
        outputs_dir: Where the aggregate ``score_southpole.tif`` lands.
        overwrite: Re-compute even when cached score COGs exist.
        echo: Logging sink, defaults to ``typer.echo``.

    Returns:
        :class:`ScoreSummary` with the aggregate's cell stats.
    """
    weights = _load_weights(weights_path)
    cfg = load_region_config(region_config)
    pixel_size_m = float(cfg["resolution_m"])
    processed_dir = Path(processed_dir)
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    score_paths: dict[str, Path] = {}
    for name, fn in CRITERION_FUNCS.items():
        path = fn(
            processed_dir,
            pixel_size_m=pixel_size_m,
            overwrite=overwrite,
            echo=echo,
        )
        if path is not None:
            score_paths[name] = path

    if not score_paths:
        raise RuntimeError(
            "no criterion score grids could be computed; run `selene preprocess` "
            "to populate data/processed/ first."
        )

    score_arrays: dict[str, xr.DataArray] = {}
    for name, path in score_paths.items():
        score_arrays[name] = rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)

    aggregate = weighted_sum(score_arrays, weights)

    out_path = outputs_dir / "score_southpole.tif"
    aggregate.rio.write_nodata(np.nan, inplace=True)
    aggregate = aggregate.rename("aggregate_score")
    cache_processed(aggregate, "score", outputs_dir, overwrite=True)
    actual_path = outputs_dir / "score_southpole_240m.tif"
    if actual_path.exists():
        actual_path.replace(out_path)

    arr = aggregate.to_numpy()
    finite_mask = np.isfinite(arr)
    finite = arr[finite_mask]
    n_finite = int(finite.size)
    if n_finite == 0:
        raise RuntimeError("aggregate score map has no finite values")

    summary = ScoreSummary(
        n_finite=n_finite,
        minimum=float(finite.min()),
        maximum=float(finite.max()),
        mean=float(finite.mean()),
        pct_above_0_7=float(100.0 * (finite > 0.7).sum() / n_finite),
        pct_nodata=float(100.0 * (~finite_mask).sum() / arr.size),
        output_path=out_path,
    )
    echo(summary.render())
    return summary
