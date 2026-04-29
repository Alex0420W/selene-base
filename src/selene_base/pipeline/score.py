"""Compute every available criterion score map and aggregate.

Pipeline shape:

1. Load weights YAML.
2. For each criterion that has the inputs it needs:
   a. derive its [0, 1] score grid,
   b. cache it as ``data/processed/scored/<name>_score_southpole_240m.tif``.
3. Aggregate the cached score grids via
   :func:`selene_base.scoring.aggregate.weighted_sum`, which warns and
   renormalises across the criteria actually present.
4. Write the final score COG to ``data/outputs/score_southpole.tif``.

Week 3 covers all six criteria; missing source data (Diviner, LEND,
scarps) makes the relevant criterion skip cleanly.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401
import typer
import xarray as xr
import yaml

from selene_base.criteria import (
    hazard as hazard_criterion,
)
from selene_base.criteria import (
    ice as ice_criterion,
)
from selene_base.criteria import (
    illumination as illumination_criterion,
)
from selene_base.criteria import (
    seismic as seismic_criterion,
)
from selene_base.criteria import (
    slope as slope_criterion,
)
from selene_base.criteria import (
    thermal as thermal_criterion,
)
from selene_base.data.reproject import cache_processed
from selene_base.pipeline.preprocess import load_region_config
from selene_base.scoring.aggregate import weighted_sum

DEFAULT_REGION_CONFIG = Path("config/region_southpole.yaml")
DEFAULT_WEIGHTS = Path("config/weights_default.yaml")
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_RAW_DIR = Path("data/raw")
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
    criteria_used: tuple[str, ...]
    criteria_skipped: tuple[str, ...]

    def render(self) -> str:
        return (
            f"score map: {self.output_path}\n"
            f"  criteria used:    {', '.join(self.criteria_used) or '(none)'}\n"
            f"  criteria skipped: {', '.join(self.criteria_skipped) or '(none)'}\n"
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


def _open_cog(path: Path) -> xr.DataArray:
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


# --------------------------------------------------------------------------
# Per-criterion helpers. Each returns Path | None — None signals "skip me".
# --------------------------------------------------------------------------
def _slope_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
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
    elevation = _open_cog(lola_cog).rename("elevation_m")
    slope_deg = slope_criterion.derive_slope_degrees(elevation, pixel_size_m=pixel_size_m)
    score = slope_criterion.compute(slope_deg)
    out_path = cache_processed(score, "slope_score", scored_dir, overwrite=overwrite)
    echo(f"[done] slope -> {out_path}")
    return out_path


def _illumination_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    illum_cog = processed_dir / "illumination_southpole_240m.tif"
    if not illum_cog.exists():
        echo(f"[skip] illumination: input COG not present at {illum_cog}")
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "illumination_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] illumination: {out_cog.name} already cached")
        return out_cog

    echo(f"[compute] illumination from {illum_cog.name}")
    illum = _open_cog(illum_cog)
    score = illumination_criterion.compute(illum)
    out_path = cache_processed(score, "illumination_score", scored_dir, overwrite=overwrite)
    echo(f"[done] illumination -> {out_path}")
    return out_path


def _hazard_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    crater_cog = processed_dir / "crater_density_southpole_240m.tif"
    if not crater_cog.exists():
        echo(f"[skip] hazard: crater density COG not present at {crater_cog}")
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "hazard_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] hazard: {out_cog.name} already cached")
        return out_cog

    echo(f"[compute] hazard from {crater_cog.name}")
    density = _open_cog(crater_cog)
    score = hazard_criterion.compute(density)
    out_path = cache_processed(score, "hazard_score", scored_dir, overwrite=overwrite)
    echo(f"[done] hazard -> {out_path}")
    return out_path


def _thermal_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    temp_avg_cog = processed_dir / "diviner_temp_avg_southpole_240m.tif"
    if not temp_avg_cog.exists():
        echo(
            "[skip] thermal: Diviner PRP temp_avg COG not present; "
            "run `selene download diviner` then `selene preprocess`"
        )
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "thermal_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] thermal: {out_cog.name} already cached")
        return out_cog

    echo("[compute] thermal from Diviner PRP temp_avg")
    score = thermal_criterion.compute(_open_cog(temp_avg_cog))
    out_path = cache_processed(score, "thermal_score", scored_dir, overwrite=overwrite)
    echo(f"[done] thermal -> {out_path}")
    return out_path


def _ice_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    ice_depth_cog = processed_dir / "diviner_ice_depth_southpole_240m.tif"
    if not ice_depth_cog.exists():
        echo(
            "[skip] ice: Diviner PRP ice_depth COG not present; "
            "run `selene download diviner` then `selene preprocess`"
        )
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "ice_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] ice: {out_cog.name} already cached")
        return out_cog

    echo(f"[compute] ice from {ice_depth_cog.name}")
    ice_depth = _open_cog(ice_depth_cog)
    psr_mask: xr.DataArray | None = None
    illum_cog = processed_dir / "illumination_southpole_240m.tif"
    if illum_cog.exists():
        echo("           (using PSR mask derived from illumination)")
        psr_mask = ice_criterion.derive_psr_mask(_open_cog(illum_cog))
    score = ice_criterion.compute(ice_depth, psr_mask=psr_mask, pixel_size_m=pixel_size_m)
    out_path = cache_processed(score, "ice_score", scored_dir, overwrite=overwrite)
    echo(f"[done] ice -> {out_path}")
    return out_path


def _seismic_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "seismic_score_southpole_240m.tif"

    candidates = [
        raw_dir / "scarps" / "scarps_southpole.geojson",
        raw_dir / "scarps" / "scarps_southpole.csv",
    ]
    scarps_path = next((p for p in candidates if p.exists()), None)
    if scarps_path is None:
        echo(
            "[skip] seismic: Watters scarp catalog not present "
            "(see download_scarps docstring; place file at "
            f"{candidates[0].as_posix()})"
        )
        return None

    if out_cog.exists() and not overwrite:
        echo(f"[skip] seismic: {out_cog.name} already cached")
        return out_cog

    # Need a target grid in the projected CRS — borrow it from the LOLA COG.
    lola_cog = processed_dir / "lola_southpole_240m.tif"
    if not lola_cog.exists():
        echo(
            "[skip] seismic: needs a target grid; run `selene preprocess` "
            "to produce lola_southpole_240m.tif first"
        )
        return None

    echo(f"[compute] seismic from {scarps_path.name}")
    if scarps_path.suffix.lower() == ".csv":
        import pandas as pd
        from shapely.geometry import Point

        df = pd.read_csv(scarps_path)
        scarps = gpd.GeoDataFrame(
            df,
            geometry=[Point(xy) for xy in zip(df["lon"], df["lat"], strict=True)],
        ).set_crs("+proj=longlat +R=1737400 +no_defs +type=crs")
    else:
        scarps = gpd.read_file(scarps_path)
    target_grid = _open_cog(lola_cog)
    distance = seismic_criterion.distance_to_scarps(scarps, target_grid)
    score = seismic_criterion.compute(distance)
    out_path = cache_processed(score, "seismic_score", scored_dir, overwrite=overwrite)
    echo(f"[done] seismic -> {out_path}")
    return out_path


CRITERION_FUNCS: dict[str, Callable[..., Path | None]] = {
    "slope": _slope_score,
    "illumination": _illumination_score,
    "hazard": _hazard_score,
    "thermal": _thermal_score,
    "ice": _ice_score,
    "seismic": _seismic_score,
}


def run(
    *,
    weights_path: Path = DEFAULT_WEIGHTS,
    region_config: Path = DEFAULT_REGION_CONFIG,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
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
        raw_dir: Where catalog files (e.g. scarps) live.
        overwrite: Re-compute even when cached score COGs exist.
        echo: Logging sink, defaults to ``typer.echo``.

    Returns:
        :class:`ScoreSummary` with the aggregate's cell stats and the
        list of criteria that ran versus skipped.
    """
    weights = _load_weights(weights_path)
    cfg = load_region_config(region_config)
    pixel_size_m = float(cfg["resolution_m"])
    processed_dir = Path(processed_dir)
    outputs_dir = Path(outputs_dir)
    raw_dir = Path(raw_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    score_paths: dict[str, Path] = {}
    skipped: list[str] = []
    for name, fn in CRITERION_FUNCS.items():
        path = fn(
            processed_dir,
            pixel_size_m=pixel_size_m,
            overwrite=overwrite,
            echo=echo,
            raw_dir=raw_dir,
        )
        if path is not None:
            score_paths[name] = path
        else:
            skipped.append(name)

    if not score_paths:
        raise RuntimeError(
            "no criterion score grids could be computed; run `selene preprocess` "
            "to populate data/processed/ first."
        )

    score_arrays: dict[str, xr.DataArray] = {
        name: _open_cog(path) for name, path in score_paths.items()
    }
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
        criteria_used=tuple(score_paths.keys()),
        criteria_skipped=tuple(skipped),
    )
    echo(summary.render())
    return summary
