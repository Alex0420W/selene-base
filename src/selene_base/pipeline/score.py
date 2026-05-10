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

v2.0 swaps the global spatial-coupling product (PSR-distance ×
ridge-distance) for the per-cell EVA-disc PSR access criterion
(:mod:`selene_base.criteria.eva_psr_access`), aligning with NASA's
documented walking-EVA framing (EVA-EXP-0070 Rev D, A3GT LPSC 2025).
v2.1 swaps the single ``ice`` criterion (PRP ice-depth + PSR-proximity
bonus) for the three-class volatile-access criterion
(:mod:`selene_base.criteria.multi_volatile`), distinguishing H₂O
(<110 K), CO₂/NH₃ (<66 K), and ultra-cold (<60 K) thermal classes
inside the same 2 km EVA disc. v2.2 swaps the v1.8 distance-to-
nearest-scarp seismic criterion for a per-cell PGV-style attenuation
kernel (:mod:`selene_base.criteria.pgv_seismic`) that aggregates
contributions from every mapped scarp within a 5L = 250 km cutoff
(L = 50 km, anchored to Watters 2024 PSJ strong-shaking distance).
All three legacy modules (coupling, ice, seismic) stay in tree for
sensitivity / comparison utilities but are no longer in the active
criterion set.

Week 3 covered the original six criteria; v1.5+ added LOS-to-Earth;
v1.8 activated seismic; v2.0 swapped coupling for eva_psr_access;
v2.1 swapped ice for multi_volatile; v2.2 swapped seismic for
pgv_seismic. Missing source data (Diviner, LEND, scarps) makes the
relevant criterion skip cleanly.
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
    eva_psr_access as eva_psr_access_criterion,
)
from selene_base.criteria import (
    hazard as hazard_criterion,
)
from selene_base.criteria import (
    illumination as illumination_criterion,
)
from selene_base.criteria import (
    los_to_earth as los_to_earth_criterion,
)
from selene_base.criteria import (
    multi_volatile as multi_volatile_criterion,
)
from selene_base.criteria import (
    pgv_seismic as pgv_seismic_criterion,
)
from selene_base.criteria import (
    slope as slope_criterion,
)
from selene_base.criteria import (
    thermal as thermal_criterion,
)
from selene_base.data.reproject import cache_processed
from selene_base.pipeline.preprocess import load_region_config
from selene_base.scoring.aggregate import AggregateMethod
from selene_base.scoring.aggregate import aggregate as _aggregate_fn

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


def _multi_volatile_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    temp_max_cog = processed_dir / "diviner_temp_max_southpole_240m.tif"
    if not temp_max_cog.exists():
        echo(
            "[skip] multi_volatile: Diviner PRP temp_max COG not present at "
            f"{temp_max_cog}; run `selene download diviner` then `selene preprocess`"
        )
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "multi_volatile_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] multi_volatile: {out_cog.name} already cached")
        return out_cog

    echo(
        f"[compute] multi_volatile from {temp_max_cog.name} "
        "(2 km EVA disc, H₂O<110 K + CO₂/NH₃<66 K + ultra-cold<60 K)"
    )
    temp_max = _open_cog(temp_max_cog)
    components = multi_volatile_criterion.compute_components(temp_max, pixel_size_m=pixel_size_m)
    # Persist all four arrays so per-class diagnostics survive without
    # rerunning preprocess+score: the combined score lands at the
    # canonical filename (the active criterion grid the aggregator
    # consumes); the three sub-scores land alongside as `_h2o`,
    # `_co2_nh3`, and `_ultracold` siblings for downstream reporting.
    out_path = cache_processed(
        components["combined_score"], "multi_volatile_score", scored_dir, overwrite=overwrite
    )
    cache_processed(
        components["h2o_score"], "multi_volatile_score_h2o", scored_dir, overwrite=overwrite
    )
    cache_processed(
        components["co2_nh3_score"],
        "multi_volatile_score_co2_nh3",
        scored_dir,
        overwrite=overwrite,
    )
    cache_processed(
        components["ultracold_score"],
        "multi_volatile_score_ultracold",
        scored_dir,
        overwrite=overwrite,
    )
    echo(f"[done] multi_volatile -> {out_path}")
    return out_path


def _eva_psr_access_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    temp_max_cog = processed_dir / "diviner_temp_max_southpole_240m.tif"
    if not temp_max_cog.exists():
        echo(
            "[skip] eva_psr_access: Diviner PRP temp_max COG not present at "
            f"{temp_max_cog}; run `selene download diviner` then `selene preprocess`"
        )
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "eva_psr_access_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] eva_psr_access: {out_cog.name} already cached")
        return out_cog

    echo(f"[compute] eva_psr_access from {temp_max_cog.name} (2 km EVA disc, T<110 K)")
    temp_max = _open_cog(temp_max_cog)
    score = eva_psr_access_criterion.compute(temp_max, pixel_size_m=pixel_size_m)
    out_path = cache_processed(score, "eva_psr_access_score", scored_dir, overwrite=overwrite)
    echo(f"[done] eva_psr_access -> {out_path}")
    return out_path


def _los_to_earth_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    visibility_cog = processed_dir / "los_visibility_fraction_southpole_240m.tif"
    if not visibility_cog.exists():
        echo(
            "[skip] los_to_earth: visibility-fraction COG not present at "
            f"{visibility_cog}; run `selene preprocess` first"
        )
        return None

    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "los_to_earth_score_southpole_240m.tif"
    if out_cog.exists() and not overwrite:
        echo(f"[skip] los_to_earth: {out_cog.name} already cached")
        return out_cog

    echo(f"[compute] los_to_earth from {visibility_cog.name}")
    visibility = _open_cog(visibility_cog)
    score = los_to_earth_criterion.compute(visibility)
    out_path = cache_processed(score, "los_to_earth_score", scored_dir, overwrite=overwrite)
    echo(f"[done] los_to_earth -> {out_path}")
    return out_path


def _pgv_seismic_score(
    processed_dir: Path,
    *,
    pixel_size_m: float,
    overwrite: bool,
    echo: Callable[[str], None],
    raw_dir: Path,
) -> Path | None:
    scored_dir = processed_dir / SCORED_SUBDIR
    out_cog = scored_dir / "pgv_seismic_score_southpole_240m.tif"

    # User-supplied overrides take precedence over the bundled default;
    # the in-repo Mishra & Kumar 2022 shapefile (v1.8+) is the fallback
    # so the criterion fires out of the box without any download.
    candidates = [
        raw_dir / "scarps" / "scarps_southpole.geojson",
        raw_dir / "scarps" / "scarps_southpole.csv",
        pgv_seismic_criterion.BUNDLED_MISHRA_KUMAR_2022,
    ]
    scarps_path = next((p for p in candidates if p.exists()), None)
    if scarps_path is None:
        echo(
            "[skip] pgv_seismic: no scarp catalog available — neither a "
            "user-supplied override at "
            f"{candidates[0].as_posix()} nor the bundled Mishra & Kumar "
            f"2022 shapefile at {candidates[2].as_posix()} could be found"
        )
        return None

    if out_cog.exists() and not overwrite:
        echo(f"[skip] pgv_seismic: {out_cog.name} already cached")
        return out_cog

    # Need a target grid in the projected CRS — borrow it from the LOLA COG.
    lola_cog = processed_dir / "lola_southpole_240m.tif"
    if not lola_cog.exists():
        echo(
            "[skip] pgv_seismic: needs a target grid; run `selene preprocess` "
            "to produce lola_southpole_240m.tif first"
        )
        return None

    label = (
        "Mishra & Kumar 2022 (bundled)"
        if scarps_path == pgv_seismic_criterion.BUNDLED_MISHRA_KUMAR_2022
        else scarps_path.name
    )
    echo(f"[compute] pgv_seismic from {label} (L=50 km, cutoff=250 km)")
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
    components = pgv_seismic_criterion.compute_components(scarps, target_grid)
    # Persist score + diagnostic siblings (cum_pgv, nearest_scarp_distance_km,
    # n_contributing_scarps) for downstream reporting.
    out_path = cache_processed(
        components["score"], "pgv_seismic_score", scored_dir, overwrite=overwrite
    )
    cache_processed(components["cum_pgv"], "pgv_seismic_cum_pgv", scored_dir, overwrite=overwrite)
    cache_processed(
        components["nearest_scarp_distance_km"],
        "pgv_seismic_nearest_scarp_distance_km",
        scored_dir,
        overwrite=overwrite,
    )
    cache_processed(
        components["n_contributing_scarps"].astype(np.float32),
        "pgv_seismic_n_contributing_scarps",
        scored_dir,
        overwrite=overwrite,
    )
    echo(f"[done] pgv_seismic -> {out_path}")
    return out_path


CRITERION_FUNCS: dict[str, Callable[..., Path | None]] = {
    "slope": _slope_score,
    "illumination": _illumination_score,
    "eva_psr_access": _eva_psr_access_score,
    "hazard": _hazard_score,
    "thermal": _thermal_score,
    "multi_volatile": _multi_volatile_score,
    "los_to_earth": _los_to_earth_score,
    "pgv_seismic": _pgv_seismic_score,
}


def run(
    *,
    weights_path: Path = DEFAULT_WEIGHTS,
    region_config: Path = DEFAULT_REGION_CONFIG,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    overwrite: bool = False,
    method: AggregateMethod = "weighted_sum",
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
    aggregate = _aggregate_fn(score_arrays, weights, method=method)

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
