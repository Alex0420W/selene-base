"""Extract the top-N candidate sites from the aggregate score map.

Loads ``data/outputs/score_southpole.tif`` plus any per-criterion score
COGs in ``data/processed/scored/``, runs non-maximum suppression at the
configured separation, and writes both a GeoJSON (full schema) and a
CSV (flat, human-friendly) to ``data/outputs/``.

Filled in week 3.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import rioxarray  # noqa: F401
import typer

from selene_base.scoring.ranking import (
    DEFAULT_CRITERIA,
    load_sub_scores,
    top_n_sites,
)

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
SCORED_SUBDIR = "scored"


def run(
    *,
    score_map_path: Path | None = None,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    top_n: int = 20,
    min_distance_km: float = 5.0,
    min_score: float = 0.5,
    echo: Callable[[str], None] = typer.echo,
) -> gpd.GeoDataFrame:
    """Rank the top sites and write GeoJSON + CSV.

    Args:
        score_map_path: Path to the aggregate score COG. Defaults to
            ``<outputs_dir>/score_southpole.tif``.
        processed_dir: Directory holding ``scored/<crit>_score_*.tif``.
        outputs_dir: Where ``top_sites.geojson`` and ``top_sites.csv``
            land.
        top_n: Maximum number of sites to extract.
        min_distance_km: Minimum pairwise separation, in km.
        min_score: Floor on candidate score.
        echo: Logging sink.

    Returns:
        The ranked GeoDataFrame (also written to disk).

    Raises:
        FileNotFoundError: When the aggregate score COG is missing.
    """
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    if score_map_path is None:
        score_map_path = outputs_dir / "score_southpole.tif"
    if not score_map_path.exists():
        raise FileNotFoundError(
            f"score map not found at {score_map_path}; run `selene score` first"
        )

    score_map = rioxarray.open_rasterio(score_map_path, masked=True).squeeze("band", drop=True)

    scored_dir = Path(processed_dir) / SCORED_SUBDIR
    sub_scores = load_sub_scores(scored_dir) if scored_dir.exists() else {}
    if sub_scores:
        echo(f"[rank] using sub-scores: {', '.join(sorted(sub_scores))}")
    else:
        echo("[rank] no per-criterion sub-scores found; geometry + total only")

    echo(f"[rank] top-{top_n}, min_distance={min_distance_km} km, min_score={min_score}")
    sites = top_n_sites(
        score_map,
        n=top_n,
        min_distance_m=min_distance_km * 1000.0,
        min_score=min_score,
        sub_scores=sub_scores,
    )

    geojson_path = outputs_dir / "top_sites.geojson"
    csv_path = outputs_dir / "top_sites.csv"

    if len(sites) == 0:
        echo("[rank] no sites met the threshold; writing empty outputs")
    sites.to_file(geojson_path, driver="GeoJSON")
    flat = sites.drop(columns="geometry").copy()
    flat.to_csv(csv_path, index=False)

    echo(f"[done] {len(sites)} site(s) -> {geojson_path}")
    echo(f"[done] {len(sites)} site(s) -> {csv_path}")
    if len(sites):
        echo("")
        echo(_format_table(sites))
    return sites


def _format_table(sites: gpd.GeoDataFrame) -> str:
    """Compact stdout table: rank | site_id | lat | lon | score | top criterion."""
    crit_cols = [f"score_{c}" for c in DEFAULT_CRITERIA if f"score_{c}" in sites.columns]
    lines = [
        f"{'rank':>4} {'site_id':<8} {'lat':>9} {'lon':>9} {'score':>6}  top criterion",
    ]
    for _, row in sites.iterrows():
        if crit_cols:
            sub = {c[len("score_") :]: row[c] for c in crit_cols}
            sub = {k: v for k, v in sub.items() if v == v}  # drop NaN
            top_crit = max(sub, key=sub.get) if sub else "—"
        else:
            top_crit = "—"
        lines.append(
            f"{int(row['rank']):>4} {row['site_id']:<8} "
            f"{row['lat']:>9.3f} {row['lon']:>9.3f} {row['score']:>6.3f}  {top_crit}"
        )
    return "\n".join(lines)
