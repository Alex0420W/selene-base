"""Quantitative comparison of selene-base sites against Wueller et al. 2026.

Wueller, F., et al. (2026), JGR Planets, doi:10.1029/2025JE009434, published
130 candidate Artemis III landing sites identified by the same outer
methodology selene-base implements: NASA HLS hard filters (slope < 8°,
100 m buffer to steeper terrain) followed by within-region selection.

This module ships the **comparison framework** as the v1.4.0 deliverable:

- :func:`load_wueller_sites` parses the CSV at
  :data:`WUELLER_SITES_CSV` (planetocentric lat/lon) and reprojects to
  any target CRS.
- :func:`compare_sites` does a pairwise nearest-neighbour match between
  the selene-base per-region sites (from
  :func:`selene_base.scoring.ranking.top_n_sites_per_region`) and the
  Wueller sites, returning per-site distances and per-region agreement
  counts.

The Wueller 2026 supplementary data is currently gated behind AGU/Wiley
and no open data release has been located (see README §"Data acquisition
status"). The CSV bundled in this repo at
``src/selene_base/validation/data/wueller_2026_sites.csv`` is a
**synthetic 5-row placeholder** (each row prefixed
``synthetic-placeholder-``) so the harness is runnable end-to-end
without external data; the CLI detects the placeholder by site_id
prefix and clearly labels its output as not a real scientific result.

When the real Table A1 / data release becomes available, replace the
CSV in place — the comparison logic, CLI, tests, and notebook will
produce real agreement numbers without further code changes.
"""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

WUELLER_SITES_CSV = Path(__file__).parent / "data" / "wueller_2026_sites.csv"
WUELLER_SOURCE_CRS = "+proj=longlat +R=1737400 +no_defs"
SYNTHETIC_PLACEHOLDER_PREFIX = "synthetic-placeholder-"
DEFAULT_MATCH_THRESHOLD_KM = 5.0

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


class PerRegionAgreement(TypedDict):
    region: str
    n_selene: int
    n_wueller: int
    n_matched: int
    median_distance_km: float
    selene_only: list[str]
    wueller_only: list[str]


class PerSeleneEntry(TypedDict):
    site_id: str
    region: str
    nearest_wueller_id: str
    distance_km: float
    matched: bool


class PerWuellerEntry(TypedDict):
    wueller_site_id: str
    region: str
    nearest_selene_id: str
    distance_km: float
    matched: bool


class WuellerComparisonResult(TypedDict):
    n_selene_sites: int
    n_wueller_sites: int
    n_selene_matched: int
    n_wueller_matched: int
    median_match_distance_km: float
    max_match_distance_km: float
    match_threshold_km: float
    using_synthetic_placeholder: bool
    per_region: list[PerRegionAgreement]
    per_selene_site: list[PerSeleneEntry]
    per_wueller_site: list[PerWuellerEntry]


def load_wueller_sites(
    csv_path: Path | None = None,
    *,
    target_crs: str | None = None,
) -> gpd.GeoDataFrame:
    """Load Wueller 2026 (or a placeholder stand-in) site coordinates.

    The CSV is read with comment-line stripping (lines starting ``#``);
    columns required: ``wueller_site_id``, ``region``, ``lat``,
    ``lon``. Coordinates are interpreted as planetocentric lat/lon on
    the lunar sphere (``R = 1737.4 km``); see :data:`WUELLER_SOURCE_CRS`.

    Args:
        csv_path: Override the bundled CSV path (mainly for tests).
        target_crs: Target CRS for the returned GeoDataFrame; defaults
            to lunar geographic lon/lat.

    Returns:
        GeoDataFrame with columns ``wueller_site_id``, ``region``,
        ``lat``, ``lon``, and a ``geometry`` Point in ``target_crs``.

    Raises:
        FileNotFoundError: If the CSV is missing.
        ValueError: If a required column is missing.
    """
    if csv_path is None:
        csv_path = WUELLER_SITES_CSV
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Wueller sites CSV not found at {csv_path}; "
            "place a file matching the schema described in "
            "selene_base.validation.wueller_comparison."
        )

    df = pd.read_csv(csv_path, comment="#")
    expected = {"wueller_site_id", "region", "lat", "lon"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"Wueller sites CSV missing required columns: {sorted(missing)}")

    df["wueller_site_id"] = df["wueller_site_id"].astype(str)
    df["region"] = df["region"].astype(str)
    geometries = [Point(lon, lat) for lon, lat in zip(df["lon"], df["lat"], strict=True)]
    gdf = gpd.GeoDataFrame(df, geometry=geometries, crs=WUELLER_SOURCE_CRS)
    if target_crs is not None and target_crs != WUELLER_SOURCE_CRS:
        gdf = gdf.to_crs(target_crs)
    return gdf


def is_synthetic_placeholder(wueller_sites: gpd.GeoDataFrame) -> bool:
    """True iff every site_id starts with ``synthetic-placeholder-``."""
    if len(wueller_sites) == 0:
        return False
    return all(
        str(sid).startswith(SYNTHETIC_PLACEHOLDER_PREFIX)
        for sid in wueller_sites["wueller_site_id"]
    )


def _project_xy(lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    xs, ys = transformer.transform(lons, lats)
    return np.column_stack([xs, ys])


def compare_sites(
    selene_sites: gpd.GeoDataFrame,
    wueller_sites: gpd.GeoDataFrame,
    *,
    match_threshold_km: float = DEFAULT_MATCH_THRESHOLD_KM,
) -> WuellerComparisonResult:
    """Pairwise nearest-neighbour comparison of selene-base vs Wueller sites.

    Distances are computed in lunar south-polar stereographic metres
    on the same projection used throughout the project, which is
    conformal at the pole and accurate to within a few percent for
    sub-100 km distances.

    The default ``match_threshold_km`` is 5 km — the upper end of the
    1-5 km granularity at which regional candidate site selection
    typically operates (NASA HLS landing accuracy is 100 m, but
    candidate site *selection* is regional). The threshold is exposed
    rather than hard-coded so tests and downstream tooling can exercise
    sensitivity, but the default is *not* tuned from the agreement
    result.

    Args:
        selene_sites: GeoDataFrame from
            :func:`selene_base.scoring.ranking.top_n_sites_per_region`,
            with at minimum ``site_id``, ``region_name``, ``lat``,
            ``lon`` columns.
        wueller_sites: GeoDataFrame from :func:`load_wueller_sites`,
            with ``wueller_site_id``, ``region``, ``lat``, ``lon``.
        match_threshold_km: Distance below which a pair is marked
            "matched". Default 5 km.

    Returns:
        :class:`WuellerComparisonResult` with headline counts,
        per-region agreement, and per-site nearest-neighbour tables.
        The ``using_synthetic_placeholder`` flag is True when every
        ``wueller_site_id`` starts with ``synthetic-placeholder-``;
        callers should treat the agreement numbers as a harness sanity
        check, not as a scientific result, when this flag is set.
    """
    n_selene = int(len(selene_sites))
    n_wueller = int(len(wueller_sites))
    using_placeholder = is_synthetic_placeholder(wueller_sites)
    threshold_m = float(match_threshold_km) * 1000.0

    if n_selene == 0 or n_wueller == 0:
        per_region: list[PerRegionAgreement] = []
        per_selene_site: list[PerSeleneEntry] = []
        per_wueller_site: list[PerWuellerEntry] = []
        for _, row in selene_sites.iterrows():
            per_selene_site.append(
                {
                    "site_id": str(row.get("site_id", "")),
                    "region": str(row.get("region_name", "")),
                    "nearest_wueller_id": "",
                    "distance_km": float("nan"),
                    "matched": False,
                }
            )
        for _, row in wueller_sites.iterrows():
            per_wueller_site.append(
                {
                    "wueller_site_id": str(row.get("wueller_site_id", "")),
                    "region": str(row.get("region", "")),
                    "nearest_selene_id": "",
                    "distance_km": float("nan"),
                    "matched": False,
                }
            )
        return {
            "n_selene_sites": n_selene,
            "n_wueller_sites": n_wueller,
            "n_selene_matched": 0,
            "n_wueller_matched": 0,
            "median_match_distance_km": float("nan"),
            "max_match_distance_km": float("nan"),
            "match_threshold_km": float(match_threshold_km),
            "using_synthetic_placeholder": using_placeholder,
            "per_region": per_region,
            "per_selene_site": per_selene_site,
            "per_wueller_site": per_wueller_site,
        }

    selene_xy = _project_xy(
        selene_sites["lon"].to_numpy(dtype=np.float64),
        selene_sites["lat"].to_numpy(dtype=np.float64),
    )
    wueller_xy = _project_xy(
        wueller_sites["lon"].to_numpy(dtype=np.float64),
        wueller_sites["lat"].to_numpy(dtype=np.float64),
    )

    selene_tree = cKDTree(selene_xy)
    wueller_tree = cKDTree(wueller_xy)

    selene_to_w_dist_m, selene_to_w_idx = wueller_tree.query(selene_xy, k=1)
    wueller_to_s_dist_m, wueller_to_s_idx = selene_tree.query(wueller_xy, k=1)

    selene_matched = selene_to_w_dist_m <= threshold_m
    wueller_matched = wueller_to_s_dist_m <= threshold_m

    # Per-site tables.
    per_selene_site = []
    selene_ids = selene_sites["site_id"].astype(str).to_numpy()
    selene_regions = selene_sites["region_name"].astype(str).to_numpy()
    wueller_ids = wueller_sites["wueller_site_id"].astype(str).to_numpy()
    wueller_regions = wueller_sites["region"].astype(str).to_numpy()
    for i in range(n_selene):
        per_selene_site.append(
            {
                "site_id": str(selene_ids[i]),
                "region": str(selene_regions[i]),
                "nearest_wueller_id": str(wueller_ids[selene_to_w_idx[i]]),
                "distance_km": float(selene_to_w_dist_m[i]) / 1000.0,
                "matched": bool(selene_matched[i]),
            }
        )
    per_wueller_site = []
    for j in range(n_wueller):
        per_wueller_site.append(
            {
                "wueller_site_id": str(wueller_ids[j]),
                "region": str(wueller_regions[j]),
                "nearest_selene_id": str(selene_ids[wueller_to_s_idx[j]]),
                "distance_km": float(wueller_to_s_dist_m[j]) / 1000.0,
                "matched": bool(wueller_matched[j]),
            }
        )

    matched_distances_km = selene_to_w_dist_m[selene_matched] / 1000.0
    median_match_km = (
        float(np.median(matched_distances_km)) if matched_distances_km.size > 0 else float("nan")
    )
    max_match_km = (
        float(np.max(matched_distances_km)) if matched_distances_km.size > 0 else float("nan")
    )

    # Per-region aggregation across the union of region names that
    # appear on either side.
    regions_in_either = sorted(set(selene_regions.tolist()) | set(wueller_regions.tolist()))
    per_region = []
    for region in regions_in_either:
        selene_mask = selene_regions == region
        wueller_mask = wueller_regions == region
        n_selene_r = int(selene_mask.sum())
        n_wueller_r = int(wueller_mask.sum())
        # Match count for this region: count selene sites in this region whose
        # nearest Wueller match (anywhere) is within threshold AND that match
        # is also in this region. This keeps "matched within region" honest;
        # cross-region matches don't count.
        n_matched_r = 0
        match_dists_r: list[float] = []
        selene_only: list[str] = []
        for i, in_region in enumerate(selene_mask):
            if not in_region:
                continue
            j = int(selene_to_w_idx[i])
            same_region_match = wueller_regions[j] == region
            close_enough = bool(selene_matched[i])
            if same_region_match and close_enough:
                n_matched_r += 1
                match_dists_r.append(float(selene_to_w_dist_m[i]) / 1000.0)
            elif n_selene_r > 0:
                selene_only.append(str(selene_ids[i]))
        wueller_only: list[str] = []
        for j, in_region in enumerate(wueller_mask):
            if not in_region:
                continue
            i_match = int(wueller_to_s_idx[j])
            same_region_match = selene_regions[i_match] == region
            close_enough = bool(wueller_matched[j])
            if not (same_region_match and close_enough):
                wueller_only.append(str(wueller_ids[j]))
        per_region.append(
            {
                "region": region,
                "n_selene": n_selene_r,
                "n_wueller": n_wueller_r,
                "n_matched": n_matched_r,
                "median_distance_km": (
                    float(np.median(match_dists_r)) if match_dists_r else float("nan")
                ),
                "selene_only": selene_only,
                "wueller_only": wueller_only,
            }
        )

    return {
        "n_selene_sites": n_selene,
        "n_wueller_sites": n_wueller,
        "n_selene_matched": int(selene_matched.sum()),
        "n_wueller_matched": int(wueller_matched.sum()),
        "median_match_distance_km": median_match_km,
        "max_match_distance_km": max_match_km,
        "match_threshold_km": float(match_threshold_km),
        "using_synthetic_placeholder": using_placeholder,
        "per_region": per_region,
        "per_selene_site": per_selene_site,
        "per_wueller_site": per_wueller_site,
    }


def render_summary(result: WuellerComparisonResult) -> str:
    """Stdout-ready three-block summary for ``selene compare-wueller``."""
    lines: list[str] = []
    n_s = result["n_selene_sites"]
    n_w = result["n_wueller_sites"]
    threshold = result["match_threshold_km"]

    if result["using_synthetic_placeholder"]:
        lines.append(
            "*** SYNTHETIC PLACEHOLDER ACTIVE ***  "
            "Comparison harness is running against a synthetic 5-row stand-in;"
        )
        lines.append(
            "   awaiting Wueller 2026 data release before this command produces "
            "a real scientific result."
        )
        lines.append("")

    lines.append(
        f"selene-base sites: {n_s}    Wueller 2026 sites: {n_w}    "
        f"match threshold: {threshold:.1f} km"
    )
    lines.append(f"selene matched (any region) : {result['n_selene_matched']:>3} / {n_s}")
    lines.append(f"Wueller matched (any region): {result['n_wueller_matched']:>3} / {n_w}")
    median_km = result["median_match_distance_km"]
    if median_km == median_km:  # not NaN
        lines.append(f"median matched-pair distance: {median_km:.2f} km")
    else:
        lines.append("median matched-pair distance: n/a (no matches)")

    lines.append("")
    lines.append("per-region agreement:")
    lines.append(
        f"{'region':<22} {'selene n':>8} {'wueller n':>9} {'matched':>7} {'median dist (km)':>16}"
    )
    for entry in result["per_region"]:
        med = entry["median_distance_km"]
        med_str = f"{med:>16.2f}" if med == med else f"{'--':>16}"
        lines.append(
            f"{entry['region']:<22} {entry['n_selene']:>8} {entry['n_wueller']:>9} "
            f"{entry['n_matched']:>7} {med_str}"
        )

    # Notable disagreements: up to 3 selene-only and 3 wueller-only entries
    # with the largest distance to the other set.
    lines.append("")
    lines.append("notable disagreements (largest nearest-distance, capped at 6):")
    sorted_selene = sorted(
        (e for e in result["per_selene_site"] if not e["matched"]),
        key=lambda e: -e["distance_km"] if e["distance_km"] == e["distance_km"] else 0,
    )[:3]
    sorted_wueller = sorted(
        (e for e in result["per_wueller_site"] if not e["matched"]),
        key=lambda e: -e["distance_km"] if e["distance_km"] == e["distance_km"] else 0,
    )[:3]
    if sorted_selene:
        lines.append("  selene-only:")
        for e in sorted_selene:
            lines.append(
                f"    site {e['site_id']:<3} ({e['region']}): "
                f"nearest Wueller {e['nearest_wueller_id']} at {e['distance_km']:.1f} km"
            )
    if sorted_wueller:
        lines.append("  Wueller-only:")
        for e in sorted_wueller:
            lines.append(
                f"    {e['wueller_site_id']} ({e['region']}): "
                f"nearest selene site {e['nearest_selene_id']} at {e['distance_km']:.1f} km"
            )

    if result["using_synthetic_placeholder"]:
        lines.append("")
        lines.append(
            "*** Reminder: numbers above are computed against the synthetic "
            "placeholder, not real Wueller 2026 data. ***"
        )
    return "\n".join(lines)
