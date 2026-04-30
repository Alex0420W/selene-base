"""Quantitative comparison of selene-base sites against Wueller et al. 2026.

Wueller, F., et al. (2026), JGR Planets, doi:10.1029/2025JE009434, published
130 candidate Artemis III landing sites identified by the same outer
methodology selene-base implements: NASA HLS hard filters (slope < 8°,
100 m buffer to steeper terrain) followed by within-region selection.

Starting in v1.4.1 the comparison runs against the **real** 130-site
catalog from the authors' Zenodo deposit
(`Complementary Data for Wueller et al. (2026)`,
doi:10.5281/zenodo.17084058, CC-BY 4.0). The bundled shapefile lives
at ``src/selene_base/validation/data/wueller_2026/LandingSites.shp``;
:func:`load_wueller_sites` reads it by default and falls back to the
v1.4.0 synthetic CSV only if the shapefile is absent (the fallback
emits a deprecation warning).

- :func:`load_wueller_sites` reads the bundled shapefile (or the legacy
  CSV fallback) and reprojects to any target CRS. The returned
  GeoDataFrame carries the project schema (``wueller_site_id``,
  ``region``, ``lat``, ``lon``) plus a ``wueller_region`` 3-letter code
  column, an ``in_usgs_scope`` bool, and all upstream attribute columns
  (e.g. ``SunDays25``–``SunDays32``, ``PSR_AREA``) preserved verbatim.
- :func:`compare_sites` does a pairwise nearest-neighbour match. By
  default it filters Wueller to the 73 sites within NASA's October 2024
  down-selected nine regions; pass ``filter_to_usgs_scope=False`` to
  compare against all 130.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TypedDict

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

WUELLER_SITES_SHAPEFILE = Path(__file__).parent / "data" / "wueller_2026" / "LandingSites.shp"
WUELLER_SITES_CSV = Path(__file__).parent / "data" / "wueller_2026_sites.csv"
WUELLER_SOURCE_CRS = "+proj=longlat +R=1737400 +no_defs"
SYNTHETIC_PLACEHOLDER_PREFIX = "synthetic-placeholder-"
DEFAULT_MATCH_THRESHOLD_KM = 5.0

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)

# 3-letter codes used in the Wueller shapefile's ``Landing_Re`` column,
# expanded to the full Wueller region name.
WUELLER_CODE_TO_NAME: dict[str, str] = {
    "AR": "Amundsen Rim",
    "CR": "Connecting Ridge",
    "CRE": "Connecting Ridge Extension",
    "FRA": "Faustini Rim A",
    "HW": "Haworth",
    "MMA": "Mons Malapert",
    "MMO": "Mons Mouton",
    "MMP": "Mons Mouton Plateau",
    "NR1": "Nobile Rim 1",
    "NR2": "Nobile Rim 2",
    "PCB": "Peak Near Cabeus B",
    "PNS": "Peak Near Shackleton",
    "SP": "Slater Plain",
    "dGKM": "de Gerlache-Kocher Massif",
    "dGR1": "de Gerlache Rim 1",
    "dGR2": "de Gerlache Rim 2",
}

# Wueller -> NASA October-2024 down-selected USGS canonical name.
# ``None`` = Wueller region was on the earlier 13-region list but is
# not in NASA's final nine. Mons Malapert is intentionally distinct
# from Malapert Massif (which is the USGS-scope region) — they sit on
# different terrain.
WUELLER_TO_USGS_REGION_MAP: dict[str, str | None] = {
    "Amundsen Rim": None,
    "Connecting Ridge": None,
    "Connecting Ridge Extension": None,
    "Faustini Rim A": None,
    "Haworth": "Haworth",
    "Mons Malapert": None,
    "Mons Mouton": "Mons Mouton",
    "Mons Mouton Plateau": "Mons Mouton Plateau",
    "Nobile Rim 1": "Nobile Rim 1",
    "Nobile Rim 2": "Nobile Rim 2",
    "Peak Near Cabeus B": "Peak Near Cabeus B",
    "Peak Near Shackleton": None,
    "Slater Plain": "Slater Plain",
    "de Gerlache-Kocher Massif": None,
    "de Gerlache Rim 1": None,
    "de Gerlache Rim 2": "de Gerlache Rim 2",
}


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
    n_wueller_total: int
    n_wueller_in_scope: int
    n_wueller_out_of_scope: int
    scope_filter_applied: bool
    out_of_scope_regions: list[str]
    n_selene_matched: int
    n_wueller_matched: int
    median_match_distance_km: float
    max_match_distance_km: float
    match_threshold_km: float
    using_synthetic_placeholder: bool
    per_region: list[PerRegionAgreement]
    per_selene_site: list[PerSeleneEntry]
    per_wueller_site: list[PerWuellerEntry]


def _load_from_shapefile(shapefile_path: Path, target_crs: str | None) -> gpd.GeoDataFrame:
    """Load the bundled Wueller 2026 shapefile and rename to project schema."""
    gdf = gpd.read_file(shapefile_path)
    if "Name" not in gdf.columns or "Landing_Re" not in gdf.columns:
        raise ValueError(
            f"Wueller shapefile at {shapefile_path} missing expected columns "
            "'Name' and/or 'Landing_Re'; got "
            f"{sorted(gdf.columns.tolist())}"
        )

    gdf = gdf.rename(
        columns={
            "Name": "wueller_site_id",
            "Landing_Re": "wueller_region",
            "Latitude": "lat",
            "Longitude": "lon",
        }
    )
    gdf["wueller_site_id"] = gdf["wueller_site_id"].astype(str)
    gdf["wueller_region"] = gdf["wueller_region"].astype(str)

    # Expand 3-letter codes to full names; map to USGS canonical name
    # where one exists; mark in_usgs_scope.
    full_names = gdf["wueller_region"].map(WUELLER_CODE_TO_NAME)
    if full_names.isna().any():
        unknown = sorted(gdf.loc[full_names.isna(), "wueller_region"].unique().tolist())
        raise ValueError(
            f"Wueller shapefile contains unknown region code(s): {unknown}. "
            "Update WUELLER_CODE_TO_NAME if this is a new release."
        )
    usgs_names = full_names.map(WUELLER_TO_USGS_REGION_MAP)
    gdf["region"] = usgs_names.where(usgs_names.notna(), full_names)
    gdf["in_usgs_scope"] = usgs_names.notna()

    if target_crs is not None:
        gdf = gdf.to_crs(target_crs)
    return gdf


def _load_from_csv(csv_path: Path, target_crs: str | None) -> gpd.GeoDataFrame:
    """Legacy CSV loader for the v1.4.0 synthetic placeholder."""
    df = pd.read_csv(csv_path, comment="#")
    expected = {"wueller_site_id", "region", "lat", "lon"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"Wueller sites CSV missing required columns: {sorted(missing)}")

    df["wueller_site_id"] = df["wueller_site_id"].astype(str)
    df["region"] = df["region"].astype(str)
    df["wueller_region"] = df["region"]
    df["in_usgs_scope"] = df["region"].isin(set(WUELLER_TO_USGS_REGION_MAP.values()) - {None})

    geometries = [Point(lon, lat) for lon, lat in zip(df["lon"], df["lat"], strict=True)]
    gdf = gpd.GeoDataFrame(df, geometry=geometries, crs=WUELLER_SOURCE_CRS)
    if target_crs is not None and target_crs != WUELLER_SOURCE_CRS:
        gdf = gdf.to_crs(target_crs)
    return gdf


def load_wueller_sites(
    source_path: Path | None = None,
    *,
    target_crs: str | None = None,
) -> gpd.GeoDataFrame:
    """Load the Wueller 2026 site catalog.

    Default load path (v1.4.1+): the bundled shapefile at
    ``src/selene_base/validation/data/wueller_2026/LandingSites.shp``,
    which carries 130 real sites from the authors' Zenodo deposit
    (CC-BY 4.0). The loader renames a small subset of upstream columns
    to match the project schema (``wueller_site_id``, ``region``,
    ``lat``, ``lon``) and adds two derived columns
    (``wueller_region`` 3-letter code, ``in_usgs_scope`` bool); all
    other upstream attribute columns are preserved.

    Fallback (deprecated): if ``source_path`` points to a ``.csv`` (or
    if the bundled shapefile is missing and only the legacy synthetic
    CSV is present), the loader reads that CSV instead and emits a
    ``DeprecationWarning``. This path remains for backward compatibility
    with anyone who has the v1.4.0 placeholder CSV in place.

    Args:
        source_path: Override path. Accepts a ``.shp`` shapefile or a
            ``.csv`` file. ``None`` (default) tries the bundled
            shapefile first, then the legacy CSV.
        target_crs: Target CRS for the returned GeoDataFrame. ``None``
            keeps the source CRS (polar stereographic for the
            shapefile; lunar geographic for the CSV).

    Returns:
        GeoDataFrame with at minimum ``wueller_site_id``, ``region``,
        ``wueller_region``, ``lat``, ``lon``, ``in_usgs_scope``, and a
        ``geometry`` Point in the chosen CRS.

    Raises:
        FileNotFoundError: If no source data is available.
        ValueError: If the source schema is unrecognised.
    """
    if source_path is None:
        if WUELLER_SITES_SHAPEFILE.exists():
            return _load_from_shapefile(WUELLER_SITES_SHAPEFILE, target_crs)
        if WUELLER_SITES_CSV.exists():
            warnings.warn(
                "Loading Wueller sites from the legacy synthetic CSV at "
                f"{WUELLER_SITES_CSV}; the v1.4.1 shapefile bundle at "
                f"{WUELLER_SITES_SHAPEFILE} is missing. The synthetic "
                "fallback exists for backward compatibility only.",
                DeprecationWarning,
                stacklevel=2,
            )
            return _load_from_csv(WUELLER_SITES_CSV, target_crs)
        raise FileNotFoundError(
            "No Wueller site data available: neither "
            f"{WUELLER_SITES_SHAPEFILE} nor {WUELLER_SITES_CSV} exists."
        )

    source_path = Path(source_path)
    if not source_path.exists():
        raise FileNotFoundError(f"Wueller sites source not found at {source_path}")

    if source_path.suffix.lower() == ".shp":
        return _load_from_shapefile(source_path, target_crs)
    if source_path.suffix.lower() == ".csv":
        warnings.warn(
            "Loading Wueller sites from a CSV path; the v1.4.1 default "
            "is the bundled shapefile. CSV mode is retained for "
            "backward compatibility only.",
            DeprecationWarning,
            stacklevel=2,
        )
        return _load_from_csv(source_path, target_crs)
    raise ValueError(
        f"Unsupported Wueller source extension {source_path.suffix!r}; expected .shp or .csv"
    )


def is_synthetic_placeholder(wueller_sites: gpd.GeoDataFrame) -> bool:
    """True iff every site_id starts with ``synthetic-placeholder-``.

    Returns False when the v1.4.1 real shapefile is loaded (real
    Wueller IDs are short codes like ``MMO01``, ``NR101``).
    """
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
    filter_to_usgs_scope: bool = True,
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

    When ``filter_to_usgs_scope`` is True (default), Wueller sites
    outside NASA's October 2024 down-selected nine regions are dropped
    before the comparison. This is the apples-to-apples mode: selene-base
    only ranks within USGS-scope regions, so comparing against Wueller's
    out-of-scope sites would be a structural mismatch, not a methodology
    disagreement. The result dict reports both the in-scope count
    actually used and the total Wueller catalog count for context.
    Pass ``filter_to_usgs_scope=False`` to use all 130 sites.

    Args:
        selene_sites: GeoDataFrame from
            :func:`selene_base.scoring.ranking.top_n_sites_per_region`,
            with at minimum ``site_id``, ``region_name``, ``lat``,
            ``lon`` columns.
        wueller_sites: GeoDataFrame from :func:`load_wueller_sites`,
            with ``wueller_site_id``, ``region``, ``lat``, ``lon``,
            and (for shapefile-sourced loads) an ``in_usgs_scope`` bool.
        match_threshold_km: Distance below which a pair is marked
            "matched". Default 5 km.
        filter_to_usgs_scope: When True (default), drop Wueller sites
            with ``in_usgs_scope == False`` before comparing.

    Returns:
        :class:`WuellerComparisonResult` with headline counts,
        per-region agreement, and per-site nearest-neighbour tables.
        ``n_wueller_total`` is the count before any scope filter;
        ``n_wueller_in_scope`` is the count after; ``n_wueller_sites``
        is the count actually used in the comparison (equal to
        ``n_wueller_in_scope`` when the filter is applied, else
        ``n_wueller_total``).
    """
    n_wueller_total = int(len(wueller_sites))
    if "in_usgs_scope" in wueller_sites.columns:
        in_scope_mask = wueller_sites["in_usgs_scope"].astype(bool).to_numpy()
    else:
        in_scope_mask = np.ones(n_wueller_total, dtype=bool)
    n_wueller_in_scope = int(in_scope_mask.sum())
    n_wueller_out_of_scope = n_wueller_total - n_wueller_in_scope

    if filter_to_usgs_scope:
        out_of_scope_regions = sorted(
            set(wueller_sites.loc[~in_scope_mask, "region"].astype(str).tolist())
        )
        wueller_used = wueller_sites.loc[in_scope_mask].reset_index(drop=True)
    else:
        out_of_scope_regions = []
        wueller_used = wueller_sites.reset_index(drop=True)

    n_selene = int(len(selene_sites))
    n_wueller = int(len(wueller_used))
    using_placeholder = is_synthetic_placeholder(wueller_used)
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
        for _, row in wueller_used.iterrows():
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
            "n_wueller_total": n_wueller_total,
            "n_wueller_in_scope": n_wueller_in_scope,
            "n_wueller_out_of_scope": n_wueller_out_of_scope,
            "scope_filter_applied": bool(filter_to_usgs_scope),
            "out_of_scope_regions": out_of_scope_regions,
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
        wueller_used["lon"].to_numpy(dtype=np.float64),
        wueller_used["lat"].to_numpy(dtype=np.float64),
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
    wueller_ids = wueller_used["wueller_site_id"].astype(str).to_numpy()
    wueller_regions = wueller_used["region"].astype(str).to_numpy()
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
        "n_wueller_total": n_wueller_total,
        "n_wueller_in_scope": n_wueller_in_scope,
        "n_wueller_out_of_scope": n_wueller_out_of_scope,
        "scope_filter_applied": bool(filter_to_usgs_scope),
        "out_of_scope_regions": out_of_scope_regions,
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
            "Comparison harness is running against a synthetic stand-in;"
        )
        lines.append(
            "   the v1.4.1 default loads the real Zenodo shapefile bundle "
            "— this run is not a real scientific result."
        )
        lines.append("")

    lines.append(
        f"selene-base sites: {n_s}    Wueller 2026 sites: {n_w}    "
        f"match threshold: {threshold:.1f} km"
    )
    if result["scope_filter_applied"]:
        lines.append(
            f"  (Wueller catalog total: {result['n_wueller_total']}; "
            f"in USGS scope: {result['n_wueller_in_scope']}; "
            f"out-of-scope dropped: {result['n_wueller_out_of_scope']})"
        )
    else:
        lines.append(
            f"  (full catalog mode: scope filter not applied; "
            f"in-scope: {result['n_wueller_in_scope']} / "
            f"{result['n_wueller_total']})"
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
            "fallback, not the real Wueller 2026 shapefile. ***"
        )
    return "\n".join(lines)
