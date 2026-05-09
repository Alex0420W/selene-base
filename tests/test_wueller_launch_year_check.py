"""Tests for the v2.0 Wueller per-launch-year HLS cross-check."""

from __future__ import annotations

import math

import geopandas as gpd
import pytest
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.wueller_comparison import (
    ARTEMIS_IV_LAUNCH_YEAR,
    HLS_CONTINUOUS_LIT_THRESHOLD_DAYS,
    LAUNCH_YEARS,
    _resolve_sundays_columns,
    evaluate_per_launch_year,
    load_wueller_sites,
)


def _empty_comparison(n_selene: int = 0) -> dict:
    """Stub a WuellerComparisonResult-shaped dict for unit-level tests."""
    return {
        "n_selene_sites": n_selene,
        "n_wueller_sites": 0,
        "n_wueller_total": 0,
        "n_wueller_in_scope": 0,
        "n_wueller_out_of_scope": 0,
        "scope_filter_applied": True,
        "out_of_scope_regions": [],
        "n_selene_matched": 0,
        "n_wueller_matched": 0,
        "median_match_distance_km": float("nan"),
        "max_match_distance_km": float("nan"),
        "match_threshold_km": 5.0,
        "using_synthetic_placeholder": False,
        "per_region": [],
        "per_selene_site": [],
        "per_wueller_site": [],
    }


def _wueller_with_sundays(
    site_id: str, sundays_per_year: dict[int, int], lat: float = -89.5, lon: float = 0.0
) -> dict:
    """Build a Wueller-row dict carrying SunDays columns under the
    same case-mixed names the upstream DBF actually uses."""
    column_for_year = {
        2025: "SunDays25",
        2026: "SunDays26",
        2027: "SunDays27",
        2028: "SunDays28",
        2029: "SunDays29",
        2030: "Sundays30",  # lower-case d in the upstream DBF
        2031: "SunDays31",
        2032: "Sundays32",  # lower-case d in the upstream DBF
    }
    row = {
        "wueller_site_id": site_id,
        "region": "Haworth",
        "lat": lat,
        "lon": lon,
        "in_usgs_scope": True,
    }
    for year, col in column_for_year.items():
        row[col] = sundays_per_year[year]
    return row


def _gdf_from_rows(rows: list[dict]) -> gpd.GeoDataFrame:
    if not rows:
        return gpd.GeoDataFrame()
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


# ---- DBF schema ----


def test_loads_bundled_sundays_columns_for_all_eight_launch_years() -> None:
    """The bundled real shapefile must carry SunDays25..32; the resolver
    must locate one column per launch year regardless of letter case."""
    gdf = load_wueller_sites()
    cols = _resolve_sundays_columns(gdf)
    assert sorted(cols.keys()) == list(LAUNCH_YEARS)
    # Every resolved column name must actually be in the GeoDataFrame.
    for year, col_name in cols.items():
        assert col_name in gdf.columns, (year, col_name)


def test_resolver_handles_case_inconsistency() -> None:
    """The upstream DBF mixes ``SunDays`` and ``Sundays`` capitalisations;
    the resolver must accept either."""
    gdf = load_wueller_sites()
    cols = _resolve_sundays_columns(gdf)
    # Upstream uses lower-case 'd' for 2030 and 2032 specifically.
    # We don't pin the exact spelling — only that *some* column is found
    # for each year.
    assert 2030 in cols
    assert 2032 in cols


def test_resolver_raises_when_year_columns_missing() -> None:
    gdf = gpd.GeoDataFrame(
        [{"wueller_site_id": "x", "region": "Haworth", "SunDays25": 12}],
        geometry=[Point(0.0, -89.5)],
        crs=LUNAR_GEOGRAPHIC_CRS,
    )
    with pytest.raises(ValueError, match="missing SunDays columns"):
        _resolve_sundays_columns(gdf)


# ---- evaluation logic ----


def test_passes_hls_threshold_correctly_evaluates_sundays() -> None:
    # Site that scores 14 days at 2028 (>= 10 threshold) -> passes 2028.
    wueller_row = _wueller_with_sundays("MMO01", {y: 14 if y == 2028 else 0 for y in LAUNCH_YEARS})
    wueller_gdf = _gdf_from_rows([wueller_row])
    comparison = _empty_comparison(n_selene=1)
    comparison["n_selene_matched"] = 1
    comparison["per_selene_site"] = [
        {
            "site_id": "selene_HW_01",
            "region": "Haworth",
            "nearest_wueller_id": "MMO01",
            "distance_km": 1.4,
            "matched": True,
        }
    ]
    result = evaluate_per_launch_year(comparison, wueller_gdf)

    rec = result["per_site_records"][0]
    assert rec["passes_hls_2028"] is True
    assert rec["passes_hls_any_year_25_32"] is True
    assert rec["wueller_sundays"]["2028"] == 14
    # 2025-2027 + 2029-2032 are zero; only 2028 passes.
    assert result["aggregate"]["per_year_pass_count"]["2028"] == 1
    assert result["aggregate"]["per_year_pass_count"]["2025"] == 0


def test_below_threshold_correctly_fails() -> None:
    # SunDays at threshold-1 (=9) for every year: site fails every year.
    wueller_row = _wueller_with_sundays("X", {y: 9 for y in LAUNCH_YEARS})
    wueller_gdf = _gdf_from_rows([wueller_row])
    comparison = _empty_comparison(n_selene=1)
    comparison["n_selene_matched"] = 1
    comparison["per_selene_site"] = [
        {
            "site_id": "selene_HW_02",
            "region": "Haworth",
            "nearest_wueller_id": "X",
            "distance_km": 0.5,
            "matched": True,
        }
    ]
    result = evaluate_per_launch_year(comparison, wueller_gdf)
    rec = result["per_site_records"][0]
    assert rec["passes_hls_2028"] is False
    assert rec["passes_hls_any_year_25_32"] is False
    for y in LAUNCH_YEARS:
        assert result["aggregate"]["per_year_pass_count"][str(y)] == 0


def test_threshold_boundary_is_inclusive() -> None:
    # SunDays exactly at threshold (10) passes — the spec is "≥ 10".
    wueller_row = _wueller_with_sundays("X", {y: 10 for y in LAUNCH_YEARS})
    wueller_gdf = _gdf_from_rows([wueller_row])
    comparison = _empty_comparison(n_selene=1)
    comparison["n_selene_matched"] = 1
    comparison["per_selene_site"] = [
        {
            "site_id": "s1",
            "region": "Haworth",
            "nearest_wueller_id": "X",
            "distance_km": 0.5,
            "matched": True,
        }
    ]
    result = evaluate_per_launch_year(comparison, wueller_gdf)
    rec = result["per_site_records"][0]
    assert rec["passes_hls_2028"] is True
    for y in LAUNCH_YEARS:
        assert result["aggregate"]["per_year_pass_count"][str(y)] == 1


def test_handles_no_match_case() -> None:
    # selene site with no Wueller match (matched=False): passes_hls_*
    # are False, wueller_sundays is empty, distance is NaN.
    wueller_gdf = _gdf_from_rows([_wueller_with_sundays("X", {y: 30 for y in LAUNCH_YEARS})])
    comparison = _empty_comparison(n_selene=2)
    comparison["per_selene_site"] = [
        {
            "site_id": "matched_site",
            "region": "Haworth",
            "nearest_wueller_id": "X",
            "distance_km": 0.5,
            "matched": True,
        },
        {
            "site_id": "unmatched_site",
            "region": "Mons Mouton",
            "nearest_wueller_id": "X",
            "distance_km": 12.5,  # outside threshold
            "matched": False,
        },
    ]
    comparison["n_selene_matched"] = 1
    result = evaluate_per_launch_year(comparison, wueller_gdf)

    matched_rec = result["per_site_records"][0]
    unmatched_rec = result["per_site_records"][1]

    assert matched_rec["selene_site_id"] == "matched_site"
    assert matched_rec["passes_hls_2028"] is True

    assert unmatched_rec["selene_site_id"] == "unmatched_site"
    assert unmatched_rec["passes_hls_2028"] is False
    assert unmatched_rec["passes_hls_any_year_25_32"] is False
    assert unmatched_rec["matched_wueller_id"] == ""
    assert unmatched_rec["wueller_sundays"] == {}
    assert math.isnan(unmatched_rec["match_distance_km"])

    # Aggregate counts only matched sites.
    assert result["aggregate"]["selene_sites_with_wueller_match"] == 1
    assert result["aggregate"]["matched_sites_passing_hls_2028"] == 1
    assert result["aggregate"]["total_selene_sites"] == 2


def test_per_year_pass_count_aggregate() -> None:
    # Two matched selene sites with different SunDays profiles.
    # Site A passes 2025/2026; site B passes 2026/2027/2028.
    rows = [
        _wueller_with_sundays(
            "wA", {2025: 12, 2026: 11, 2027: 0, 2028: 0, 2029: 0, 2030: 0, 2031: 0, 2032: 0}
        ),
        _wueller_with_sundays(
            "wB",
            {2025: 0, 2026: 15, 2027: 13, 2028: 14, 2029: 0, 2030: 0, 2031: 0, 2032: 0},
            lat=-88.0,
        ),
    ]
    wueller_gdf = _gdf_from_rows(rows)
    comparison = _empty_comparison(n_selene=2)
    comparison["n_selene_matched"] = 2
    comparison["per_selene_site"] = [
        {
            "site_id": "sA",
            "region": "Haworth",
            "nearest_wueller_id": "wA",
            "distance_km": 0.5,
            "matched": True,
        },
        {
            "site_id": "sB",
            "region": "Haworth",
            "nearest_wueller_id": "wB",
            "distance_km": 0.5,
            "matched": True,
        },
    ]
    result = evaluate_per_launch_year(comparison, wueller_gdf)

    by_year = result["aggregate"]["per_year_pass_count"]
    assert by_year["2025"] == 1  # only A
    assert by_year["2026"] == 2  # both
    assert by_year["2027"] == 1  # only B
    assert by_year["2028"] == 1  # only B
    assert by_year["2029"] == 0
    assert result["aggregate"]["matched_sites_passing_hls_2028"] == 1
    assert result["aggregate"]["matched_sites_passing_hls_any_year"] == 2


def test_artemis_iv_2028_target_specifically_extracted() -> None:
    """``passes_hls_2028`` is the headline metric and must lock to 2028
    even if the user tweaks the threshold or LAUNCH_YEARS."""
    assert ARTEMIS_IV_LAUNCH_YEAR == 2028
    assert HLS_CONTINUOUS_LIT_THRESHOLD_DAYS == 10
    assert 2028 in LAUNCH_YEARS

    # Site that ONLY passes 2028 — every other year is below threshold.
    wueller_row = _wueller_with_sundays(
        "MMO01", {y: (15 if y == 2028 else 0) for y in LAUNCH_YEARS}
    )
    wueller_gdf = _gdf_from_rows([wueller_row])
    comparison = _empty_comparison(n_selene=1)
    comparison["n_selene_matched"] = 1
    comparison["per_selene_site"] = [
        {
            "site_id": "s",
            "region": "Haworth",
            "nearest_wueller_id": "MMO01",
            "distance_km": 1.0,
            "matched": True,
        }
    ]
    result = evaluate_per_launch_year(comparison, wueller_gdf)
    rec = result["per_site_records"][0]
    assert rec["passes_hls_2028"] is True
    assert rec["passes_hls_any_year_25_32"] is True
    assert result["aggregate"]["matched_sites_passing_hls_2028"] == 1


def test_metadata_block_is_populated() -> None:
    wueller_gdf = _gdf_from_rows([_wueller_with_sundays("X", {y: 0 for y in LAUNCH_YEARS})])
    comparison = _empty_comparison(n_selene=0)
    result = evaluate_per_launch_year(comparison, wueller_gdf, selene_version="v2.0")
    meta = result["metadata"]
    assert meta["selene_version"] == "v2.0"
    assert meta["hls_continuous_lit_threshold_days"] == 10
    assert meta["match_radius_km"] == 5.0
    # ISO 8601 UTC stamp (literal Z suffix from strftime).
    assert meta["comparison_date"].endswith("Z")
    assert "T" in meta["comparison_date"]


def test_negative_threshold_rejected() -> None:
    wueller_gdf = _gdf_from_rows([_wueller_with_sundays("X", {y: 0 for y in LAUNCH_YEARS})])
    with pytest.raises(ValueError, match="hls_threshold_days"):
        evaluate_per_launch_year(_empty_comparison(), wueller_gdf, hls_threshold_days=0)


def test_invalid_artemis_year_rejected() -> None:
    wueller_gdf = _gdf_from_rows([_wueller_with_sundays("X", {y: 0 for y in LAUNCH_YEARS})])
    with pytest.raises(ValueError, match="artemis_iv_year"):
        evaluate_per_launch_year(_empty_comparison(), wueller_gdf, artemis_iv_year=2099)
