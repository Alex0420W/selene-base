"""Tests for :mod:`selene_base.validation.wueller_comparison` (v1.4.1).

Starting in v1.4.1 the Wueller 2026 catalog ships in-repo as a real
130-site shapefile bundle from the authors' Zenodo deposit
(doi:10.5281/zenodo.17084058, CC-BY 4.0). Two formerly-skipped tests
(``test_real_wueller_comparison_against_v13_per_region_sites`` and
``test_per_region_wueller_counts_match_published_artemis_breakdown``)
are now active and exercise the real bundle.

The synthetic-input tests below remain — they verify the comparison
logic on inputs the test owns in-memory.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.wueller_comparison import (
    SYNTHETIC_PLACEHOLDER_PREFIX,
    WUELLER_CODE_TO_NAME,
    WUELLER_TO_USGS_REGION_MAP,
    compare_sites,
    is_synthetic_placeholder,
    load_wueller_sites,
)


def _selene_gdf(rows: list[dict]) -> gpd.GeoDataFrame:
    if not rows:
        return gpd.GeoDataFrame(
            columns=["site_id", "region_name", "lat", "lon"],
            geometry=[],
            crs=LUNAR_GEOGRAPHIC_CRS,
        )
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


def _wueller_gdf(rows: list[dict]) -> gpd.GeoDataFrame:
    if not rows:
        return gpd.GeoDataFrame(
            columns=["wueller_site_id", "region", "lat", "lon", "in_usgs_scope"],
            geometry=[],
            crs=LUNAR_GEOGRAPHIC_CRS,
        )
    for r in rows:
        r.setdefault("in_usgs_scope", True)
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


# ----------------------- load_wueller_sites + placeholder -------------------


def test_bundled_shapefile_loads_with_expected_schema() -> None:
    gdf = load_wueller_sites()
    assert {
        "wueller_site_id",
        "region",
        "wueller_region",
        "lat",
        "lon",
        "in_usgs_scope",
    }.issubset(gdf.columns)
    assert (gdf["lat"] < -80).all(), "Wueller sites should sit near the south pole"


def test_bundled_shapefile_is_not_synthetic_placeholder() -> None:
    """v1.4.1: real shapefile load is no longer flagged as a placeholder."""
    gdf = load_wueller_sites()
    assert not is_synthetic_placeholder(gdf)
    # Real Wueller site IDs are short codes like MMO01, NR101 — none
    # should carry the legacy synthetic prefix.
    assert not any(
        str(sid).startswith(SYNTHETIC_PLACEHOLDER_PREFIX) for sid in gdf["wueller_site_id"]
    )


def test_load_missing_source_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_wueller_sites(tmp_path / "nonexistent.shp")


def test_load_missing_columns_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.csv"
    bad.write_text("foo,bar\n1,2\n", encoding="utf-8")
    with pytest.warns(DeprecationWarning):
        with pytest.raises(ValueError, match="missing required columns"):
            load_wueller_sites(bad)


def test_load_reprojects_to_target_crs() -> None:
    polar = (
        "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 "
        "+R=1737400 +no_defs +type=crs"
    )
    gdf = load_wueller_sites(target_crs=polar)
    assert "stere" in str(gdf.crs.to_proj4())


def test_legacy_csv_load_emits_deprecation_warning() -> None:
    """The legacy synthetic CSV path remains for backward compat but warns."""
    from selene_base.validation.wueller_comparison import WUELLER_SITES_CSV

    if not WUELLER_SITES_CSV.exists():
        pytest.skip("Legacy CSV not present in this checkout")
    with pytest.warns(DeprecationWarning):
        gdf = load_wueller_sites(WUELLER_SITES_CSV)
    assert is_synthetic_placeholder(gdf)


# ----------------------- compare_sites: synthetic-only ----------------------


class TestEmptyInputs:
    def test_empty_selene(self) -> None:
        wueller = _wueller_gdf(
            [{"wueller_site_id": "synthetic-placeholder-1", "region": "R", "lat": -85, "lon": 0}]
        )
        result = compare_sites(_selene_gdf([]), wueller, filter_to_usgs_scope=False)
        assert result["n_selene_sites"] == 0
        assert result["n_wueller_sites"] == 1
        assert result["n_selene_matched"] == 0
        assert result["n_wueller_matched"] == 0
        assert np.isnan(result["median_match_distance_km"])

    def test_empty_wueller(self) -> None:
        selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85, "lon": 0}])
        result = compare_sites(selene, _wueller_gdf([]), filter_to_usgs_scope=False)
        assert result["n_selene_sites"] == 1
        assert result["n_wueller_sites"] == 0
        assert result["n_selene_matched"] == 0


class TestPerfectAgreement:
    def test_identical_sets(self) -> None:
        # Three sites at the exact same coordinates: every pair is a
        # zero-distance match and median distance is 0.
        rows = [
            {"site_id": 1, "region_name": "Mons Mouton", "lat": -85.3, "lon": 30.0},
            {"site_id": 2, "region_name": "Haworth", "lat": -86.5, "lon": -25.0},
            {"site_id": 3, "region_name": "Nobile Rim 1", "lat": -85.5, "lon": 37.0},
        ]
        wueller_rows = [
            {"wueller_site_id": "W1", "region": "Mons Mouton", "lat": -85.3, "lon": 30.0},
            {"wueller_site_id": "W2", "region": "Haworth", "lat": -86.5, "lon": -25.0},
            {"wueller_site_id": "W3", "region": "Nobile Rim 1", "lat": -85.5, "lon": 37.0},
        ]
        result = compare_sites(_selene_gdf(rows), _wueller_gdf(wueller_rows))
        assert result["n_selene_matched"] == 3
        assert result["n_wueller_matched"] == 3
        assert result["median_match_distance_km"] == pytest.approx(0.0, abs=1e-3)
        # Every per-region entry should be "n_matched == n_selene" with no
        # selene-only or wueller-only stragglers.
        for entry in result["per_region"]:
            assert entry["n_matched"] == entry["n_selene"] == entry["n_wueller"]
            assert entry["selene_only"] == []
            assert entry["wueller_only"] == []


class TestNoAgreement:
    def test_completely_disjoint(self) -> None:
        # selene site at south pole (-89, 0); Wueller site at the
        # opposite hemisphere far away. Nearest match is hundreds of km;
        # at the default 5 km threshold, nothing matches.
        selene = _selene_gdf([{"site_id": 1, "region_name": "Disjoint", "lat": -89.0, "lon": 0.0}])
        wueller = _wueller_gdf(
            [{"wueller_site_id": "W1", "region": "Disjoint", "lat": -82.0, "lon": 180.0}]
        )
        result = compare_sites(selene, wueller)
        assert result["n_selene_matched"] == 0
        assert result["n_wueller_matched"] == 0
        site = result["per_selene_site"][0]
        assert not site["matched"]
        assert site["distance_km"] > 100.0


class TestMatchThreshold:
    def test_threshold_boundary(self) -> None:
        # selene at (lat=-89, lon=0); two Wueller sites — one ~3 km
        # away, one ~8 km away. At threshold=5, only the closer one
        # matches; at threshold=10, both match.
        selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -89.0, "lon": 0.0}])
        # 0.1 deg of latitude near the pole ~= 3 km on a R=1737.4 km sphere.
        wueller = _wueller_gdf(
            [
                {"wueller_site_id": "W_close", "region": "R", "lat": -89.099, "lon": 0.0},
                {"wueller_site_id": "W_far", "region": "R", "lat": -89.27, "lon": 0.0},
            ]
        )
        # Distance from (lat -89, lon 0) to (lat -89.27, lon 0) is
        # ~0.27 deg * (pi/180) * 1737.4 km ~= 8.2 km on the lunar sphere.
        result_5 = compare_sites(selene, wueller, match_threshold_km=5.0)
        assert result_5["n_selene_matched"] == 1
        assert result_5["per_selene_site"][0]["nearest_wueller_id"] == "W_close"

        result_10 = compare_sites(selene, wueller, match_threshold_km=10.0)
        assert result_10["n_selene_matched"] == 1  # selene side is still 1 site
        assert result_10["n_wueller_matched"] == 2  # both Wueller sites within 10 km

    def test_distance_computation_known_offset(self) -> None:
        # Construct a pair we KNOW are 3 km apart along the meridian
        # at lat=-89: 3 km on a R=1737.4 km sphere = 3000 / 1737400 rad
        # = 0.0017269 rad = 0.09894 deg of latitude.
        delta_deg = 3000.0 / 1737400.0 * (180.0 / np.pi)
        selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -89.0, "lon": 0.0}])
        wueller = _wueller_gdf(
            [
                {
                    "wueller_site_id": "W",
                    "region": "R",
                    "lat": -89.0 - delta_deg,
                    "lon": 0.0,
                }
            ]
        )
        result = compare_sites(selene, wueller, match_threshold_km=5.0)
        # Stereographic at the pole is conformal — distance error vs
        # great-circle should be sub-percent for 3 km offsets.
        assert result["per_selene_site"][0]["distance_km"] == pytest.approx(3.0, rel=0.02)


class TestPerRegionAggregation:
    def test_cross_region_match_does_not_count_as_in_region(self) -> None:
        # selene site in region A is geometrically 1 km from a Wueller
        # site that is *labelled* region B. The pair is a "matched"
        # selene site in the global count (it's within threshold) but
        # the per-region "n_matched" must NOT count it for region A.
        # Use the same lon=0 meridian; 1 km offset.
        delta_deg = 1000.0 / 1737400.0 * (180.0 / np.pi)
        selene = _selene_gdf([{"site_id": 1, "region_name": "Region A", "lat": -89.0, "lon": 0.0}])
        wueller = _wueller_gdf(
            [
                {
                    "wueller_site_id": "W",
                    "region": "Region B",
                    "lat": -89.0 - delta_deg,
                    "lon": 0.0,
                }
            ]
        )
        result = compare_sites(selene, wueller, match_threshold_km=5.0)
        # Global match is True (the geometric pair is < 1.5 km apart).
        assert result["n_selene_matched"] == 1
        # Per-region: no region has both an A-selene and an A-Wueller,
        # so n_matched is 0 for every region.
        for entry in result["per_region"]:
            assert entry["n_matched"] == 0


class TestScopeFilter:
    def test_filter_drops_out_of_scope(self) -> None:
        # Two Wueller sites: one in-scope (Haworth), one out (Amundsen Rim).
        # With filter on, only the Haworth one is compared; n_wueller_sites
        # reflects the post-filter count and n_wueller_total carries the original.
        selene = _selene_gdf([{"site_id": 1, "region_name": "Haworth", "lat": -86.5, "lon": -25.0}])
        wueller = _wueller_gdf(
            [
                {
                    "wueller_site_id": "W_in",
                    "region": "Haworth",
                    "lat": -86.5,
                    "lon": -25.0,
                    "in_usgs_scope": True,
                },
                {
                    "wueller_site_id": "W_out",
                    "region": "Amundsen Rim",
                    "lat": -84.0,
                    "lon": 70.0,
                    "in_usgs_scope": False,
                },
            ]
        )
        result_filtered = compare_sites(selene, wueller, filter_to_usgs_scope=True)
        assert result_filtered["n_wueller_total"] == 2
        assert result_filtered["n_wueller_in_scope"] == 1
        assert result_filtered["n_wueller_out_of_scope"] == 1
        assert result_filtered["n_wueller_sites"] == 1
        assert result_filtered["scope_filter_applied"] is True
        assert result_filtered["out_of_scope_regions"] == ["Amundsen Rim"]

        result_unfiltered = compare_sites(selene, wueller, filter_to_usgs_scope=False)
        assert result_unfiltered["n_wueller_sites"] == 2
        assert result_unfiltered["scope_filter_applied"] is False


# ------------- placeholder-flag plumbing -------------


def test_synthetic_placeholder_flag_clear_for_non_placeholder_ids() -> None:
    selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85, "lon": 0}])
    wueller = _wueller_gdf([{"wueller_site_id": "W1", "region": "R", "lat": -85.0, "lon": 0.0}])
    result = compare_sites(selene, wueller)
    assert result["using_synthetic_placeholder"] is False


def test_synthetic_placeholder_flag_set_when_legacy_csv_loaded() -> None:
    """The legacy synthetic CSV still gets flagged when explicitly loaded."""
    from selene_base.validation.wueller_comparison import WUELLER_SITES_CSV

    if not WUELLER_SITES_CSV.exists():
        pytest.skip("Legacy CSV not present in this checkout")
    selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85.0, "lon": 0.0}])
    with pytest.warns(DeprecationWarning):
        wueller = load_wueller_sites(WUELLER_SITES_CSV)
    result = compare_sites(selene, wueller, filter_to_usgs_scope=False)
    assert result["using_synthetic_placeholder"] is True


# ---- real-data integration tests (v1.4.1 — bundle now ships) ----


def test_real_wueller_comparison_against_v13_per_region_sites() -> None:
    """Real Wueller 2026 catalog vs selene per-region sites.

    Loads selene per-region sites from the canonical
    ``data/outputs/per_region/sites.geojson`` artefact (skips the test
    if the artefact is not present, since it's gitignored). Asserts:

    - n_wueller_total == 130 (Wueller's headline count)
    - n_wueller_in_scope is bracketed at 60-80 (actual: 73)
    - n_selene_matched <= n_selene_sites
    - no NaN values in per-site distances
    - using_synthetic_placeholder is False
    """
    sites_path = Path("data/outputs/per_region/sites.geojson")
    if not sites_path.exists():
        pytest.skip(
            f"selene per-region sites not found at {sites_path}; "
            "run `selene rank-per-region` to generate them."
        )
    selene_sites = gpd.read_file(sites_path)
    wueller_sites = load_wueller_sites()
    result = compare_sites(
        selene_sites, wueller_sites, match_threshold_km=5.0, filter_to_usgs_scope=True
    )

    assert result["n_wueller_total"] == 130
    assert 60 <= result["n_wueller_in_scope"] <= 80, (
        f"in-scope count {result['n_wueller_in_scope']} outside expected range"
    )
    assert result["n_selene_matched"] <= result["n_selene_sites"]
    assert 0 <= result["n_selene_matched"] <= 70

    # No NaN distances on the per-site tables.
    for entry in result["per_selene_site"]:
        assert not np.isnan(entry["distance_km"])
    for entry in result["per_wueller_site"]:
        assert not np.isnan(entry["distance_km"])

    assert result["using_synthetic_placeholder"] is False


def test_per_region_wueller_counts_match_published_artemis_breakdown() -> None:
    """Per-region site counts as bundled in the v1.4.1 shapefile.

    Encodes the per-region in-scope counts from the Zenodo deposit so
    that schema drift in a future re-bundle is caught.
    """
    wueller_sites = load_wueller_sites()
    in_scope = wueller_sites[wueller_sites["in_usgs_scope"]]
    counts = in_scope.groupby("region")["wueller_site_id"].count().to_dict()
    expected = {
        "Haworth": 11,
        "Mons Mouton": 10,
        "Mons Mouton Plateau": 11,
        "Nobile Rim 1": 9,
        "Nobile Rim 2": 9,
        "Peak Near Cabeus B": 5,
        "Slater Plain": 11,
        "de Gerlache Rim 2": 7,
    }
    assert counts == expected, f"per-region in-scope counts drifted: got {counts}"
    assert sum(expected.values()) == 73


def test_region_code_mapping_round_trip() -> None:
    """Every shapefile region code maps to a Wueller name, which in turn
    is in WUELLER_TO_USGS_REGION_MAP."""
    wueller_sites = load_wueller_sites()
    codes = set(wueller_sites["wueller_region"].unique())
    for code in codes:
        assert code in WUELLER_CODE_TO_NAME, f"code {code!r} missing from WUELLER_CODE_TO_NAME"
        full = WUELLER_CODE_TO_NAME[code]
        assert full in WUELLER_TO_USGS_REGION_MAP, (
            f"name {full!r} missing from WUELLER_TO_USGS_REGION_MAP"
        )


def test_render_summary_does_not_flag_real_data() -> None:
    from selene_base.validation.wueller_comparison import render_summary

    selene = _selene_gdf([{"site_id": 1, "region_name": "Haworth", "lat": -86.5, "lon": -25.0}])
    wueller = load_wueller_sites()
    text = render_summary(compare_sites(selene, wueller))
    assert "SYNTHETIC PLACEHOLDER ACTIVE" not in text


def _smoke_use_pandas() -> None:
    """Sanity that pandas is importable for this test module."""
    pd.DataFrame({"a": [1]})
