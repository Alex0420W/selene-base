"""Tests for :mod:`selene_base.validation.wueller_comparison` (week 12 / v1.4.0).

The Wueller 2026 supplementary data is currently gated behind AGU /
Wiley and no open data release has been located. The bundled
``wueller_2026_sites.csv`` is a 5-row synthetic placeholder. These
tests exercise the comparison logic on synthetic inputs the test owns
in-memory; tests that would *only* be meaningful against the real
Wueller catalog are explicitly skipped with
``reason="awaiting upstream data"``.

Synthetic inputs cover:

- the empty-input edge cases (zero selene sites or zero Wueller sites),
- perfect-agreement (identical site sets at identical coordinates),
- a no-agreement pair (sites in completely different hemispheres),
- the match-threshold boundary (one pair just inside, one just outside),
- distance computation against a known km offset,
- the synthetic-placeholder detection on the bundled CSV.
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
            columns=["wueller_site_id", "region", "lat", "lon"],
            geometry=[],
            crs=LUNAR_GEOGRAPHIC_CRS,
        )
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


# ----------------------- load_wueller_sites + placeholder -------------------


def test_bundled_csv_loads_with_expected_schema() -> None:
    gdf = load_wueller_sites()
    assert {"wueller_site_id", "region", "lat", "lon"}.issubset(gdf.columns)
    assert (gdf["lat"] < -80).all(), "Wueller sites should sit near the south pole"


def test_bundled_csv_is_synthetic_placeholder() -> None:
    gdf = load_wueller_sites()
    assert is_synthetic_placeholder(gdf)
    assert all(sid.startswith(SYNTHETIC_PLACEHOLDER_PREFIX) for sid in gdf["wueller_site_id"])


def test_load_missing_csv_raises() -> None:
    with pytest.raises(FileNotFoundError):
        load_wueller_sites(Path("/nonexistent/wueller.csv"))


def test_load_missing_columns_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.csv"
    bad.write_text("foo,bar\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required columns"):
        load_wueller_sites(bad)


def test_load_reprojects_to_target_crs() -> None:
    polar = (
        "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 "
        "+R=1737400 +no_defs +type=crs"
    )
    gdf = load_wueller_sites(target_crs=polar)
    assert "stere" in str(gdf.crs.to_proj4())


# ----------------------- compare_sites: synthetic-only ----------------------


class TestEmptyInputs:
    def test_empty_selene(self) -> None:
        wueller = _wueller_gdf(
            [{"wueller_site_id": "synthetic-placeholder-1", "region": "R", "lat": -85, "lon": 0}]
        )
        result = compare_sites(_selene_gdf([]), wueller)
        assert result["n_selene_sites"] == 0
        assert result["n_wueller_sites"] == 1
        assert result["n_selene_matched"] == 0
        assert result["n_wueller_matched"] == 0
        assert np.isnan(result["median_match_distance_km"])

    def test_empty_wueller(self) -> None:
        selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85, "lon": 0}])
        result = compare_sites(selene, _wueller_gdf([]))
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


# ------------- placeholder-flag plumbing -------------


def test_synthetic_placeholder_flag_set_on_bundled_csv() -> None:
    selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85.0, "lon": 0.0}])
    wueller = load_wueller_sites()
    result = compare_sites(selene, wueller)
    assert result["using_synthetic_placeholder"] is True


def test_synthetic_placeholder_flag_clear_for_non_placeholder_ids() -> None:
    selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85, "lon": 0}])
    wueller = _wueller_gdf([{"wueller_site_id": "W1", "region": "R", "lat": -85.0, "lon": 0.0}])
    result = compare_sites(selene, wueller)
    assert result["using_synthetic_placeholder"] is False


# ---- gated-on-real-data tests: skipped until upstream data is available ----


@pytest.mark.skip(reason="awaiting upstream data: real Wueller 2026 catalog not yet acquired")
def test_real_wueller_comparison_against_v13_per_region_sites() -> None:
    """Placeholder for the real Wueller 2026 quantitative comparison.

    Once the AGU/Wiley supplementary Table A1 (or any open data release)
    is acquired and the bundled CSV is replaced with real coordinates,
    this test should:

    - Load selene v1.3 per-region sites from data/outputs/per_region/sites.geojson
    - Load real Wueller sites via load_wueller_sites
    - Call compare_sites at threshold=5 km
    - Assert: n_selene_sites > 0, n_wueller_sites == 130, and the
      result is *not* flagged as a synthetic placeholder.
    """


@pytest.mark.skip(
    reason="awaiting upstream data: per-region n_wueller counts depend on real catalog"
)
def test_per_region_wueller_counts_match_published_artemis_breakdown() -> None:
    """Once real data is in, verify the per-region site distribution matches
    the breakdown reported in the Wueller 2026 paper (e.g. how many of the
    130 sites fall in each Artemis CLR). This is a data-integrity smoke
    test that catches CSV-load errors after a future data refresh."""


# Sanity: keep an empty "real data" parquet shape as a sentinel so that
# the test discovery picks up the skipped placeholders even on fresh
# clones.
def test_skipped_real_data_tests_exist() -> None:
    from inspect import getmembers, isfunction

    import tests.test_wueller_comparison as module

    skipped = [
        name
        for name, fn in getmembers(module, isfunction)
        if getattr(fn, "pytestmark", [])
        and any(getattr(m, "name", "") == "skip" for m in fn.pytestmark)
    ]
    assert "test_real_wueller_comparison_against_v13_per_region_sites" in skipped
    assert "test_per_region_wueller_counts_match_published_artemis_breakdown" in skipped


def test_render_summary_flags_placeholder_run() -> None:
    from selene_base.validation.wueller_comparison import render_summary

    selene = _selene_gdf([{"site_id": 1, "region_name": "R", "lat": -85, "lon": 0}])
    wueller = load_wueller_sites()
    text = render_summary(compare_sites(selene, wueller))
    assert "SYNTHETIC PLACEHOLDER ACTIVE" in text
    assert "awaiting Wueller 2026 data release" in text


def _smoke_use_pandas() -> None:
    """Sanity that pandas is importable for this test module."""
    pd.DataFrame({"a": [1]})
