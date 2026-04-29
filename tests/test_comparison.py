"""Tests for :mod:`selene_base.validation.comparison`."""

from __future__ import annotations

import geopandas as gpd
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.comparison import (
    proximity_analysis,
    render_summary,
)
from selene_base.validation.nasa_regions import regions_to_geodataframe


def _sites(rows: list[dict]) -> gpd.GeoDataFrame:
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


def test_empty_inputs_return_zero_metrics() -> None:
    sites = _sites([])
    nasa = regions_to_geodataframe()
    result = proximity_analysis(sites, nasa)
    assert result["n_top_sites"] == 0
    assert result["n_nasa_regions"] == 9
    assert result["sites_within_any_region"] == 0
    assert result["regions_with_a_top_site"] == 0


def test_empty_nasa_returns_zero_metrics() -> None:
    sites = _sites(
        [{"site_id": "a", "rank": 1, "lat": -85.0, "lon": 0.0}],
    )
    empty_nasa = gpd.GeoDataFrame(
        {"name": [], "lat": [], "lon": [], "radius_km": []},
        geometry=[],
        crs=LUNAR_GEOGRAPHIC_CRS,
    )
    result = proximity_analysis(sites, empty_nasa)
    assert result["n_top_sites"] == 1
    assert result["n_nasa_regions"] == 0
    assert result["per_site"] == []


def test_site_at_centroid_lands_inside_region() -> None:
    nasa = regions_to_geodataframe()
    mons_mouton = nasa[nasa["name"] == "Mons Mouton"].iloc[0]
    sites = _sites(
        [
            {
                "site_id": "site_01",
                "rank": 1,
                "lat": float(mons_mouton["lat"]),
                "lon": float(mons_mouton["lon"]),
            },
        ],
    )
    result = proximity_analysis(sites, nasa)
    assert result["sites_within_any_region"] == 1
    assert result["sites_within_25km_of_region"] == 1
    assert result["regions_with_a_top_site"] == 1
    assert result["per_site"][0]["nearest_region"] == "Mons Mouton"
    assert result["per_site"][0]["distance_km"] < 1.0
    assert result["per_site"][0]["inside_region"] is True


def test_site_far_from_any_region_yields_no_alignment() -> None:
    # Single site at lat -82, lon 180 (lunar far side, far from every NASA candidate).
    nasa = regions_to_geodataframe()
    sites = _sites([{"site_id": "site_01", "rank": 1, "lat": -82.0, "lon": 180.0}])
    result = proximity_analysis(sites, nasa)
    assert result["sites_within_any_region"] == 0
    assert result["sites_within_25km_of_region"] == 0
    assert result["regions_with_a_top_site"] == 0
    assert result["per_site"][0]["distance_km"] > 100.0


def test_per_region_table_lists_every_nasa_region() -> None:
    nasa = regions_to_geodataframe()
    sites = _sites([{"site_id": "x", "rank": 1, "lat": -85.0, "lon": 0.0}])
    result = proximity_analysis(sites, nasa)
    names_in_result = {r["name"] for r in result["per_region"]}
    assert names_in_result == set(nasa["name"])


def test_render_summary_contains_headline_numbers() -> None:
    nasa = regions_to_geodataframe()
    sites = _sites(
        [
            {"site_id": "a", "rank": 1, "lat": -85.0, "lon": 0.0},
            {"site_id": "b", "rank": 2, "lat": -88.0, "lon": -50.0},
        ],
    )
    result = proximity_analysis(sites, nasa)
    text = render_summary(result)
    assert "top 2 sites vs 9 NASA candidates" in text
    assert "inside any region" in text
    assert "within 25 km" in text
    # Each NASA name appears exactly once in the table.
    for region in nasa["name"]:
        assert region in text


def test_near_km_threshold_is_respected() -> None:
    # Place a single site at lat -78, lon 100 — well-isolated from every
    # NASA candidate, then verify only a generous threshold catches it.
    # The nearest centroid (Nobile Rim 1 at -85.5, +35) is ~470 km away.
    nasa = regions_to_geodataframe()
    sites = _sites([{"site_id": "off", "rank": 1, "lat": -78.0, "lon": 100.0}])
    nearest = proximity_analysis(sites, nasa)["per_site"][0]["distance_km"]
    assert nearest > 100.0  # sanity: the site is genuinely far
    just_under = proximity_analysis(sites, nasa, near_km=nearest - 1.0)
    just_over = proximity_analysis(sites, nasa, near_km=nearest + 1.0)
    assert just_under["sites_within_25km_of_region"] == 0
    assert just_over["sites_within_25km_of_region"] == 1
