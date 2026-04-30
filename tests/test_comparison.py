"""Tests for :mod:`selene_base.validation.comparison`."""

from __future__ import annotations

import geopandas as gpd
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.comparison import (
    proximity_analysis,
    render_summary,
)
from selene_base.validation.nasa_regions import (
    regions_polygons_to_geodataframe,
    regions_to_geodataframe,
)


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
    # week 8 polygon metrics
    assert result["sites_inside_any_region"] == 1
    assert result["regions_containing_top_site"] == 1
    assert result["regions_with_top_site_within_disk_radius"] >= 1
    site = result["per_site"][0]
    assert site["inside_any_region"] is True
    # Centroid of a 15 km disk is roughly the disk centre, ~15 km from
    # the boundary; signed distance must be negative (inside) and around
    # -15 km.
    assert site["distance_to_edge_km"] < -10.0
    assert site["distance_to_edge_km"] > -16.0


def test_site_far_from_any_region_yields_no_alignment() -> None:
    # Single site at lat -82, lon 180 (lunar far side, far from every NASA candidate).
    nasa = regions_to_geodataframe()
    sites = _sites([{"site_id": "site_01", "rank": 1, "lat": -82.0, "lon": 180.0}])
    result = proximity_analysis(sites, nasa)
    assert result["sites_within_any_region"] == 0
    assert result["sites_within_25km_of_region"] == 0
    assert result["regions_with_a_top_site"] == 0
    assert result["per_site"][0]["distance_km"] > 100.0
    # week 8 polygon metrics
    assert result["sites_inside_any_region"] == 0
    assert result["regions_containing_top_site"] == 0
    assert result["regions_with_top_site_within_disk_radius"] == 0
    site = result["per_site"][0]
    assert site["inside_any_region"] is False
    assert site["distance_to_edge_km"] > 0  # outside every disk
    # signed distance to nearest edge ~= centroid distance - disk radius
    assert site["distance_to_edge_km"] >= site["distance_km"] - 16.0


def test_polygon_inside_distinct_from_centroid_proximity() -> None:
    """A site within 25 km of a centroid but outside the 15 km disk should
    register on the legacy centroid metric but NOT on the polygon metric.
    """
    nasa = regions_to_geodataframe()
    haworth = nasa[nasa["name"] == "Haworth"].iloc[0]
    # Bump latitude by ~0.7 deg ≈ 21 km — outside the 15 km disk but
    # within 25 km of the centroid.
    sites = _sites(
        [
            {
                "site_id": "off",
                "rank": 1,
                "lat": float(haworth["lat"]) + 0.7,
                "lon": float(haworth["lon"]),
            }
        ]
    )
    result = proximity_analysis(sites, nasa, near_km=25.0)
    # Centroid-distance metric counts it (≤ 25 km of Haworth).
    assert result["sites_within_25km_of_region"] == 1
    # Polygon-inside metric does NOT (outside the 15 km disk).
    assert result["sites_inside_any_region"] == 0
    site = result["per_site"][0]
    assert site["inside_any_region"] is False
    # Edge distance should be ~6 km (centroid distance ~21 - radius 15).
    assert 0.0 < site["distance_to_edge_km"] < 12.0


def test_render_summary_prints_both_tables() -> None:
    nasa = regions_to_geodataframe()
    sites = _sites([{"site_id": "x", "rank": 1, "lat": -85.0, "lon": 0.0}])
    text = render_summary(proximity_analysis(sites, nasa))
    assert "centroid-distance metrics" in text
    assert "15 km disk metrics" in text
    assert "inside any 15 km disk" in text
    assert "regions containing a top site" in text
    assert "dist-edge" in text


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
    assert "inside any 15 km disk" in text
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


# ----------------------- USGS polygon pathway (week 10) -----------------------


def test_usgs_metrics_absent_when_polygons_not_provided() -> None:
    nasa = regions_to_geodataframe()
    sites = _sites([{"site_id": "a", "rank": 1, "lat": -85.0, "lon": 0.0}])
    result = proximity_analysis(sites, nasa)
    # Legacy keys present:
    assert "sites_inside_any_region" in result
    # USGS keys absent when polygons not given:
    assert "sites_inside_any_usgs_polygon" not in result
    assert "per_usgs_region" not in result
    assert "per_site_usgs" not in result


def test_usgs_metrics_present_when_polygons_provided() -> None:
    nasa = regions_to_geodataframe()
    nasa_polys = regions_polygons_to_geodataframe()
    sites = _sites([{"site_id": "a", "rank": 1, "lat": -85.0, "lon": 0.0}])
    result = proximity_analysis(sites, nasa, nasa_regions_polygons=nasa_polys)
    assert result["n_usgs_regions"] == 9
    assert "sites_inside_any_usgs_polygon" in result
    assert "regions_with_top_site_inside_usgs_polygon" in result
    assert "median_distance_to_nearest_usgs_polygon_km" in result
    assert len(result["per_usgs_region"]) == 9
    assert len(result["per_site_usgs"]) == 1


def test_site_inside_usgs_polygon_is_detected() -> None:
    """Place a top site at the centre of Mons Mouton Plateau's USGS
    polygon (computed from the polygon's own centroid) and verify the
    USGS pathway flags it as inside.
    """
    nasa = regions_to_geodataframe()
    nasa_polys = regions_polygons_to_geodataframe()
    plateau = nasa_polys[nasa_polys["Region"] == "Mons Mouton Plateau"].iloc[0]
    centroid = plateau.geometry.centroid
    sites = _sites(
        [{"site_id": "inside_mp", "rank": 1, "lat": float(centroid.y), "lon": float(centroid.x)}]
    )
    result = proximity_analysis(sites, nasa, nasa_regions_polygons=nasa_polys)
    assert result["sites_inside_any_usgs_polygon"] == 1
    assert result["regions_with_top_site_inside_usgs_polygon"] == 1
    site = result["per_site_usgs"][0]
    assert site["inside_any_usgs_polygon"] is True
    assert site["containing_polygon_name"] == "Mons Mouton Plateau"
    assert site["distance_to_nearest_polygon_km"] == 0.0


def test_site_far_from_usgs_polygons_records_distance() -> None:
    nasa = regions_to_geodataframe()
    nasa_polys = regions_polygons_to_geodataframe()
    # Place a site at lat -82, lon 180 — opposite hemisphere from every
    # USGS polygon, so distance must be large and inside flags False.
    sites = _sites([{"site_id": "off", "rank": 1, "lat": -82.0, "lon": 180.0}])
    result = proximity_analysis(sites, nasa, nasa_regions_polygons=nasa_polys)
    assert result["sites_inside_any_usgs_polygon"] == 0
    assert result["regions_with_top_site_inside_usgs_polygon"] == 0
    site = result["per_site_usgs"][0]
    assert site["inside_any_usgs_polygon"] is False
    assert site["containing_polygon_name"] is None
    assert site["distance_to_nearest_polygon_km"] > 100.0


def test_per_usgs_region_table_lists_every_region() -> None:
    nasa = regions_to_geodataframe()
    nasa_polys = regions_polygons_to_geodataframe()
    sites = _sites([{"site_id": "a", "rank": 1, "lat": -85.0, "lon": 0.0}])
    result = proximity_analysis(sites, nasa, nasa_regions_polygons=nasa_polys)
    names = {row["name"] for row in result["per_usgs_region"]}
    assert names == set(nasa_polys["Region"])
    # Each row carries a code and area.
    for row in result["per_usgs_region"]:
        assert isinstance(row["code"], str) and row["code"]
        assert row["area_km2"] > 0


def test_render_summary_prints_three_tables() -> None:
    nasa = regions_to_geodataframe()
    nasa_polys = regions_polygons_to_geodataframe()
    sites = _sites([{"site_id": "x", "rank": 1, "lat": -85.0, "lon": 0.0}])
    text = render_summary(proximity_analysis(sites, nasa, nasa_regions_polygons=nasa_polys))
    assert "centroid-distance metrics" in text
    assert "15 km disk metrics" in text
    assert "USGS polygon metrics" in text
    assert "inside any USGS polygon" in text
    assert "USGS polygon per-region" in text


def test_empty_sites_with_usgs_polygons_returns_zeroed_keys() -> None:
    sites = _sites([])
    nasa = regions_to_geodataframe()
    nasa_polys = regions_polygons_to_geodataframe()
    result = proximity_analysis(sites, nasa, nasa_regions_polygons=nasa_polys)
    assert result["n_usgs_regions"] == 9
    assert result["sites_inside_any_usgs_polygon"] == 0
    assert result["regions_with_top_site_inside_usgs_polygon"] == 0
    assert result["per_usgs_region"] == []
    assert result["per_site_usgs"] == []
