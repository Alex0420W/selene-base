"""Tests for :func:`selene_base.validation.comparison.per_region_compliance_analysis`."""

from __future__ import annotations

import math

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.comparison import (
    per_region_compliance_analysis,
    render_per_region_compliance_summary,
)


def _polygons() -> gpd.GeoDataFrame:
    """Two synthetic regions with published Area_km2."""
    return gpd.GeoDataFrame(
        {
            "Region": ["Alpha", "Bravo"],
            "RegionCode": ["AA", "BB"],
            "Area_km2": [400.0, 256.0],
        },
        geometry=[
            Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1)]),
            Polygon([(2, 2), (4, 2), (4, 4), (2, 4), (2, 2)]),
        ],
        crs=LUNAR_GEOGRAPHIC_CRS,
    )


def _sites(rows: list[dict]) -> gpd.GeoDataFrame:
    if rows:
        geoms = [Point(r["lon"], r["lat"]) for r in rows]
    else:
        geoms = []
    df = pd.DataFrame(
        rows,
        columns=[
            "site_id",
            "region_name",
            "region_code",
            "rank_in_region",
            "score",
            "lat",
            "lon",
        ],
    )
    return gpd.GeoDataFrame(df, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


def test_counts_sites_per_region() -> None:
    polys = _polygons()
    sites = _sites(
        [
            {
                "site_id": 1,
                "region_name": "Alpha",
                "region_code": "AA",
                "rank_in_region": 1,
                "score": 0.9,
                "lat": 0.0,
                "lon": 0.0,
            },
            {
                "site_id": 2,
                "region_name": "Alpha",
                "region_code": "AA",
                "rank_in_region": 2,
                "score": 0.7,
                "lat": 0.5,
                "lon": 0.5,
            },
            {
                "site_id": 3,
                "region_name": "Bravo",
                "region_code": "BB",
                "rank_in_region": 1,
                "score": 0.6,
                "lat": 3.0,
                "lon": 3.0,
            },
        ]
    )
    result = per_region_compliance_analysis(sites, polys)
    assert result["n_sites_total"] == 3
    assert result["n_regions_total"] == 2
    assert result["n_regions_with_sites"] == 2
    assert result["n_regions_with_no_compliant_cells"] == 0
    assert result["regions_with_no_compliant_cells"] == []

    alpha = next(r for r in result["per_region"] if r["name"] == "Alpha")
    assert alpha["n_sites"] == 2
    assert alpha["best_score"] == 0.9
    assert alpha["mean_score"] == 0.8
    assert alpha["best_site_id"] == 1


def test_records_regions_with_no_sites() -> None:
    polys = _polygons()
    sites = _sites(
        [
            {
                "site_id": 1,
                "region_name": "Alpha",
                "region_code": "AA",
                "rank_in_region": 1,
                "score": 0.9,
                "lat": 0.0,
                "lon": 0.0,
            }
        ]
    )
    result = per_region_compliance_analysis(sites, polys)
    assert result["n_regions_with_sites"] == 1
    assert result["n_regions_with_no_compliant_cells"] == 1
    assert result["regions_with_no_compliant_cells"] == ["Bravo"]
    bravo = next(r for r in result["per_region"] if r["name"] == "Bravo")
    assert bravo["n_sites"] == 0
    assert math.isnan(bravo["best_score"])
    assert math.isnan(bravo["mean_score"])
    assert bravo["best_site_id"] == 0


def test_eligible_area_fraction_uses_published_when_no_cell_total_given() -> None:
    polys = _polygons()
    sites = _sites([])
    result = per_region_compliance_analysis(
        sites,
        polys,
        eligible_area_km2={"Alpha": 40.0, "Bravo": 12.8},
    )
    alpha = next(r for r in result["per_region"] if r["name"] == "Alpha")
    bravo = next(r for r in result["per_region"] if r["name"] == "Bravo")
    assert alpha["eligible_area_fraction"] == 40.0 / 400.0
    assert bravo["eligible_area_fraction"] == 12.8 / 256.0


def test_eligible_area_fraction_zero_when_no_polygon_area() -> None:
    polys = _polygons().drop(columns="Area_km2")
    sites = _sites([])
    result = per_region_compliance_analysis(sites, polys)
    for entry in result["per_region"]:
        assert entry["eligible_area_fraction"] == 0.0


def test_render_summary_contains_each_region_and_totals() -> None:
    polys = _polygons()
    sites = _sites(
        [
            {
                "site_id": 1,
                "region_name": "Alpha",
                "region_code": "AA",
                "rank_in_region": 1,
                "score": 0.5,
                "lat": 0.0,
                "lon": 0.0,
            }
        ]
    )
    result = per_region_compliance_analysis(
        sites, polys, eligible_area_km2={"Alpha": 40.0, "Bravo": 0.0}
    )
    text = render_per_region_compliance_summary(result)
    assert "1 total sites across 1 / 2 regions" in text
    assert "Alpha" in text
    assert "Bravo" in text
    assert "Bravo" in text  # no-compliant region must still be listed
    assert "regions with zero HLS-compliant cells: Bravo" in text


def test_empty_inputs() -> None:
    polys = _polygons()
    sites = _sites([])
    result = per_region_compliance_analysis(sites, polys)
    assert result["n_sites_total"] == 0
    assert result["n_regions_with_sites"] == 0
    assert result["n_regions_with_no_compliant_cells"] == 2
    assert sorted(result["regions_with_no_compliant_cells"]) == ["Alpha", "Bravo"]
