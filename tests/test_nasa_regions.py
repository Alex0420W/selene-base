"""Tests for :mod:`selene_base.validation.nasa_regions`."""

from __future__ import annotations

import numpy as np
import pytest
from pyproj import Transformer
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.nasa_regions import (
    ARTEMIS_III_CANDIDATE_REGIONS,
    DEFAULT_RADIUS_KM,
    USGS_POLYGONS_GEOJSON,
    USGS_REGION_NAMES,
    regions_polygons_to_geodataframe,
    regions_to_geodataframe,
)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def test_nine_canonical_regions() -> None:
    names = {r["name"] for r in ARTEMIS_III_CANDIDATE_REGIONS}
    expected = {
        "Cabeus B",
        "Haworth",
        "Malapert Massif",
        "Mons Mouton",
        "Mons Mouton Plateau",
        "Nobile Rim 1",
        "Nobile Rim 2",
        "de Gerlache Rim 2",
        "Slater Plain",
    }
    assert names == expected
    assert len(ARTEMIS_III_CANDIDATE_REGIONS) == 9


def test_all_regions_in_southern_hemisphere() -> None:
    for region in ARTEMIS_III_CANDIDATE_REGIONS:
        assert region["lat"] < -80, f"{region['name']} not at south pole"


def test_default_geodataframe_loads_in_geographic_crs() -> None:
    gdf = regions_to_geodataframe()
    assert len(gdf) == 9
    assert gdf.crs == LUNAR_GEOGRAPHIC_CRS
    assert set(gdf.columns) == {"name", "lat", "lon", "radius_km", "geometry"}


def test_reproject_to_polar_succeeds() -> None:
    gdf = regions_to_geodataframe(target_crs=POLAR_PROJ)
    assert str(gdf.crs.to_proj4()).strip() == POLAR_PROJ.strip() or "stere" in str(
        gdf.crs.to_proj4()
    )
    assert len(gdf) == 9


def test_disk_radius_close_to_15km_after_reprojection() -> None:
    # Each disk's projected geometry should bound at roughly the requested
    # radius. Distortion is small at the pole; allow ±5 % tolerance.
    gdf = regions_to_geodataframe(target_crs=POLAR_PROJ)
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    for _, row in gdf.iterrows():
        cx, cy = transformer.transform(row["lon"], row["lat"])
        coords = list(row.geometry.exterior.coords)
        radii_m = [np.hypot(x - cx, y - cy) for x, y in coords]
        observed_km = float(np.mean(radii_m)) / 1000.0
        assert observed_km == pytest.approx(DEFAULT_RADIUS_KM, rel=0.05)


def test_geometry_contains_centroid() -> None:
    gdf = regions_to_geodataframe()
    for _, row in gdf.iterrows():
        assert row.geometry.contains(Point(row["lon"], row["lat"])), row["name"]


# ----------------------------- USGS polygons (week 10) -----------------------------


def test_usgs_geojson_ships_with_package() -> None:
    assert USGS_POLYGONS_GEOJSON.exists(), USGS_POLYGONS_GEOJSON


def test_usgs_polygons_load_with_expected_columns_and_count() -> None:
    gdf = regions_polygons_to_geodataframe()
    assert len(gdf) == 9
    assert {"Region", "RegionCode", "Area_km2", "geometry"}.issubset(gdf.columns)


def test_usgs_polygons_canonical_names_match_constant() -> None:
    gdf = regions_polygons_to_geodataframe()
    assert set(gdf["Region"]) == set(USGS_REGION_NAMES)
    # Sanity: the USGS canon uses the prefixed names that distinguish
    # them from the legacy disk centroid names ("Cabeus B" → "Peak Near
    # Cabeus B"; "de Gerlache Rim 2" stays the same).
    assert "Peak Near Cabeus B" in set(gdf["Region"])
    assert "de Gerlache Rim 2" in set(gdf["Region"])
    assert "Mons Mouton Plateau" in set(gdf["Region"])


def test_usgs_total_area_matches_published() -> None:
    gdf = regions_polygons_to_geodataframe()
    # Sum of published Area_km2 values; Mons Mouton Plateau dominates
    # at 4452 km² and the eight other regions average ~400 km² for a
    # total around 8000 km². Allow generous tolerance because the
    # property is only intended as a sanity check on data integrity.
    total = float(gdf["Area_km2"].sum())
    assert 7000.0 < total < 9000.0
    # Mons Mouton Plateau alone is over half of the total.
    plateau = float(gdf.loc[gdf["Region"] == "Mons Mouton Plateau", "Area_km2"].iloc[0])
    assert plateau > 4000.0


def test_usgs_polygons_valid_in_native_and_polar_crs() -> None:
    gdf = regions_polygons_to_geodataframe()
    assert gdf.geometry.is_valid.all()
    polar = (
        "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 "
        "+R=1737400 +no_defs +type=crs"
    )
    gdf_p = regions_polygons_to_geodataframe(target_crs=polar)
    assert gdf_p.geometry.is_valid.all()
    # After reprojection to polar metres the geometric area should
    # roughly match the published Area_km2 (modest distortion expected
    # because USGS computes Shape_Area in a different equal-area frame).
    for _, row in gdf_p.iterrows():
        observed_km2 = float(row.geometry.area) / 1e6
        published = float(row["Area_km2"])
        # Allow ±20 % — the simplified envelopes plus projection
        # distortion give a few-percent gap for most regions and a
        # larger gap for the very-large Mons Mouton Plateau because
        # USGS computed its area in a different projection.
        assert observed_km2 == pytest.approx(published, rel=0.20), row["Region"]
