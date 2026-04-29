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
