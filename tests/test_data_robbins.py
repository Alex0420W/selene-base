"""Smoke test for the Robbins south-polar crater catalog.

Skipped automatically when the filtered CSV is not present, so CI stays
green without the data files.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS, load_crater_catalog

ROBBINS_PATH = Path("data/raw/robbins/robbins_southpole.csv.gz")

pytestmark = pytest.mark.skipif(
    not ROBBINS_PATH.exists(),
    reason=f"data not downloaded: {ROBBINS_PATH}",
)


def test_returns_geodataframe() -> None:
    gdf = load_crater_catalog(ROBBINS_PATH)
    assert isinstance(gdf, gpd.GeoDataFrame)


def test_non_empty_and_expected_columns() -> None:
    gdf = load_crater_catalog(ROBBINS_PATH)
    assert len(gdf) > 0
    for col in ("lat", "lon", "diam_km"):
        assert col in gdf.columns


def test_filtered_to_south_polar() -> None:
    gdf = load_crater_catalog(ROBBINS_PATH)
    assert (gdf["lat"] <= -75.0).all()


def test_diameters_are_positive_finite() -> None:
    gdf = load_crater_catalog(ROBBINS_PATH)
    diams = gdf["diam_km"].to_numpy()
    assert np.isfinite(diams).any()
    assert (diams[np.isfinite(diams)] > 0).all()


def test_crs_is_lunar_geographic() -> None:
    gdf = load_crater_catalog(ROBBINS_PATH)
    assert gdf.crs == LUNAR_GEOGRAPHIC_CRS
