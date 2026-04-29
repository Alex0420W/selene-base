"""Tests for :mod:`selene_base.criteria.seismic`."""

from __future__ import annotations

import math

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import LineString, Point

from selene_base.criteria.seismic import compute, distance_to_scarps

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _polar_grid(width: int, height: int, pixel_m: float) -> xr.DataArray:
    half_x = (width / 2.0) * pixel_m
    half_y = (height / 2.0) * pixel_m
    da = xr.DataArray(
        np.zeros((height, width), dtype=np.float32),
        dims=("y", "x"),
        coords={
            "y": np.linspace(half_y - pixel_m / 2, -half_y + pixel_m / 2, height),
            "x": np.linspace(-half_x + pixel_m / 2, half_x - pixel_m / 2, width),
        },
    )
    return da.rio.write_crs(LUNAR_SOUTH_POLAR, inplace=False)


class TestDistanceToScarps:
    def test_empty_catalog_returns_inf(self) -> None:
        grid = _polar_grid(5, 5, pixel_m=1000.0)
        empty = gpd.GeoDataFrame({"id": []}, geometry=[], crs=LUNAR_SOUTH_POLAR)
        out = distance_to_scarps(empty, grid).to_numpy()
        assert np.all(np.isinf(out))

    def test_point_at_origin(self) -> None:
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        scarps = gpd.GeoDataFrame({"id": [0]}, geometry=[Point(0.0, 0.0)], crs=LUNAR_SOUTH_POLAR)
        out = distance_to_scarps(scarps, grid).to_numpy()
        # Centre pixel (5, 5) is at (0, 0) — distance ~ 0.
        assert out[5, 5] == pytest.approx(0.0, abs=0.05)
        # Diagonal corner at (~5 km, 5 km) → ~7.07 km
        diag = out[0, 0]
        assert 6.5 < diag < 8.0

    def test_line_distance(self) -> None:
        # Vertical scarp along x=0; pixel (row=5, col=8) is 3 pixels east.
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        line = LineString([(0.0, -10_000.0), (0.0, 10_000.0)])
        scarps = gpd.GeoDataFrame({"id": [0]}, geometry=[line], crs=LUNAR_SOUTH_POLAR)
        out = distance_to_scarps(scarps, grid).to_numpy()
        assert out[5, 5] == pytest.approx(0.0, abs=0.05)
        assert out[5, 8] == pytest.approx(3.0, abs=0.1)
        assert out[5, 2] == pytest.approx(3.0, abs=0.1)


class TestCompute:
    def test_far_distance_scores_one(self) -> None:
        d = xr.DataArray(np.array([[60.0, 100.0]]), dims=("y", "x"))
        out = compute(d, safe_distance_km=50.0).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_zero_distance_scores_zero(self) -> None:
        d = xr.DataArray(np.array([[0.0]]), dims=("y", "x"))
        out = compute(d, safe_distance_km=50.0).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_linear_in_between(self) -> None:
        d = xr.DataArray(np.array([[10.0, 25.0, 50.0]]), dims=("y", "x"))
        out = compute(d, safe_distance_km=50.0).to_numpy()
        np.testing.assert_allclose(out, [[0.2, 0.5, 1.0]], atol=1e-9)

    def test_inf_distance_treated_as_safe(self) -> None:
        d = xr.DataArray(np.array([[np.inf]]), dims=("y", "x"))
        out = compute(d).to_numpy()
        assert out[0, 0] == pytest.approx(1.0)

    def test_nan_propagates(self) -> None:
        d = xr.DataArray(np.array([[np.nan, 50.0]]), dims=("y", "x"))
        out = compute(d).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    def test_safe_distance_must_be_positive(self) -> None:
        d = xr.DataArray(np.array([[10.0]]), dims=("y", "x"))
        with pytest.raises(ValueError, match="safe_distance_km"):
            compute(d, safe_distance_km=0.0)
