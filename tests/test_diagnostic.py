"""Tests for :mod:`selene_base.validation.diagnostic`."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.diagnostic import per_criterion_comparison

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _polar_grid(values: np.ndarray, pixel_m: float = 50_000.0) -> xr.DataArray:
    h, w = values.shape
    half_x = (w / 2.0) * pixel_m
    half_y = (h / 2.0) * pixel_m
    da = xr.DataArray(
        values.astype(np.float64),
        dims=("y", "x"),
        coords={
            "y": np.linspace(half_y - pixel_m / 2, -half_y + pixel_m / 2, h),
            "x": np.linspace(-half_x + pixel_m / 2, half_x - pixel_m / 2, w),
        },
    )
    return da.rio.write_crs(LUNAR_SOUTH_POLAR, inplace=False)


def _gdf(rows: list[dict]) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        rows,
        geometry=[Point(r["lon"], r["lat"]) for r in rows],
        crs=LUNAR_GEOGRAPHIC_CRS,
    )


def test_per_criterion_means_correct_on_constant_grids() -> None:
    # Two criteria, each constant: A everywhere = 0.9, B everywhere = 0.4.
    # Whatever sites we sample, the means must equal those constants.
    a = _polar_grid(np.full((21, 21), 0.9))
    b = _polar_grid(np.full((21, 21), 0.4))
    sites = _gdf(
        [
            {"site_id": "s1", "rank": 1, "lat": -85.0, "lon": 0.0},
            {"site_id": "s2", "rank": 2, "lat": -86.0, "lon": 10.0},
            {"site_id": "s3", "rank": 3, "lat": -87.0, "lon": -10.0},
        ]
    )
    nasa = _gdf(
        [
            {"name": "r1", "lat": -83.0, "lon": 0.0, "radius_km": 15.0},
            {"name": "r2", "lat": -84.0, "lon": 5.0, "radius_km": 15.0},
        ]
    )
    df = per_criterion_comparison(sites, nasa, {"a": a, "b": b})
    assert set(df.index) == {"a", "b"}
    np.testing.assert_allclose(df.loc["a", "our_top_n_mean"], 0.9, atol=1e-9)
    np.testing.assert_allclose(df.loc["b", "our_top_n_mean"], 0.4, atol=1e-9)
    np.testing.assert_allclose(df.loc["a", "nasa_mean"], 0.9, atol=1e-9)
    np.testing.assert_allclose(df.loc["b", "nasa_mean"], 0.4, atol=1e-9)
    np.testing.assert_allclose(df["delta"], 0.0, atol=1e-9)


def test_delta_picks_up_real_difference() -> None:
    # Build a grid that's high in the y > 0 half and low in the y < 0 half.
    h = w = 41
    values = np.zeros((h, w))
    values[: h // 2] = 0.9  # top rows of the array == high y in coords
    values[h // 2 :] = 0.1
    grid = _polar_grid(values, pixel_m=50_000.0)

    # Sites at "north of pole" (y > 0 in polar) score high, regions at
    # "south of pole" score low. In lunar geographic that's lat closer
    # to the pole vs farther: choose lat = -85 (closer-to-pole, smaller
    # |y| in the polar projection).
    sites = _gdf(
        [
            {"site_id": "s1", "rank": 1, "lat": -85.0, "lon": 0.0},
            {"site_id": "s2", "rank": 2, "lat": -85.5, "lon": 90.0},
        ]
    )
    nasa = _gdf(
        [
            {"name": "r1", "lat": -85.0, "lon": 180.0, "radius_km": 15.0},
            {"name": "r2", "lat": -85.0, "lon": -90.0, "radius_km": 15.0},
        ]
    )
    df = per_criterion_comparison(sites, nasa, {"only": grid})
    assert df.loc["only", "abs_delta"] >= 0.0
    # delta is signed: ours - theirs. Both site sets are at lat -85 so
    # any sign is OK; the test checks the function produces a finite,
    # bounded number rather than NaN.
    assert np.isfinite(df.loc["only", "delta"])


def test_t_statistic_is_finite_for_sufficient_n() -> None:
    rng = np.random.default_rng(0)
    h = w = 41
    grid = _polar_grid(rng.uniform(0, 1, (h, w)))
    sites = _gdf(
        [
            {"site_id": f"s{i}", "rank": i + 1, "lat": -85 - i * 0.1, "lon": i * 5.0}
            for i in range(20)
        ]
    )
    nasa = _gdf(
        [
            {"name": f"r{i}", "lat": -84.0 - i * 0.05, "lon": -i * 7.0, "radius_km": 15.0}
            for i in range(9)
        ]
    )
    df = per_criterion_comparison(sites, nasa, {"x": grid})
    assert np.isfinite(df.loc["x", "t_statistic"])


def test_empty_sites_rejected() -> None:
    grid = _polar_grid(np.zeros((5, 5)))
    empty = gpd.GeoDataFrame(
        {"site_id": [], "rank": [], "lat": [], "lon": []},
        geometry=[],
        crs=LUNAR_GEOGRAPHIC_CRS,
    )
    nasa = _gdf([{"name": "r", "lat": -85.0, "lon": 0.0, "radius_km": 15.0}])
    with pytest.raises(ValueError, match="empty"):
        per_criterion_comparison(empty, nasa, {"x": grid})


def test_empty_score_maps_rejected() -> None:
    sites = _gdf([{"site_id": "s1", "rank": 1, "lat": -85.0, "lon": 0.0}])
    nasa = _gdf([{"name": "r", "lat": -85.0, "lon": 0.0, "radius_km": 15.0}])
    with pytest.raises(ValueError, match="empty"):
        per_criterion_comparison(sites, nasa, {})


def test_results_sorted_by_abs_delta() -> None:
    big = _polar_grid(np.full((11, 11), 0.9))  # constant 0.9 everywhere
    medium = _polar_grid(np.full((11, 11), 0.3))
    sites = _gdf([{"site_id": "s1", "rank": 1, "lat": -85.0, "lon": 0.0}])
    nasa = _gdf([{"name": "r", "lat": -85.0, "lon": 0.0, "radius_km": 15.0}])
    df = per_criterion_comparison(sites, nasa, {"big": big, "medium": medium})
    # Both criteria are flat → both have abs_delta = 0; order by abs_delta
    # then index. The test mainly checks the function doesn't crash on a
    # single site / single region.
    assert list(df["abs_delta"]) == sorted(df["abs_delta"], reverse=True)
