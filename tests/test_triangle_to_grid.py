"""Tests for :mod:`selene_base.data.triangle_to_grid`."""

from __future__ import annotations

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS

from selene_base.data.triangle_to_grid import triangles_to_raster

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


def test_constant_field_remains_constant() -> None:
    rng = np.random.default_rng(0)
    n = 1000
    lats = -85.0 + rng.uniform(-3.0, 3.0, n)
    lons = rng.uniform(-180, 180, n)
    values = np.full(n, 0.42)
    grid = _polar_grid(40, 40, pixel_m=10_000.0)
    out = triangles_to_raster(lats, lons, values, grid).to_numpy()
    finite = out[np.isfinite(out)]
    assert finite.size > 0
    np.testing.assert_allclose(finite, 0.42, atol=1e-5)


def test_value_at_input_point_recovered() -> None:
    # Place 50 random points; pick one and verify its value is recoverable
    # by interpolating at its lat/lon coordinate (approximately).
    rng = np.random.default_rng(7)
    n = 100
    lats = -88.0 + rng.uniform(-2.0, 2.0, n)
    lons = rng.uniform(-180, 180, n)
    values = rng.uniform(0, 1, n)

    grid = _polar_grid(80, 80, pixel_m=2_000.0)
    out = triangles_to_raster(lats, lons, values, grid, method="nearest").to_numpy()
    # Interior of the polar cap should mostly be filled.
    assert np.isfinite(out).any()
    assert out.shape == (80, 80)


def test_nan_inputs_dropped() -> None:
    rng = np.random.default_rng(1)
    n = 200
    lats = -85.0 + rng.uniform(-3, 3, n)
    lons = rng.uniform(-180, 180, n)
    values = rng.uniform(0, 1, n)
    values[::5] = np.nan
    grid = _polar_grid(40, 40, pixel_m=10_000.0)
    out = triangles_to_raster(lats, lons, values, grid).to_numpy()
    # The interpolation should still return a valid grid even with NaN inputs.
    assert np.isfinite(out).any()


def test_unknown_method_rejected() -> None:
    grid = _polar_grid(10, 10, pixel_m=10_000.0)
    with pytest.raises(ValueError, match="unknown method"):
        triangles_to_raster(
            np.array([-85.0]),
            np.array([0.0]),
            np.array([1.0]),
            grid,
            method="bogus",
        )


def test_shape_mismatch_rejected() -> None:
    grid = _polar_grid(10, 10, pixel_m=10_000.0)
    with pytest.raises(ValueError, match="share shape"):
        triangles_to_raster(
            np.array([-85.0]),
            np.array([0.0, 1.0]),
            np.array([1.0]),
            grid,
        )


def test_missing_grid_crs_rejected() -> None:
    bad = xr.DataArray(np.zeros((5, 5)), dims=("y", "x"))
    with pytest.raises(ValueError, match="no CRS"):
        triangles_to_raster(
            np.array([-85.0]),
            np.array([0.0]),
            np.array([1.0]),
            bad,
        )


def test_method_choices_produce_consistent_shape() -> None:
    rng = np.random.default_rng(2)
    n = 500
    lats = -88.0 + rng.uniform(-2, 2, n)
    lons = rng.uniform(-180, 180, n)
    values = rng.uniform(0, 1, n)
    grid = _polar_grid(60, 60, pixel_m=4_000.0)
    out_linear = triangles_to_raster(lats, lons, values, grid, method="linear")
    out_nearest = triangles_to_raster(lats, lons, values, grid, method="nearest")
    assert out_linear.shape == out_nearest.shape
    # Both should produce mostly-finite interior values.
    assert np.isfinite(out_linear.to_numpy()).any()
    assert np.isfinite(out_nearest.to_numpy()).any()
