"""Tests for :mod:`selene_base.criteria.slope`.

Synthetic DEMs with closed-form slope make the algorithmic tests
deterministic; they run in CI without any downloaded data.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.slope import compute, derive_slope_degrees


def _planar_dem(shape: tuple[int, int], pixel_size_m: float, slope_deg: float) -> xr.DataArray:
    """A perfectly tilted plane: dz/dx = tan(slope_deg), dz/dy = 0."""
    h, w = shape
    grad = math.tan(math.radians(slope_deg))
    x = np.arange(w) * pixel_size_m
    y = np.arange(h) * pixel_size_m
    z = np.broadcast_to(x * grad, (h, w)).astype(np.float64)
    return xr.DataArray(
        z,
        dims=("y", "x"),
        coords={"y": y, "x": x},
        name="elevation_m",
    )


class TestDeriveSlopeDegrees:
    @pytest.mark.parametrize("known_slope", [0.0, 5.0, 10.0, 30.0])
    def test_planar_surface_recovers_known_slope(self, known_slope: float) -> None:
        dem = _planar_dem((50, 50), pixel_size_m=240.0, slope_deg=known_slope)
        out = derive_slope_degrees(dem, pixel_size_m=240.0)
        # interior should match within numeric tolerance
        interior = out.to_numpy()[1:-1, 1:-1]
        np.testing.assert_allclose(interior, known_slope, atol=1e-6)

    def test_edges_are_nan(self) -> None:
        dem = _planar_dem((20, 20), pixel_size_m=240.0, slope_deg=10.0)
        out = derive_slope_degrees(dem, pixel_size_m=240.0).to_numpy()
        assert np.all(np.isnan(out[0, :]))
        assert np.all(np.isnan(out[-1, :]))
        assert np.all(np.isnan(out[:, 0]))
        assert np.all(np.isnan(out[:, -1]))

    def test_input_nan_propagates(self) -> None:
        dem = _planar_dem((20, 20), pixel_size_m=240.0, slope_deg=5.0)
        z = dem.to_numpy().copy()
        z[10, 10] = np.nan
        dem = xr.DataArray(z, dims=("y", "x"), coords=dem.coords)
        out = derive_slope_degrees(dem, pixel_size_m=240.0).to_numpy()
        assert math.isnan(out[10, 10])

    def test_pixel_size_must_be_positive(self) -> None:
        dem = _planar_dem((10, 10), pixel_size_m=240.0, slope_deg=0.0)
        with pytest.raises(ValueError, match="pixel_size_m"):
            derive_slope_degrees(dem, pixel_size_m=0.0)

    def test_missing_dims_rejected(self) -> None:
        bad = xr.DataArray(np.zeros((10, 10)), dims=("row", "col"))
        with pytest.raises(ValueError, match="\\('y', 'x'\\)"):
            derive_slope_degrees(bad, pixel_size_m=240.0)

    def test_diagonal_slope_recovered(self) -> None:
        # 45° tilt along the diagonal: dz/dx = dz/dy = 1, slope = arctan(sqrt(2))
        h = w = 30
        x = np.arange(w) * 100.0
        y = np.arange(h) * 100.0
        xs, ys = np.meshgrid(x, y)
        z = xs + ys  # gradient (100, 100) per metre direction means tan(45)
        dem = xr.DataArray(z, dims=("y", "x"), coords={"y": y, "x": x})
        out = derive_slope_degrees(dem, pixel_size_m=100.0).to_numpy()
        expected = math.degrees(math.atan(math.sqrt(2)))  # ~54.7°
        np.testing.assert_allclose(out[1:-1, 1:-1], expected, atol=1e-6)

    def test_pixel_size_units_matter(self) -> None:
        # Same elevation array, different pixel size -> different slope.
        h = w = 20
        z = (np.arange(w) * 1.0).reshape(1, w).repeat(h, axis=0)
        dem = xr.DataArray(z, dims=("y", "x"), coords={"y": np.arange(h), "x": np.arange(w)})
        s_at_1 = derive_slope_degrees(dem, pixel_size_m=1.0).to_numpy()[10, 10]
        s_at_10 = derive_slope_degrees(dem, pixel_size_m=10.0).to_numpy()[10, 10]
        # 10x finer dx → 10x smaller slope at the same elevation step.
        assert s_at_1 > s_at_10


class TestCompute:
    def test_zero_slope_scores_one(self) -> None:
        slope = xr.DataArray(np.zeros((5, 5)), dims=("y", "x"))
        out = compute(slope).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_threshold_slope_scores_zero(self) -> None:
        slope = xr.DataArray(np.full((3, 3), 15.0), dims=("y", "x"))
        out = compute(slope).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_above_threshold_scores_zero(self) -> None:
        slope = xr.DataArray(np.array([[30.0, 90.0], [25.0, 50.0]]), dims=("y", "x"))
        out = compute(slope).to_numpy()
        np.testing.assert_array_equal(out, np.zeros((2, 2)))

    def test_linear_ramp_in_between(self) -> None:
        slope = xr.DataArray(
            np.array([[0.0, 3.0, 7.5, 10.0, 12.0, 15.0]]),
            dims=("y", "x"),
        )
        out = compute(slope, max_slope_deg=15.0).to_numpy()
        np.testing.assert_allclose(out, [[1.0, 0.8, 0.5, 1.0 / 3.0, 0.2, 0.0]], atol=1e-9)

    def test_nan_propagates(self) -> None:
        slope = xr.DataArray(np.array([[np.nan, 0.0]]), dims=("y", "x"))
        out = compute(slope).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    def test_max_slope_deg_must_be_positive(self) -> None:
        slope = xr.DataArray(np.zeros((2, 2)), dims=("y", "x"))
        with pytest.raises(ValueError, match="max_slope_deg"):
            compute(slope, max_slope_deg=0.0)

    def test_custom_threshold_changes_score(self) -> None:
        slope = xr.DataArray(np.array([[5.0]]), dims=("y", "x"))
        s10 = compute(slope, max_slope_deg=10.0).to_numpy()[0, 0]
        s20 = compute(slope, max_slope_deg=20.0).to_numpy()[0, 0]
        assert s20 > s10  # gentler cutoff → more lenient score
