"""Tests for :mod:`selene_base.criteria.los_to_earth` (week 9)."""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.los_to_earth import (
    DEFAULT_MOON_RADIUS_M,
    compute,
    compute_earth_visibility_fraction,
    derive_horizon_profile,
)


def _elev(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(
        np.asarray(values, dtype=np.float32),
        dims=("y", "x"),
        name="elevation_m",
    )


class TestDeriveHorizonProfile:
    def test_flat_surface_only_curvature(self) -> None:
        # On a perfectly flat surface, the apparent horizon angle at
        # distance d is -atan(d^2 / (2R) / d) = -atan(d / (2R)).
        # The maximum over a log-spaced ray from 240 m to 100 km is
        # the *least negative* of those values, which is at the
        # smallest distance: -atan(120 / 1737400) ≈ -0.0040°.
        # The test verifies the shape of the result (negative
        # everywhere, very close to the curvature dip at one of the
        # sampled distances).
        elev = _elev([[0.0] * 9 for _ in range(9)])
        profile = derive_horizon_profile(
            elev,
            n_azimuths=4,
            max_horizon_km=10.0,  # smaller to speed the test
            pixel_size_m=240.0,
            n_distance_samples=20,
        )
        # The center pixel has rays of length up to ~5 pixels before
        # they leave the grid. The off-grid contribution is NaN, so
        # only on-grid distances contribute. Either way, the horizon
        # angle must be negative (curvature dip) and shallow.
        center = profile.values[:, 4, 4]
        assert (center < 0.0).all()
        assert (center > -1.0).all()  # curvature dip at <=10 km is well under 1°

    def test_flat_surface_curvature_at_full_range(self) -> None:
        # With a 100 km cutoff and a uniformly-zero elevation grid, the
        # most-negative horizon angle from a centre-of-mass pixel is
        # -atan(d_max / (2R)) when the ray reaches d_max. Since rays
        # extend off our small grid, we use a large grid here.
        elev = _elev([[0.0] * 1001 for _ in range(1001)])
        profile = derive_horizon_profile(
            elev,
            n_azimuths=4,
            max_horizon_km=100.0,
            pixel_size_m=240.0,
            n_distance_samples=20,
        )
        center = profile.values[:, 500, 500]
        # All four azimuths should agree by symmetry; minimum-most-negative
        # value should be close to -atan(d/(2R)) at smallest sampled
        # distance (240 m): -atan(120/1737400) ≈ -0.004°.
        assert center.max() == pytest.approx(
            np.degrees(np.arctan(-120.0 / DEFAULT_MOON_RADIUS_M)), abs=1e-3
        )
        assert (center < 0.0).all()

    def test_known_peak_dominates_one_azimuth(self) -> None:
        # 9x9 grid with a 200 m peak two pixels EAST of the centre.
        # From the centre, the +x_grid azimuth (= π/2 = 90°) should see
        # a high horizon angle; the opposite azimuth (270°) should not.
        grid = [[0.0] * 9 for _ in range(9)]
        grid[4][6] = 200.0
        elev = _elev(grid)

        profile = derive_horizon_profile(
            elev,
            n_azimuths=4,  # azimuths 0, 90, 180, 270
            max_horizon_km=2.0,
            pixel_size_m=240.0,
            n_distance_samples=20,
        )

        # Azimuth bucket 1 (90° = +x_grid) crosses the peak.
        # Azimuth bucket 3 (270° = -x_grid) does not.
        center = profile.values[:, 4, 4]
        # Peak is 2 pixels = 480 m east, height 200 m. Apparent angle:
        # atan(200 / 480) ~= 22.6°, minus a tiny curvature drop.
        assert center[1] > 20.0  # toward the peak
        assert center[3] < -0.001  # opposite direction — only curvature

    def test_param_validation(self) -> None:
        elev = _elev([[0.0]])
        with pytest.raises(ValueError, match="n_azimuths"):
            derive_horizon_profile(elev, n_azimuths=0)
        with pytest.raises(ValueError, match="max_horizon_km"):
            derive_horizon_profile(elev, max_horizon_km=0.0)
        with pytest.raises(ValueError, match="pixel_size_m"):
            derive_horizon_profile(elev, pixel_size_m=0.0)
        with pytest.raises(ValueError, match="moon_radius_m"):
            derive_horizon_profile(elev, moon_radius_m=0.0)

    def test_rejects_non_2d(self) -> None:
        with pytest.raises(ValueError, match="must be 2D"):
            derive_horizon_profile(
                xr.DataArray(np.zeros((2, 2, 2), dtype=np.float32), dims=("a", "b", "c"))
            )


class TestComputeEarthVisibilityFraction:
    def _coords_2x2(
        self, lat: float, lon: float
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        lat_da = xr.DataArray(
            np.full((2, 2), lat, dtype=np.float64),
            dims=("y", "x"),
        )
        lon_da = xr.DataArray(
            np.full((2, 2), lon, dtype=np.float64),
            dims=("y", "x"),
        )
        gamma = xr.DataArray(
            np.zeros((2, 2), dtype=np.float64),
            dims=("y", "x"),
        )
        return lat_da, lon_da, gamma

    def test_high_horizon_blocks_earth(self) -> None:
        # Pixel at (lat=-89, lon=0). Earth's elevation at peak libration
        # is roughly 6.5°; at the unfavorable extreme it's roughly -6.5°.
        # If the local horizon is 30° everywhere, Earth never visible.
        lat_da, lon_da, gamma = self._coords_2x2(-89.0, 0.0)
        horizon = xr.DataArray(
            np.full((36, 2, 2), 30.0, dtype=np.float32),
            dims=("azimuth", "y", "x"),
        )
        frac = compute_earth_visibility_fraction(
            horizon, lat_da, lon_da, gamma, n_libration_samples=24
        )
        assert (frac.to_numpy() == 0.0).all()

    def test_low_horizon_admits_earth(self) -> None:
        # Same pixel, but horizon flat at -10° (below local horizontal).
        # Earth oscillates between ~+7° and ~-6° (half cycle visible
        # when it's above 0°, all of cycle visible when horizon is at
        # -10°). So visibility should be 1.0.
        lat_da, lon_da, gamma = self._coords_2x2(-89.0, 0.0)
        horizon = xr.DataArray(
            np.full((36, 2, 2), -10.0, dtype=np.float32),
            dims=("azimuth", "y", "x"),
        )
        frac = compute_earth_visibility_fraction(
            horizon, lat_da, lon_da, gamma, n_libration_samples=24
        )
        assert (frac.to_numpy() == 1.0).all()

    def test_libration_phase_ordering(self) -> None:
        # At lat=-89, lon=0, with a flat horizon at 0°, Earth is above
        # the horizon when sub-Earth lat is negative (the favorable
        # half of the libration cycle). This should give ~50% visibility.
        # Roughness from quantising azimuth + libration is tolerated.
        lat_da, lon_da, gamma = self._coords_2x2(-89.0, 0.0)
        horizon = xr.DataArray(
            np.zeros((36, 2, 2), dtype=np.float32),
            dims=("azimuth", "y", "x"),
        )
        frac = compute_earth_visibility_fraction(
            horizon, lat_da, lon_da, gamma, n_libration_samples=200
        )
        v = float(frac.to_numpy()[0, 0])
        assert 0.4 < v < 0.6

    def test_longitude_asymmetry(self) -> None:
        # At lat=-89: the libration ellipse is centred on (0, 0). For
        # lon=0 (Earth-facing meridian) Earth's elevation reaches
        # ~+7.5° at peak; for lon=180 (anti-Earth meridian) it reaches
        # only ~+5.5°. So with a flat horizon at +6°, lon=0 sees Earth
        # for some libration samples; lon=180 sees Earth for none.
        lat_da_a, lon_da_a, gamma_a = self._coords_2x2(-89.0, 0.0)
        lat_da_b, lon_da_b, gamma_b = self._coords_2x2(-89.0, 180.0)
        horizon = xr.DataArray(
            np.full((36, 2, 2), 6.0, dtype=np.float32),
            dims=("azimuth", "y", "x"),
        )
        frac_a = compute_earth_visibility_fraction(
            horizon, lat_da_a, lon_da_a, gamma_a, n_libration_samples=200
        )
        frac_b = compute_earth_visibility_fraction(
            horizon, lat_da_b, lon_da_b, gamma_b, n_libration_samples=200
        )
        assert float(frac_a.to_numpy()[0, 0]) > 0.0
        assert float(frac_b.to_numpy()[0, 0]) == 0.0

    def test_param_validation(self) -> None:
        lat_da, lon_da, gamma = self._coords_2x2(-89.0, 0.0)
        horizon = xr.DataArray(
            np.zeros((4, 2, 2), dtype=np.float32),
            dims=("azimuth", "y", "x"),
        )
        with pytest.raises(ValueError, match="n_libration_samples"):
            compute_earth_visibility_fraction(horizon, lat_da, lon_da, gamma, n_libration_samples=0)
        with pytest.raises(ValueError, match="libration semi-axes"):
            compute_earth_visibility_fraction(
                horizon,
                lat_da,
                lon_da,
                gamma,
                earth_max_libration_lat_deg=0.0,
            )

    def test_shape_mismatch(self) -> None:
        lat_da = xr.DataArray(np.zeros((3, 3), dtype=np.float64), dims=("y", "x"))
        lon_da = xr.DataArray(np.zeros((2, 2), dtype=np.float64), dims=("y", "x"))
        gamma = xr.DataArray(np.zeros((2, 2), dtype=np.float64), dims=("y", "x"))
        horizon = xr.DataArray(
            np.zeros((4, 2, 2), dtype=np.float32),
            dims=("azimuth", "y", "x"),
        )
        with pytest.raises(ValueError, match="shape mismatch"):
            compute_earth_visibility_fraction(horizon, lat_da, lon_da, gamma)


class TestCompute:
    def _frac(self, values: list[list[float]]) -> xr.DataArray:
        return xr.DataArray(
            np.asarray(values, dtype=np.float64),
            dims=("y", "x"),
            name="los_visibility_fraction",
        )

    def test_below_min_scores_zero(self) -> None:
        out = compute(self._frac([[0.0, 0.1, 0.19999]])).to_numpy()
        np.testing.assert_allclose(out, 0.0, atol=1e-9)

    def test_at_min_scores_zero(self) -> None:
        out = compute(self._frac([[0.20]])).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_at_target_scores_one(self) -> None:
        out = compute(self._frac([[0.50]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_above_target_scores_one(self) -> None:
        out = compute(self._frac([[0.51, 0.75, 1.0]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_midpoint_scores_half(self) -> None:
        # min=0.2, target=0.5; midpoint = 0.35 -> 0.5
        out = compute(self._frac([[0.35]])).to_numpy()[0, 0]
        assert out == pytest.approx(0.5)

    def test_nan_propagates(self) -> None:
        out = compute(self._frac([[np.nan, 0.5]])).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    def test_custom_thresholds(self) -> None:
        out = compute(
            self._frac([[0.0, 0.1, 0.5, 0.9, 1.0]]),
            min_visibility=0.1,
            target_visibility=0.9,
        ).to_numpy()[0]
        assert out[0] == pytest.approx(0.0)
        assert out[1] == pytest.approx(0.0)
        assert out[2] == pytest.approx(0.5)
        assert out[3] == pytest.approx(1.0)
        assert out[4] == pytest.approx(1.0)

    def test_invalid_thresholds(self) -> None:
        with pytest.raises(ValueError, match="min_visibility"):
            compute(self._frac([[0.5]]), min_visibility=0.6, target_visibility=0.5)
        with pytest.raises(ValueError, match="min_visibility"):
            compute(self._frac([[0.5]]), min_visibility=-0.1, target_visibility=0.5)
        with pytest.raises(ValueError, match="min_visibility"):
            compute(self._frac([[0.5]]), min_visibility=0.2, target_visibility=1.5)
