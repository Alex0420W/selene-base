"""Tests for :mod:`selene_base.criteria.coupling`."""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.coupling import (
    compute,
    derive_distance_to_psr,
    derive_distance_to_sunlit_ridge,
)


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims=("y", "x"))


class TestDeriveDistanceToPsr:
    def test_psr_pixel_returns_zero(self) -> None:
        # 5x5 grid; one PSR pixel in the middle.
        illum = np.full((5, 5), 0.5, dtype=np.float64)
        illum[2, 2] = 0.0
        out = derive_distance_to_psr(_da(illum), pixel_size_m=240.0).to_numpy()
        assert out[2, 2] == pytest.approx(0.0)

    def test_distance_is_pixel_units_times_metres(self) -> None:
        illum = np.full((11, 11), 0.5)
        illum[5, 5] = 0.0
        out = derive_distance_to_psr(_da(illum), pixel_size_m=1000.0).to_numpy()
        # Cardinal neighbour at pixel distance 1 -> 1000 m.
        assert out[5, 6] == pytest.approx(1000.0)
        # Diagonal neighbour at pixel distance sqrt(2) -> ~1414 m.
        assert out[4, 4] == pytest.approx(1000.0 * math.sqrt(2.0), rel=1e-6)
        # Two cardinal steps away -> 2000 m.
        assert out[5, 7] == pytest.approx(2000.0)

    def test_no_psr_returns_all_nan(self) -> None:
        illum = np.full((5, 5), 0.9)
        out = derive_distance_to_psr(_da(illum)).to_numpy()
        assert np.isnan(out).all()

    def test_threshold_validated(self) -> None:
        with pytest.raises(ValueError, match="psr_threshold"):
            derive_distance_to_psr(_da(np.zeros((3, 3))), psr_threshold=0.0)
        with pytest.raises(ValueError, match="psr_threshold"):
            derive_distance_to_psr(_da(np.zeros((3, 3))), psr_threshold=1.0)

    def test_pixel_size_validated(self) -> None:
        with pytest.raises(ValueError, match="pixel_size_m"):
            derive_distance_to_psr(_da(np.zeros((3, 3))), pixel_size_m=0.0)


class TestDeriveDistanceToSunlitRidge:
    def test_ridge_pixel_returns_zero(self) -> None:
        # 5x5 grid; one cell satisfies all three thresholds.
        illum = np.full((5, 5), 0.3)  # below 0.70 threshold
        slope = np.full((5, 5), 1.0)  # below 5 deg
        illum[2, 2] = 0.85
        slope[2, 2] = 12.0
        out = derive_distance_to_sunlit_ridge(_da(illum), _da(slope), pixel_size_m=240.0).to_numpy()
        assert out[2, 2] == pytest.approx(0.0)

    def test_low_illumination_excluded(self) -> None:
        illum = np.full((3, 3), 0.5)  # too dark
        slope = np.full((3, 3), 12.0)
        out = derive_distance_to_sunlit_ridge(_da(illum), _da(slope), pixel_size_m=240.0).to_numpy()
        assert np.isnan(out).all()

    def test_too_steep_excluded(self) -> None:
        illum = np.full((3, 3), 0.85)
        slope = np.full((3, 3), 40.0)  # cliff
        out = derive_distance_to_sunlit_ridge(_da(illum), _da(slope), pixel_size_m=240.0).to_numpy()
        assert np.isnan(out).all()

    def test_too_flat_excluded(self) -> None:
        illum = np.full((3, 3), 0.85)
        slope = np.full((3, 3), 1.0)  # plain
        out = derive_distance_to_sunlit_ridge(_da(illum), _da(slope), pixel_size_m=240.0).to_numpy()
        assert np.isnan(out).all()

    def test_distance_is_metric(self) -> None:
        illum = np.full((11, 11), 0.85)
        slope = np.full((11, 11), 30.0)  # too steep everywhere except one cell
        slope[5, 5] = 12.0
        out = derive_distance_to_sunlit_ridge(_da(illum), _da(slope), pixel_size_m=500.0).to_numpy()
        assert out[5, 5] == pytest.approx(0.0)
        assert out[5, 6] == pytest.approx(500.0)
        assert out[4, 4] == pytest.approx(500.0 * math.sqrt(2.0), rel=1e-6)

    def test_shape_mismatch_rejected(self) -> None:
        with pytest.raises(ValueError, match="shape"):
            derive_distance_to_sunlit_ridge(
                _da(np.zeros((3, 3))),
                _da(np.zeros((4, 4))),
            )

    @pytest.mark.parametrize(
        "kwarg,bad",
        [
            ("illumination_threshold", 1.5),
            ("illumination_threshold", -0.1),
            ("slope_min_deg", -1.0),
            ("pixel_size_m", 0.0),
        ],
    )
    def test_param_validation(self, kwarg: str, bad: float) -> None:
        with pytest.raises(ValueError):
            derive_distance_to_sunlit_ridge(
                _da(np.full((3, 3), 0.85)),
                _da(np.full((3, 3), 12.0)),
                **{kwarg: bad},
            )


class TestCompute:
    def test_zero_distances_score_one(self) -> None:
        psr = _da(np.zeros((3, 3)))
        ridge = _da(np.zeros((3, 3)))
        out = compute(psr, ridge, coupling_distance_km=5.0).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_psr_zero_ridge_far_scores_zero(self) -> None:
        # Cell at PSR (distance 0) but 10 km from ridge with cap = 5 km.
        psr = _da(np.zeros((1, 1)))
        ridge = _da(np.full((1, 1), 10_000.0))
        out = compute(psr, ridge, coupling_distance_km=5.0).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_two_km_each_with_five_km_cap(self) -> None:
        # Spec example: 2 km from each, cap 5 km -> (1-0.4)*(1-0.4) = 0.36
        psr = _da(np.full((1, 1), 2_000.0))
        ridge = _da(np.full((1, 1), 2_000.0))
        out = compute(psr, ridge, coupling_distance_km=5.0).to_numpy()
        np.testing.assert_allclose(out, 0.36, atol=1e-9)

    def test_product_is_below_arithmetic_mean(self) -> None:
        # The structural property: when both falloffs are < 1, the product
        # is strictly below their arithmetic mean. That's why a sum-based
        # aggregator misses the conjunction.
        psr = _da(np.full((1, 1), 1_000.0))  # falloff: 1 - 0.2 = 0.8
        ridge = _da(np.full((1, 1), 3_000.0))  # falloff: 1 - 0.6 = 0.4
        out = compute(psr, ridge, coupling_distance_km=5.0).to_numpy()[0, 0]
        product = 0.8 * 0.4  # = 0.32
        mean = (0.8 + 0.4) / 2.0  # = 0.60
        assert out == pytest.approx(product)
        assert product < mean

    def test_nan_propagates(self) -> None:
        psr = _da(np.array([[np.nan, 0.0]]))
        ridge = _da(np.array([[0.0, np.nan]]))
        out = compute(psr, ridge).to_numpy()
        assert math.isnan(out[0, 0])
        assert math.isnan(out[0, 1])

    def test_above_cap_clipped_to_zero(self) -> None:
        # 100 km > 5 km cap -> falloff is 0, product is 0.
        psr = _da(np.full((1, 1), 100_000.0))
        ridge = _da(np.full((1, 1), 100_000.0))
        out = compute(psr, ridge, coupling_distance_km=5.0).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_shape_mismatch_rejected(self) -> None:
        with pytest.raises(ValueError, match="shape"):
            compute(_da(np.zeros((3, 3))), _da(np.zeros((4, 4))))

    def test_non_positive_coupling_rejected(self) -> None:
        with pytest.raises(ValueError, match="coupling_distance_km"):
            compute(_da(np.zeros((1, 1))), _da(np.zeros((1, 1))), coupling_distance_km=0.0)

    def test_smaller_cap_makes_score_drop_faster(self) -> None:
        psr = _da(np.full((1, 1), 2_000.0))
        ridge = _da(np.full((1, 1), 2_000.0))
        s_5 = compute(psr, ridge, coupling_distance_km=5.0).to_numpy()[0, 0]
        s_3 = compute(psr, ridge, coupling_distance_km=3.0).to_numpy()[0, 0]
        # 2 km < both caps. Smaller cap -> lower score (steeper falloff).
        assert s_3 < s_5
