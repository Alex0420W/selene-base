"""Tests for :mod:`selene_base.scoring.normalize`.

These functions are the only real implementations in the v0 scaffold,
so the tests need to actually exercise their behaviour: the basic
mapping, the edge cases (constant input, all-NaN, threshold input
validation), and NaN preservation.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from selene_base.scoring.normalize import inverse_threshold, min_max, optimal_range


class TestMinMax:
    def test_basic_linear_rescale(self) -> None:
        out = min_max(np.array([0.0, 5.0, 10.0]))
        np.testing.assert_allclose(out, [0.0, 0.5, 1.0])

    def test_negative_and_positive_range(self) -> None:
        out = min_max(np.array([-2.0, 0.0, 2.0]))
        np.testing.assert_allclose(out, [0.0, 0.5, 1.0])

    def test_all_zeros_returns_zeros(self) -> None:
        out = min_max(np.zeros(5))
        np.testing.assert_array_equal(out, np.zeros(5))

    def test_constant_array_returns_zeros(self) -> None:
        out = min_max(np.full(4, 7.5))
        np.testing.assert_array_equal(out, np.zeros(4))

    def test_nan_is_preserved(self) -> None:
        out = min_max(np.array([0.0, np.nan, 10.0]))
        assert math.isnan(out[1])
        assert out[0] == 0.0
        assert out[2] == 1.0

    def test_all_nan_input(self) -> None:
        out = min_max(np.array([np.nan, np.nan]))
        assert np.all(np.isnan(out))

    def test_accepts_list_input(self) -> None:
        out = min_max([0.0, 1.0, 2.0])
        np.testing.assert_allclose(out, [0.0, 0.5, 1.0])

    def test_2d_array_preserves_shape(self) -> None:
        out = min_max(np.array([[0.0, 2.0], [4.0, 6.0]]))
        assert out.shape == (2, 2)
        np.testing.assert_allclose(out, [[0.0, 1 / 3], [2 / 3, 1.0]])


class TestOptimalRange:
    def test_peak_at_target(self) -> None:
        out = optimal_range(np.array([5.0]), target=5.0, sigma=1.0)
        np.testing.assert_allclose(out, [1.0])

    def test_symmetric_decay(self) -> None:
        out = optimal_range(np.array([3.0, 5.0, 7.0]), target=5.0, sigma=1.0)
        assert out[0] == pytest.approx(out[2])
        assert out[1] > out[0]

    def test_far_from_target_approaches_zero(self) -> None:
        out = optimal_range(np.array([100.0]), target=0.0, sigma=1.0)
        assert out[0] < 1e-10

    def test_output_in_unit_range(self) -> None:
        x = np.linspace(-10, 10, 100)
        out = optimal_range(x, target=2.0, sigma=3.0)
        assert out.min() >= 0.0
        assert out.max() <= 1.0

    def test_nan_propagates(self) -> None:
        out = optimal_range(np.array([np.nan, 0.0]), target=0.0, sigma=1.0)
        assert math.isnan(out[0])
        assert out[1] == pytest.approx(1.0)

    def test_zero_sigma_rejected(self) -> None:
        with pytest.raises(ValueError, match="sigma must be"):
            optimal_range(np.array([0.0]), target=0.0, sigma=0.0)

    def test_negative_sigma_rejected(self) -> None:
        with pytest.raises(ValueError, match="sigma must be"):
            optimal_range(np.array([0.0]), target=0.0, sigma=-1.0)


class TestInverseThreshold:
    def test_zero_maps_to_one(self) -> None:
        out = inverse_threshold(np.array([0.0]), threshold=10.0)
        np.testing.assert_allclose(out, [1.0])

    def test_threshold_maps_to_zero(self) -> None:
        out = inverse_threshold(np.array([10.0]), threshold=10.0)
        np.testing.assert_allclose(out, [0.0])

    def test_linear_in_between(self) -> None:
        out = inverse_threshold(np.array([0.0, 2.5, 5.0, 7.5, 10.0]), threshold=10.0)
        np.testing.assert_allclose(out, [1.0, 0.75, 0.5, 0.25, 0.0])

    def test_above_threshold_clipped_to_zero(self) -> None:
        out = inverse_threshold(np.array([15.0, 100.0]), threshold=10.0)
        np.testing.assert_array_equal(out, [0.0, 0.0])

    def test_negative_input_clipped_to_one(self) -> None:
        out = inverse_threshold(np.array([-1.0, -100.0]), threshold=10.0)
        np.testing.assert_array_equal(out, [1.0, 1.0])

    def test_nan_is_preserved(self) -> None:
        out = inverse_threshold(np.array([np.nan, 0.0, 5.0]), threshold=10.0)
        assert math.isnan(out[0])
        assert out[1] == pytest.approx(1.0)
        assert out[2] == pytest.approx(0.5)

    def test_zero_threshold_rejected(self) -> None:
        with pytest.raises(ValueError, match="threshold must be"):
            inverse_threshold(np.array([0.0]), threshold=0.0)

    def test_negative_threshold_rejected(self) -> None:
        with pytest.raises(ValueError, match="threshold must be"):
            inverse_threshold(np.array([0.0]), threshold=-1.0)
