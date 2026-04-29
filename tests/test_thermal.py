"""Tests for :mod:`selene_base.criteria.thermal` (PRP-based, week 6+).

Defaults corrected in week 8: ``target_temp_k`` 230 -> 140,
``sigma_k`` 50 -> 30. See criteria/thermal.py module docstring for the
full rationale.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.thermal import compute


def _grid(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


class TestCompute:
    def test_target_temperature_scores_one(self) -> None:
        # Default target is 140 K (week 8 correction).
        out = compute(_grid([[140.0]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_polar_data_distribution_now_in_signal_band(self) -> None:
        # The corrected defaults put the actual data distribution
        # (median ~131 K, peak ~211 K) inside the responsive band of
        # the criterion. A cell at the data median should score above
        # 0.9; a cell at the data peak should score in the tail (~0.05).
        out_median = compute(_grid([[131.0]])).to_numpy()[0, 0]
        out_peak = compute(_grid([[211.0]])).to_numpy()[0, 0]
        assert out_median > 0.9
        assert 0.0 < out_peak < 0.15

    def test_thirty_kelvin_offset_scores_around_six_tenths(self) -> None:
        # The narrower sigma=30 K means a 30 K offset scores about
        # exp(-0.5) ~= 0.6065 — the design contrast for the criterion.
        out = compute(_grid([[170.0]])).to_numpy()[0, 0]
        assert 0.55 < out < 0.65

    def test_far_from_target_drops_to_tail(self) -> None:
        # 230 K (the OLD target) is now 90 K from the new peak; with
        # sigma=30 the score should be deep in the tail.
        out = compute(_grid([[230.0]])).to_numpy()[0, 0]
        assert out < 0.02

    def test_score_decays_with_offset(self) -> None:
        s_target = compute(_grid([[140.0]])).to_numpy()[0, 0]
        s_near = compute(_grid([[150.0]])).to_numpy()[0, 0]
        s_far = compute(_grid([[200.0]])).to_numpy()[0, 0]
        assert s_target == pytest.approx(1.0)
        assert s_near > s_far
        assert s_near < s_target

    def test_nan_propagates(self) -> None:
        out = compute(_grid([[np.nan, 140.0]])).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    @pytest.mark.parametrize("kwarg", ["target_temp_k", "sigma_k"])
    def test_non_positive_params_rejected(self, kwarg: str) -> None:
        with pytest.raises(ValueError, match=kwarg):
            compute(_grid([[140.0]]), **{kwarg: 0.0})

    def test_custom_target_overrides_default(self) -> None:
        out = compute(_grid([[200.0]]), target_temp_k=200.0).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_custom_sigma_overrides_default(self) -> None:
        # With sigma=50 (the old default), a 50 K offset still scores ~0.6.
        out = compute(_grid([[190.0]]), target_temp_k=140.0, sigma_k=50.0).to_numpy()[0, 0]
        assert 0.55 < out < 0.65
