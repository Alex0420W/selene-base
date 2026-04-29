"""Tests for :mod:`selene_base.criteria.thermal`."""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.thermal import compute


def _grid(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


class TestCompute:
    def test_target_with_zero_range_scores_one(self) -> None:
        # Tmax == Tmin == 180 → mean 180, range 0 → score 1.0
        out = compute(_grid([[180.0]]), _grid([[180.0]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_far_from_target_drops_mean_score(self) -> None:
        # mean = 300 K, range 0 — mean Gaussian dominates
        out = compute(_grid([[300.0]]), _grid([[300.0]])).to_numpy()[0, 0]
        # exp(-((300-180)^2) / (2*50^2)) = exp(-2.88) ~ 0.056
        assert 0.0 <= out < 0.1

    def test_max_range_zero_score(self) -> None:
        # mean = target, but range = max_range -> score 0
        out = compute(_grid([[280.0]]), _grid([[80.0]])).to_numpy()[0, 0]
        assert out == 0.0

    def test_excess_range_clipped_to_zero(self) -> None:
        out = compute(_grid([[400.0]]), _grid([[100.0]])).to_numpy()[0, 0]
        assert out == 0.0

    def test_nan_propagates_from_either_input(self) -> None:
        out = compute(_grid([[np.nan, 180.0]]), _grid([[180.0, np.nan]])).to_numpy()
        assert math.isnan(out[0, 0])
        assert math.isnan(out[0, 1])

    def test_shape_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="shape"):
            compute(_grid([[1.0, 2.0]]), _grid([[3.0]]))

    @pytest.mark.parametrize("kwarg", ["target_temp_k", "sigma_k", "max_range_k"])
    def test_non_positive_params_rejected(self, kwarg: str) -> None:
        with pytest.raises(ValueError, match=kwarg):
            compute(_grid([[180.0]]), _grid([[180.0]]), **{kwarg: 0.0})

    def test_score_decays_with_mean_offset(self) -> None:
        s_target = compute(_grid([[180.0]]), _grid([[180.0]])).to_numpy()[0, 0]
        s_offset = compute(_grid([[230.0]]), _grid([[230.0]])).to_numpy()[0, 0]
        assert s_target == pytest.approx(1.0)
        assert s_offset < s_target

    def test_score_decays_with_range(self) -> None:
        # Same mean (180 K), different ranges; narrower range scores higher.
        s_narrow = compute(_grid([[200.0]]), _grid([[160.0]])).to_numpy()[0, 0]
        s_wide = compute(_grid([[250.0]]), _grid([[110.0]])).to_numpy()[0, 0]
        assert s_narrow > s_wide
