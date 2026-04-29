"""Tests for :mod:`selene_base.criteria.thermal` (PRP-based, week 6+)."""

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
        out = compute(_grid([[230.0]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_far_from_target_drops_score(self) -> None:
        out = compute(_grid([[300.0]])).to_numpy()[0, 0]
        # exp(-((300-230)^2) / (2*50^2)) = exp(-0.98) ~ 0.375
        assert 0.3 < out < 0.45

    def test_score_decays_with_offset(self) -> None:
        s_target = compute(_grid([[230.0]])).to_numpy()[0, 0]
        s_offset = compute(_grid([[180.0]])).to_numpy()[0, 0]
        assert s_target == pytest.approx(1.0)
        assert s_offset < s_target
        # exp(-((180-230)^2)/(2*50^2)) = exp(-0.5) ~ 0.6065
        assert 0.55 < s_offset < 0.65

    def test_nan_propagates(self) -> None:
        out = compute(_grid([[np.nan, 230.0]])).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    @pytest.mark.parametrize("kwarg", ["target_temp_k", "sigma_k"])
    def test_non_positive_params_rejected(self, kwarg: str) -> None:
        with pytest.raises(ValueError, match=kwarg):
            compute(_grid([[230.0]]), **{kwarg: 0.0})

    def test_custom_target(self) -> None:
        # If we re-target the criterion to 200 K, that's where it scores 1.
        out = compute(_grid([[200.0]]), target_temp_k=200.0).to_numpy()
        np.testing.assert_allclose(out, 1.0)
