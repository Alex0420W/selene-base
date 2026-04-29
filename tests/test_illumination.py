"""Tests for :mod:`selene_base.criteria.illumination`."""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.illumination import compute


def _grid(values: list[float]) -> xr.DataArray:
    arr = np.asarray(values, dtype=np.float64).reshape(1, -1)
    return xr.DataArray(arr, dims=("y", "x"))


class TestCompute:
    def test_zero_illumination_scores_zero(self) -> None:
        out = compute(_grid([0.0])).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_target_pct_scores_one(self) -> None:
        out = compute(_grid([0.70]), target_pct=0.70).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_above_target_clipped_to_one(self) -> None:
        out = compute(_grid([0.90, 1.0]), target_pct=0.70).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_linear_below_target(self) -> None:
        out = compute(_grid([0.0, 0.35, 0.70]), target_pct=0.70).to_numpy()
        np.testing.assert_allclose(out, [[0.0, 0.5, 1.0]], atol=1e-9)

    def test_custom_target_changes_score(self) -> None:
        s_default = compute(_grid([0.5]), target_pct=0.70).to_numpy()[0, 0]
        s_strict = compute(_grid([0.5]), target_pct=1.0).to_numpy()[0, 0]
        # Stricter target → lower score for the same input
        assert s_strict < s_default

    def test_nan_propagates(self) -> None:
        out = compute(_grid([np.nan, 0.5])).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(0.5 / 0.70)

    def test_target_pct_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="target_pct"):
            compute(_grid([0.5]), target_pct=0.0)
