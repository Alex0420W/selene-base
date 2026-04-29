"""Tests for :mod:`selene_base.criteria.hazard`."""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.hazard import compute


def _da(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


class TestCompute:
    def test_zero_density_scores_one(self) -> None:
        out = compute(_da([[0.0, 0.0]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_saturation_scores_zero(self) -> None:
        out = compute(_da([[50.0]]), saturation_count=50.0).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_above_saturation_clipped_to_zero(self) -> None:
        out = compute(_da([[100.0, 1000.0]]), saturation_count=50.0).to_numpy()
        np.testing.assert_array_equal(out, [[0.0, 0.0]])

    def test_linear_ramp_in_between(self) -> None:
        out = compute(_da([[0.0, 12.5, 25.0, 50.0]]), saturation_count=50.0).to_numpy()
        np.testing.assert_allclose(out, [[1.0, 0.75, 0.5, 0.0]], atol=1e-9)

    def test_nan_propagates(self) -> None:
        out = compute(_da([[np.nan, 0.0]])).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    def test_saturation_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="saturation_count"):
            compute(_da([[0.0]]), saturation_count=0.0)
