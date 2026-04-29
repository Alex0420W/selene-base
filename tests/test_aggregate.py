"""Tests for :func:`selene_base.scoring.aggregate.weighted_sum`.

Covers the renormalisation behaviour added in week 2: when ``weights``
references criteria that aren't yet available in ``scores``, a warning
fires and the remaining weights are renormalised to sum to 1.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import xarray as xr

from selene_base.scoring.aggregate import weighted_sum


def _grid(value: float, shape: tuple[int, int] = (4, 4)) -> xr.DataArray:
    return xr.DataArray(np.full(shape, value), dims=("y", "x"))


class TestWeightedSum:
    def test_two_criteria_weighted_correctly(self) -> None:
        scores = {"a": _grid(1.0), "b": _grid(0.0)}
        out = weighted_sum(scores, {"a": 0.3, "b": 0.7})
        np.testing.assert_allclose(out.to_numpy(), 0.3)

    def test_renormalises_when_weights_dont_sum_to_one(self) -> None:
        scores = {"a": _grid(1.0), "b": _grid(0.0)}
        out = weighted_sum(scores, {"a": 1.0, "b": 4.0})
        np.testing.assert_allclose(out.to_numpy(), 0.2)  # 1/(1+4)

    def test_missing_criterion_in_weights_raises(self) -> None:
        with pytest.raises(KeyError, match="missing weights"):
            weighted_sum({"a": _grid(1.0)}, {"b": 1.0})

    def test_empty_scores_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            weighted_sum({}, {"a": 1.0})

    def test_negative_weight_rejected(self) -> None:
        with pytest.raises(ValueError, match="negative"):
            weighted_sum({"a": _grid(1.0)}, {"a": -0.5})

    def test_zero_total_weight_rejected(self) -> None:
        with pytest.raises(ValueError, match="zero"):
            weighted_sum({"a": _grid(1.0)}, {"a": 0.0})

    # --- New behaviour for week 2: tolerate weights for absent criteria ---

    def test_extra_weight_keys_warn_but_do_not_raise(self) -> None:
        scores = {"slope": _grid(0.5)}
        weights = {"slope": 0.15, "ice": 0.25, "illumination": 0.30}
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = weighted_sum(scores, weights)
        # exactly one warning that names ice and illumination
        assert any("ice" in str(w.message) and "illumination" in str(w.message) for w in caught), (
            f"expected warning naming ice + illumination, got: {[str(w.message) for w in caught]}"
        )
        # only slope is present, so its weight renormalises to 1
        np.testing.assert_allclose(out.to_numpy(), 0.5)

    def test_partial_subset_renormalises(self) -> None:
        # default-ish weights, but only two criteria present
        weights = {
            "illumination": 0.30,
            "ice": 0.25,
            "slope": 0.15,
            "thermal": 0.10,
            "hazard": 0.10,
            "seismic": 0.10,
        }
        scores = {"slope": _grid(1.0), "illumination": _grid(0.5)}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            out = weighted_sum(scores, weights)
        # renormalised: slope_w = 0.15/(0.15+0.30) = 1/3, illum_w = 2/3
        expected = (1 / 3) * 1.0 + (2 / 3) * 0.5
        np.testing.assert_allclose(out.to_numpy(), expected)
