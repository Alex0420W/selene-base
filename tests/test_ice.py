"""Tests for :mod:`selene_base.criteria.ice`.

Two compute paths:

- :func:`compute` — PRP-based (default since week 6).
- :func:`compute_from_lend` — kept as a drop-in for when LEND CSETN
  flux maps land.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.ice import compute, compute_from_lend, derive_psr_mask


def _da(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


class TestDerivePsrMask:
    def test_zero_illumination_is_psr(self) -> None:
        out = derive_psr_mask(_da([[0.0, 0.5, 1.0]])).to_numpy()
        assert out.tolist() == [[True, False, False]]

    def test_threshold_boundary(self) -> None:
        out = derive_psr_mask(_da([[0.0, 0.0009, 0.001, 0.002]])).to_numpy()
        assert out.tolist() == [[True, True, False, False]]

    def test_nan_is_not_psr(self) -> None:
        out = derive_psr_mask(_da([[np.nan, 0.0]])).to_numpy()
        assert out.tolist() == [[False, True]]

    @pytest.mark.parametrize("bad", [0.0, 1.0, -0.1, 1.5])
    def test_threshold_validated(self, bad: float) -> None:
        with pytest.raises(ValueError, match="threshold"):
            derive_psr_mask(_da([[0.0]]), threshold=bad)


class TestComputePRP:
    """The default PRP-based compute(): ice_depth_m -> [0, 1] score."""

    def test_surface_ice_max_score(self) -> None:
        # ice_depth = 0 -> base 1.0 + surface bonus 0.5 -> clipped to 1.0
        out = compute(_da([[0.0]])).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_no_ice_scores_zero(self) -> None:
        # NaN ice_depth (sentinel -999 in source) -> 0
        out = compute(_da([[np.nan]])).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_linear_ramp_with_depth(self) -> None:
        # No bonuses; pure base = 1 - depth/2.87
        depths = _da([[0.0, 0.7175, 1.435, 2.1525, 2.87]])
        out = compute(depths, surface_ice_bonus=0.0).to_numpy()
        np.testing.assert_allclose(out, [[1.0, 0.75, 0.5, 0.25, 0.0]], atol=1e-9)

    def test_psr_proximity_bonus(self) -> None:
        depths = _da([[0.5, 0.5, 0.5, 0.5]])  # base 1 - 0.5/2.87 ~ 0.826
        psr = xr.DataArray(np.array([[True, False, False, False]]), dims=("y", "x"))
        # Pixel size 240 m, 5 km radius -> all 4 cells in range.
        out = compute(
            depths,
            psr_mask=psr,
            surface_ice_bonus=0.0,
            near_psr_bonus=0.1,
            near_psr_radius_km=5.0,
            pixel_size_m=240.0,
        ).to_numpy()
        # base ~0.826; +0.1 PSR bonus = 0.926 (PSR cell itself also gets the bonus)
        np.testing.assert_allclose(out, 0.826 + 0.1, atol=5e-3)

    def test_score_bounded_to_unit(self) -> None:
        # Generous bonuses + surface ice should clip at 1.0, never exceed.
        depths = _da([[0.0]])
        psr = xr.DataArray(np.array([[True]]), dims=("y", "x"))
        out = compute(
            depths,
            psr_mask=psr,
            surface_ice_bonus=0.8,
            near_psr_bonus=0.5,
            pixel_size_m=10.0,
        ).to_numpy()
        assert out[0, 0] == pytest.approx(1.0)

    def test_psr_shape_mismatch_rejected(self) -> None:
        bad = xr.DataArray(np.array([[True, True]]), dims=("y", "x"))
        with pytest.raises(ValueError, match="shape"):
            compute(_da([[0.0]]), psr_mask=bad)

    @pytest.mark.parametrize(
        "kwarg,bad",
        [
            ("surface_ice_bonus", 1.5),
            ("surface_ice_bonus", -0.1),
            ("near_psr_bonus", 1.5),
            ("near_psr_radius_km", 0.0),
            ("pixel_size_m", 0.0),
            ("max_depth_m", 0.0),
        ],
    )
    def test_param_validation(self, kwarg: str, bad: float) -> None:
        with pytest.raises(ValueError, match=kwarg):
            compute(_da([[0.0]]), **{kwarg: bad})


class TestComputeFromLEND:
    """Backward-compat path; activates when LEND data finally lands."""

    def test_base_score_is_inverse_min_max(self) -> None:
        flux = _da([[0.0, 5.0, 10.0]])
        out = compute_from_lend(flux).to_numpy()
        np.testing.assert_allclose(out, [[1.0, 0.5, 0.0]], atol=1e-9)

    def test_nan_propagates(self) -> None:
        flux = _da([[np.nan, 0.0, 10.0]])
        out = compute_from_lend(flux).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)

    def test_psr_bonus_applied(self) -> None:
        flux = _da([[10.0, 5.0, 5.0, 0.0]])
        psr = xr.DataArray(np.array([[True, False, False, False]]), dims=("y", "x"))
        out = compute_from_lend(
            flux,
            psr_mask=psr,
            near_psr_bonus=0.3,
            pixel_size_m=240.0,
        ).to_numpy()
        # base = [0, 0.5, 0.5, 1]; +0.3 clipped -> [0.3, 0.8, 0.8, 1.0]
        np.testing.assert_allclose(out, [[0.3, 0.8, 0.8, 1.0]], atol=1e-9)

    @pytest.mark.parametrize(
        "kwarg,bad",
        [
            ("near_psr_radius_km", 0.0),
            ("near_psr_bonus", 1.5),
            ("near_psr_bonus", -0.1),
            ("pixel_size_m", 0.0),
        ],
    )
    def test_param_validation(self, kwarg: str, bad: float) -> None:
        with pytest.raises(ValueError, match=kwarg):
            compute_from_lend(_da([[1.0]]), **{kwarg: bad})
