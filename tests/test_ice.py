"""Tests for :mod:`selene_base.criteria.ice`."""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.ice import compute, derive_psr_mask


def _da(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


class TestDerivePsrMask:
    def test_zero_illumination_is_psr(self) -> None:
        out = derive_psr_mask(_da([[0.0, 0.5, 1.0]])).to_numpy()
        assert out.tolist() == [[True, False, False]]

    def test_threshold_boundary(self) -> None:
        # threshold default 0.001 — anything below is PSR, anything ≥ is not
        out = derive_psr_mask(_da([[0.0, 0.0009, 0.001, 0.002]])).to_numpy()
        assert out.tolist() == [[True, True, False, False]]

    def test_nan_is_not_psr(self) -> None:
        out = derive_psr_mask(_da([[np.nan, 0.0]])).to_numpy()
        assert out.tolist() == [[False, True]]

    @pytest.mark.parametrize("bad", [0.0, 1.0, -0.1, 1.5])
    def test_threshold_validated(self, bad: float) -> None:
        with pytest.raises(ValueError, match="threshold"):
            derive_psr_mask(_da([[0.0]]), threshold=bad)


class TestCompute:
    def test_base_score_is_inverse_min_max(self) -> None:
        # flux ranges 0..10 → score ranges 1..0
        flux = _da([[0.0, 5.0, 10.0]])
        out = compute(flux).to_numpy()
        np.testing.assert_allclose(out, [[1.0, 0.5, 0.0]], atol=1e-9)

    def test_nan_propagates(self) -> None:
        flux = _da([[np.nan, 0.0, 10.0]])
        out = compute(flux).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(1.0)
        assert out[0, 2] == pytest.approx(0.0)

    def test_psr_bonus_increases_score_near_psr(self) -> None:
        # Varying flux so base score isn't pinned at one extreme:
        # flux 10..0 → base 0..1 across the row.
        flux = _da([[10.0, 5.0, 5.0, 0.0]])
        psr_mask = xr.DataArray(np.array([[True, False, False, False]]), dims=("y", "x"))
        # All four cols are within 5 km of col 0 (since pixel_size = 240 m).
        out = compute(
            flux,
            psr_mask=psr_mask,
            near_psr_radius_km=5.0,
            near_psr_bonus=0.3,
            pixel_size_m=240.0,
        ).to_numpy()
        # base = 1 - min_max(flux) = [0, 0.5, 0.5, 1]; +0.3 clipped to [0.3, 0.8, 0.8, 1.0]
        np.testing.assert_allclose(out, [[0.3, 0.8, 0.8, 1.0]], atol=1e-9)

    def test_psr_bonus_decays_outside_radius(self) -> None:
        # 1 pixel = 1000 m; radius 1.5 km means rows 0,1 in radius, rows 2+ out.
        # Use varying flux so base is well-defined per row:
        # flux 10, 0, 5, 5 → base 0, 1, 0.5, 0.5
        flux = _da([[10.0], [0.0], [5.0], [5.0]])
        psr_mask = xr.DataArray(np.array([[True], [False], [False], [False]]), dims=("y", "x"))
        out = compute(
            flux,
            psr_mask=psr_mask,
            near_psr_radius_km=1.5,
            near_psr_bonus=0.3,
            pixel_size_m=1000.0,
        ).to_numpy()
        assert out[0, 0] == pytest.approx(0.3)  # base 0 + bonus 0.3
        assert out[1, 0] == pytest.approx(1.0)  # base 1 + bonus 0.3 clipped to 1
        assert out[2, 0] == pytest.approx(0.5)  # outside radius, base 0.5 only
        assert out[3, 0] == pytest.approx(0.5)

    def test_bonus_clipped_to_one(self) -> None:
        flux = _da([[0.0, 10.0]])  # base 1.0, 0.0
        psr_mask = xr.DataArray(np.array([[True, True]]), dims=("y", "x"))
        out = compute(flux, psr_mask=psr_mask, near_psr_bonus=0.5, pixel_size_m=10.0).to_numpy()
        # left base=1.0 + bonus 0.5 → clipped to 1.0
        assert out[0, 0] == pytest.approx(1.0)
        assert out[0, 1] == pytest.approx(0.5)

    def test_psr_shape_mismatch_raises(self) -> None:
        flux = _da([[1.0, 2.0]])
        bad_mask = xr.DataArray(np.array([[True]]), dims=("y", "x"))
        with pytest.raises(ValueError, match="shape"):
            compute(flux, psr_mask=bad_mask)

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
            compute(_da([[1.0]]), **{kwarg: bad})
