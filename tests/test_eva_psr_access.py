"""Tests for :mod:`selene_base.criteria.eva_psr_access` (v2.0+)."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.eva_psr_access import _build_disc_kernel, compute


def _grid(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


class TestComputeBasics:
    def test_score_is_in_range_zero_to_one(self) -> None:
        # 50/50 mix of warm and cold cells across a 9x9 grid; every score
        # must land in [0, 1] regardless of the mix.
        rng = np.random.default_rng(0)
        temps = rng.uniform(50.0, 200.0, size=(9, 9))
        out = compute(_grid(temps.tolist()), eva_radius_km=2.0, pixel_size_m=240.0)
        arr = out.to_numpy()
        assert np.all(arr >= 0.0)
        assert np.all(arr <= 1.0)

    def test_score_zero_when_no_cold_neighbours(self) -> None:
        # Uniformly warm grid (T = 200 K well above 110 K cap): every
        # disc contains zero cold cells, so every score is 0.
        warm = np.full((11, 11), 200.0)
        out = compute(_grid(warm.tolist())).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_score_one_when_entire_disc_is_cold(self) -> None:
        # Uniformly cold grid (T = 50 K, well below 110 K): every disc
        # is 100% cold, score is 1.0 everywhere.
        cold = np.full((11, 11), 50.0)
        out = compute(_grid(cold.tolist())).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_partial_disc_scores_between_zero_and_one(self) -> None:
        # Half-cold, half-warm vertical split: cells far from the boundary
        # should see ~50% cold neighbours; the centre column straddles the
        # boundary but the disc is symmetric so the score is around 0.5.
        # Build a 21x21 grid so the 9-pixel-radius (2 km / 240 m ≈ 8.33)
        # disc fits cleanly in the interior.
        arr = np.full((21, 21), 50.0)
        arr[:, 11:] = 200.0  # right half warm
        out = compute(_grid(arr.tolist()), eva_radius_km=2.0, pixel_size_m=240.0).to_numpy()
        # Centre column straddles the boundary -> ~50% cold.
        assert 0.4 < out[10, 10] < 0.6
        # Far-left column entirely inside the cold half (disc radius 9, col 0
        # has cold neighbours within the disc only).
        assert out[10, 0] >= 0.95
        # Far-right column entirely inside the warm half.
        assert out[10, 20] <= 0.05

    def test_cell_in_psr_with_warm_surroundings_does_not_force_one(self) -> None:
        # A single cold cell in a warm sea: the cold cell scores
        # 1/disc-cells, not 1.0 — the criterion measures *area*, not
        # mere presence.
        arr = np.full((21, 21), 200.0)
        arr[10, 10] = 50.0
        out = compute(_grid(arr.tolist()), eva_radius_km=2.0, pixel_size_m=240.0).to_numpy()
        assert 0.0 < out[10, 10] < 0.05


class TestEdges:
    def test_truncated_disc_at_boundary_uses_actual_neighbours(self) -> None:
        # All-cold grid; corner cells have ~1/4 of the disc cut off but
        # all surviving neighbours are cold, so the ratio is still 1.0.
        cold = np.full((11, 11), 50.0)
        out = compute(_grid(cold.tolist())).to_numpy()
        # Corner cells.
        assert out[0, 0] == pytest.approx(1.0)
        assert out[10, 10] == pytest.approx(1.0)

    def test_nan_input_propagates_to_nan_output(self) -> None:
        arr = np.full((11, 11), 50.0)
        arr[5, 5] = np.nan
        out = compute(_grid(arr.tolist())).to_numpy()
        assert np.isnan(out[5, 5])
        # Surrounding cells, while their disc has one fewer valid
        # neighbour, are still finite.
        assert np.isfinite(out[5, 4])
        assert np.isfinite(out[4, 5])

    def test_all_nan_grid_is_all_nan(self) -> None:
        arr = np.full((7, 7), np.nan)
        out = compute(_grid(arr.tolist())).to_numpy()
        assert np.all(np.isnan(out))


class TestParameters:
    def test_eva_radius_parameter_respected(self) -> None:
        # Single cold cell at the centre; larger discs include it across
        # progressively more candidate cells, so the count of cells with
        # any non-zero score grows monotonically with radius.
        arr = np.full((41, 41), 200.0)
        arr[20, 20] = 50.0
        out_1km = compute(_grid(arr.tolist()), eva_radius_km=1.0, pixel_size_m=240.0).to_numpy()
        out_2km = compute(_grid(arr.tolist()), eva_radius_km=2.0, pixel_size_m=240.0).to_numpy()
        out_5km = compute(_grid(arr.tolist()), eva_radius_km=5.0, pixel_size_m=240.0).to_numpy()
        n_nonzero_1 = int((out_1km > 0).sum())
        n_nonzero_2 = int((out_2km > 0).sum())
        n_nonzero_5 = int((out_5km > 0).sum())
        assert n_nonzero_1 < n_nonzero_2 < n_nonzero_5

    def test_cold_threshold_parameter_respected(self) -> None:
        # A cell at 100 K is cold-class at 110 K cap but not at 90 K.
        arr = np.full((11, 11), 200.0)
        arr[5, 5] = 100.0
        out_110 = compute(
            _grid(arr.tolist()), cold_threshold_k=110.0, pixel_size_m=240.0
        ).to_numpy()
        out_90 = compute(_grid(arr.tolist()), cold_threshold_k=90.0, pixel_size_m=240.0).to_numpy()
        assert out_110[5, 5] > 0
        assert out_90[5, 5] == 0

    def test_score_consistency_across_resolutions(self) -> None:
        # A 50/50 split at 240 m and the same split rendered at 80 m
        # (3x finer pixel size, same physical extent) should produce
        # similar scores at corresponding geographic points. We compare
        # the centre-column score on each grid; both should be near 0.5.
        # 240 m grid: 21x21
        arr_240 = np.full((21, 21), 50.0)
        arr_240[:, 11:] = 200.0
        out_240 = compute(_grid(arr_240.tolist()), pixel_size_m=240.0).to_numpy()
        # 80 m grid: 63x63 (3x linear), same boundary at column 32.
        arr_80 = np.full((63, 63), 50.0)
        arr_80[:, 32:] = 200.0
        out_80 = compute(_grid(arr_80.tolist()), eva_radius_km=2.0, pixel_size_m=80.0).to_numpy()
        # Centre score must be near 0.5 on both grids; tolerance accounts
        # for kernel-discretisation differences.
        assert abs(out_240[10, 10] - 0.5) < 0.1
        assert abs(out_80[31, 32] - 0.5) < 0.1

    def test_invalid_eva_radius_raises(self) -> None:
        with pytest.raises(ValueError, match="eva_radius_km"):
            compute(_grid([[100.0]]), eva_radius_km=0.0)

    def test_invalid_cold_threshold_raises(self) -> None:
        with pytest.raises(ValueError, match="cold_threshold_k"):
            compute(_grid([[100.0]]), cold_threshold_k=-5.0)

    def test_invalid_pixel_size_raises(self) -> None:
        with pytest.raises(ValueError, match="pixel_size_m"):
            compute(_grid([[100.0]]), pixel_size_m=0.0)


class TestKernel:
    def test_disc_kernel_is_circular(self) -> None:
        # A radius-3 kernel covers all pixels with dy^2+dx^2 <= 9; the
        # corners (3, 3) at distance sqrt(18) > 3 must be excluded.
        kernel = _build_disc_kernel(3.0)
        assert kernel.shape == (7, 7)
        # Centre included.
        assert kernel[3, 3] == 1.0
        # Corner excluded (sqrt(18) > 3).
        assert kernel[0, 0] == 0.0
        # Cardinal-direction edge included.
        assert kernel[0, 3] == 1.0

    def test_zero_radius_kernel_is_single_cell(self) -> None:
        kernel = _build_disc_kernel(0.0)
        assert kernel.shape == (1, 1)
        assert kernel[0, 0] == 1.0
