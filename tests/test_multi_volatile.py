"""Tests for :mod:`selene_base.criteria.multi_volatile` (v2.1+)."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.multi_volatile import (
    CO2_NH3_THRESHOLD_K,
    DEFAULT_SUB_WEIGHTS,
    H2O_THRESHOLD_K,
    ULTRACOLD_THRESHOLD_K,
    _build_disc_kernel,
    compute,
    compute_components,
)


def _grid(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


# ----- thresholds -----


class TestThresholds:
    def test_thresholds_strictly_decreasing(self) -> None:
        # Nesting: ultracold < co2_nh3 < h2o means ultracold ⊂ co2_nh3 ⊂ h2o.
        assert ULTRACOLD_THRESHOLD_K < CO2_NH3_THRESHOLD_K < H2O_THRESHOLD_K

    def test_default_sub_weights_equal_thirds(self) -> None:
        assert DEFAULT_SUB_WEIGHTS == (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)

    def test_h2o_threshold_is_110(self) -> None:
        assert H2O_THRESHOLD_K == 110.0

    def test_co2_nh3_threshold_is_66(self) -> None:
        assert CO2_NH3_THRESHOLD_K == 66.0

    def test_ultracold_threshold_is_60(self) -> None:
        assert ULTRACOLD_THRESHOLD_K == 60.0


class TestThermalClassMembership:
    def test_109k_is_h2o_only(self) -> None:
        # A cell at exactly 109 K and a tight EVA radius (so only the cell
        # itself is inside the disc): h2o yes, co2 no, ultracold no.
        out = compute_components(_grid([[109.0]]), eva_radius_km=0.1, pixel_size_m=240.0)
        assert out["h2o_score"].to_numpy()[0, 0] == pytest.approx(1.0)
        assert out["co2_nh3_score"].to_numpy()[0, 0] == pytest.approx(0.0)
        assert out["ultracold_score"].to_numpy()[0, 0] == pytest.approx(0.0)
        # combined = (1 + 0 + 0) / 3
        assert out["combined_score"].to_numpy()[0, 0] == pytest.approx(1.0 / 3.0)

    def test_65k_is_h2o_and_co2(self) -> None:
        out = compute_components(_grid([[65.0]]), eva_radius_km=0.1, pixel_size_m=240.0)
        assert out["h2o_score"].to_numpy()[0, 0] == pytest.approx(1.0)
        assert out["co2_nh3_score"].to_numpy()[0, 0] == pytest.approx(1.0)
        assert out["ultracold_score"].to_numpy()[0, 0] == pytest.approx(0.0)
        assert out["combined_score"].to_numpy()[0, 0] == pytest.approx(2.0 / 3.0)

    def test_59k_is_all_three(self) -> None:
        out = compute_components(_grid([[59.0]]), eva_radius_km=0.1, pixel_size_m=240.0)
        assert out["h2o_score"].to_numpy()[0, 0] == pytest.approx(1.0)
        assert out["co2_nh3_score"].to_numpy()[0, 0] == pytest.approx(1.0)
        assert out["ultracold_score"].to_numpy()[0, 0] == pytest.approx(1.0)
        assert out["combined_score"].to_numpy()[0, 0] == pytest.approx(1.0)

    def test_threshold_boundaries_are_strict(self) -> None:
        # The thresholds use ``<`` (strict). At exactly 110 K, H₂O class
        # does not pass; at exactly 66 K, CO₂ does not pass; at 60 K,
        # ultracold does not pass.
        out_110 = compute_components(_grid([[110.0]]), eva_radius_km=0.1, pixel_size_m=240.0)
        assert out_110["h2o_score"].to_numpy()[0, 0] == pytest.approx(0.0)
        out_66 = compute_components(_grid([[66.0]]), eva_radius_km=0.1, pixel_size_m=240.0)
        assert out_66["co2_nh3_score"].to_numpy()[0, 0] == pytest.approx(0.0)
        out_60 = compute_components(_grid([[60.0]]), eva_radius_km=0.1, pixel_size_m=240.0)
        assert out_60["ultracold_score"].to_numpy()[0, 0] == pytest.approx(0.0)


# ----- score range and basics -----


class TestScoreRange:
    def test_score_in_range_zero_to_one(self) -> None:
        rng = np.random.default_rng(0)
        temps = rng.uniform(40.0, 200.0, size=(15, 15))
        out = compute(_grid(temps.tolist()), eva_radius_km=2.0, pixel_size_m=240.0)
        arr = out.to_numpy()
        assert np.all(arr >= 0.0)
        assert np.all(arr <= 1.0)

    def test_score_zero_no_cold_neighbors(self) -> None:
        # Uniformly warm grid (T = 200 K, well above all thresholds):
        # every disc fully fails every class, so every score is 0.
        warm = np.full((11, 11), 200.0)
        out = compute(_grid(warm.tolist())).to_numpy()
        np.testing.assert_allclose(out, 0.0)

    def test_score_one_all_classes_covered(self) -> None:
        # Uniformly ultracold (T = 50 K, below every threshold): every
        # class passes everywhere, every score is 1.
        cold = np.full((11, 11), 50.0)
        out = compute(_grid(cold.tolist())).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_h2o_only_partial_score(self) -> None:
        # Uniformly H₂O-but-not-CO₂ (T = 80 K). Combined = 1/3.
        warm_h2o = np.full((11, 11), 80.0)
        out = compute(_grid(warm_h2o.tolist())).to_numpy()
        np.testing.assert_allclose(out, 1.0 / 3.0)

    def test_h2o_and_co2_partial_score(self) -> None:
        # T = 65 K (passes H₂O and CO₂, fails ultracold). Combined = 2/3.
        cool = np.full((11, 11), 65.0)
        out = compute(_grid(cool.tolist())).to_numpy()
        np.testing.assert_allclose(out, 2.0 / 3.0)


class TestEdges:
    def test_truncated_disc_at_boundary_uses_actual_neighbours(self) -> None:
        # All-ultracold grid; corner cells have ~1/4 of the disc cut off
        # but every surviving neighbour passes every class, so the ratio
        # remains 1.0 on every sub-score.
        cold = np.full((11, 11), 50.0)
        out = compute_components(_grid(cold.tolist()))
        for key in ("h2o_score", "co2_nh3_score", "ultracold_score", "combined_score"):
            arr = out[key].to_numpy()
            assert arr[0, 0] == pytest.approx(1.0)
            assert arr[10, 10] == pytest.approx(1.0)

    def test_nan_input_propagates_to_nan_output(self) -> None:
        arr = np.full((11, 11), 50.0)
        arr[5, 5] = np.nan
        out = compute(_grid(arr.tolist())).to_numpy()
        assert np.isnan(out[5, 5])
        # Surrounding cells (with one fewer valid neighbour) still
        # finite.
        assert np.isfinite(out[5, 4])

    def test_all_nan_grid_is_all_nan(self) -> None:
        arr = np.full((7, 7), np.nan)
        out = compute(_grid(arr.tolist())).to_numpy()
        assert np.all(np.isnan(out))


# ----- parameters -----


class TestParameters:
    def test_components_returns_four_keys(self) -> None:
        out = compute_components(_grid([[50.0]]))
        assert set(out.keys()) == {
            "h2o_score",
            "co2_nh3_score",
            "ultracold_score",
            "combined_score",
        }

    def test_custom_sub_weights_applied(self) -> None:
        # NASA-priority-style weights (h2o heavier than ultracold). At a
        # uniformly H₂O-only cell (80 K), only the H₂O sub-score
        # contributes, so combined = w_h2o (after renormalisation).
        warm_h2o = np.full((9, 9), 80.0)
        out = compute(_grid(warm_h2o.tolist()), sub_weights=(0.5, 0.3, 0.2)).to_numpy()
        # h2o_score = 1, co2 = 0, ultra = 0; combined = 0.5/(0.5+0.3+0.2) = 0.5
        np.testing.assert_allclose(out, 0.5)

    def test_sub_weights_renormalised_when_unnormalised(self) -> None:
        # (1, 0, 0) means "score on H₂O only". Internally renormalises
        # to (1, 0, 0); at a H₂O-only cell combined = 1.
        warm_h2o = np.full((9, 9), 80.0)
        out = compute(_grid(warm_h2o.tolist()), sub_weights=(1.0, 0.0, 0.0)).to_numpy()
        np.testing.assert_allclose(out, 1.0)

    def test_eva_radius_parameter_respected(self) -> None:
        # Single ultracold cell at the centre of a warm grid; larger
        # discs incorporate it across more candidate cells, so the
        # number of cells with non-zero score grows monotonically.
        arr = np.full((41, 41), 200.0)
        arr[20, 20] = 50.0
        n_at = lambda r: int(  # noqa: E731
            (compute(_grid(arr.tolist()), eva_radius_km=r, pixel_size_m=240.0).to_numpy() > 0).sum()
        )
        assert n_at(1.0) < n_at(2.0) < n_at(5.0)

    def test_pixel_size_parameter_respected(self) -> None:
        # Same physical EVA radius (2 km) at two different pixel sizes
        # gives roughly the same fraction of cold-class disc cells in
        # the interior — the score is scale-invariant.
        arr_240 = np.full((21, 21), 50.0)
        arr_240[:, 11:] = 200.0
        out_240 = compute(_grid(arr_240.tolist()), pixel_size_m=240.0).to_numpy()
        arr_80 = np.full((63, 63), 50.0)
        arr_80[:, 32:] = 200.0
        out_80 = compute(_grid(arr_80.tolist()), pixel_size_m=80.0).to_numpy()
        # Centre score ~0.5 on both grids; tolerance for kernel
        # discretisation differences.
        assert abs(out_240[10, 10] - 0.5) < 0.1
        assert abs(out_80[31, 32] - 0.5) < 0.1


class TestParameterValidation:
    def test_invalid_eva_radius_raises(self) -> None:
        with pytest.raises(ValueError, match="eva_radius_km"):
            compute(_grid([[100.0]]), eva_radius_km=0.0)

    def test_invalid_pixel_size_raises(self) -> None:
        with pytest.raises(ValueError, match="pixel_size_m"):
            compute(_grid([[100.0]]), pixel_size_m=0.0)

    def test_non_monotonic_thresholds_raise(self) -> None:
        # Misordered thresholds (ultracold > co2) violate the nesting
        # contract.
        with pytest.raises(ValueError, match="ultracold ≤ co2_nh3 ≤ h2o"):
            compute(
                _grid([[50.0]]),
                ultracold_threshold_k=80.0,
                co2_nh3_threshold_k=66.0,
                h2o_threshold_k=110.0,
            )

    def test_negative_sub_weight_raises(self) -> None:
        with pytest.raises(ValueError, match="sub_weights must be non-negative"):
            compute(_grid([[50.0]]), sub_weights=(-0.1, 0.5, 0.6))

    def test_zero_sum_sub_weights_raises(self) -> None:
        with pytest.raises(ValueError, match="positive sum"):
            compute(_grid([[50.0]]), sub_weights=(0.0, 0.0, 0.0))


class TestKernel:
    def test_disc_kernel_is_circular(self) -> None:
        kernel = _build_disc_kernel(3.0)
        assert kernel.shape == (7, 7)
        assert kernel[3, 3] == 1.0  # centre
        assert kernel[0, 0] == 0.0  # corner outside disc
        assert kernel[0, 3] == 1.0  # cardinal edge

    def test_zero_radius_kernel_is_single_cell(self) -> None:
        kernel = _build_disc_kernel(0.0)
        assert kernel.shape == (1, 1)


# ----- nesting invariant -----


class TestNestingInvariant:
    def test_per_cell_ultracold_le_co2_le_h2o(self) -> None:
        # Random temperatures across a 31x31 grid; for every cell in
        # the interior (no edge effects), ultracold_score must be
        # ≤ co2_nh3_score ≤ h2o_score because the masks are nested.
        rng = np.random.default_rng(42)
        temps = rng.uniform(30.0, 250.0, size=(31, 31))
        out = compute_components(_grid(temps.tolist()))
        h2o = out["h2o_score"].to_numpy()
        co2 = out["co2_nh3_score"].to_numpy()
        ultra = out["ultracold_score"].to_numpy()
        # Use a small tolerance for float32 round-off.
        assert np.all(ultra <= co2 + 1e-6)
        assert np.all(co2 <= h2o + 1e-6)


class TestThreeBandSyntheticScene:
    def test_three_distinct_classes_score_correctly(self) -> None:
        # 21x21 grid split into three temperature bands corresponding to
        # the three thermal classes. Verifies the criterion correctly
        # distinguishes them.
        arr = np.full((21, 21), 200.0)
        arr[:7, :] = 50.0  # ultracold band → all three classes pass
        arr[7:14, :] = 80.0  # H₂O-only band
        arr[14:, :] = 200.0  # warm band → no class passes
        # Tight EVA radius so each cell sees only itself + immediate neighbours.
        out = compute_components(_grid(arr.tolist()), eva_radius_km=0.1, pixel_size_m=240.0)
        # Centres of each band: ultracold band scores 1.0, H₂O-only 1/3,
        # warm 0.0.
        assert float(out["combined_score"][3, 10]) == pytest.approx(1.0)
        assert float(out["combined_score"][10, 10]) == pytest.approx(1.0 / 3.0)
        assert float(out["combined_score"][17, 10]) == pytest.approx(0.0)
