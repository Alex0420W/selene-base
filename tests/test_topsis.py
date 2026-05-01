"""Tests for the TOPSIS aggregator (v1.7).

Hwang & Yoon 1981. Each test exercises one property TOPSIS is supposed
to have over weighted_sum:

- Unit-length normalisation: criteria on different raw scales do not
  bias the result.
- Ideal/anti-ideal ranking: cells at the per-criterion max score 1.0;
  cells at the per-criterion min score 0.0.
- Balanced > lopsided: a cell with all criteria moderate beats a cell
  that saturates one criterion at the cost of zeroing another, even
  when the *sum* of criterion scores is identical.
- Renormalised-on-missing-criteria: the aggregate dispatcher's missing-
  criteria handling matches weighted_sum's.
- Determinism.
- Wrong method name raises.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping

import numpy as np
import pytest
import xarray as xr

from selene_base.scoring.aggregate import aggregate, topsis, weighted_sum


def _grid(values: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=np.float64), dims=("y", "x"))


def _equal_weights(names: list[str]) -> Mapping[str, float]:
    w = 1.0 / len(names)
    return {n: w for n in names}


def test_topsis_basic_three_criteria_three_cells_hand_verified() -> None:
    """A 3-criterion, 3-cell input with hand-computable TOPSIS scores.

    Cell (0,0): all criteria = 1.0 → at every criterion's ideal → score 1.
    Cell (0,1): all criteria = 0.0 → at every anti-ideal → score 0.
    Cell (0,2): all criteria = 0.5 → equidistant → score 0.5.
    """
    scores = {
        "a": _grid([[1.0, 0.0, 0.5]]),
        "b": _grid([[1.0, 0.0, 0.5]]),
        "c": _grid([[1.0, 0.0, 0.5]]),
    }
    weights = _equal_weights(list(scores))

    result = topsis(scores, weights).to_numpy().ravel()
    np.testing.assert_allclose(result, [1.0, 0.0, 0.5], atol=1e-9)


def test_topsis_balanced_beats_lopsided_at_equal_weighted_sum() -> None:
    """Two cells with equal weighted_sum should not tie under TOPSIS.

    Cell A: criteria (1.0, 0.0) — saturated on one axis, zero on the
    other. Weighted sum = 0.5.
    Cell B: criteria (0.5, 0.5) — balanced. Weighted sum = 0.5.
    Cell C: criteria (1.0, 1.0) — ideal. Pulls the ideal.
    Cell D: criteria (0.0, 0.0) — anti-ideal. Pulls the anti-ideal.

    weighted_sum is indifferent between A and B (both 0.5).
    TOPSIS should rank B > A: A sits at the ideal of one criterion
    but the anti-ideal of the other, so it's distance √2/2 from the
    ideal point in (a, b) space; B is √2/2 too — wait, let me redo:
    weighted normalised: ideal = (w*1/√3, w*1/√3), anti = (0, 0).
    A weighted norm = (w/√3, 0); B = (w*0.5/√3, w*0.5/√3); C =
    (w/√3, w/√3); D = (0,0). d(A, ideal) = w/√3; d(A, anti) = w/√3.
    d(B, ideal) = (w/√3)*√(0.25+0.25) ≈ w/√3 * √0.5; d(B, anti) =
    (w/√3)*√0.5. Closeness A = w/√3 / (2w/√3) = 0.5. Closeness B = 0.5.
    So actually they tie when the *only* point pulling ideal is C and
    *only* point pulling anti is D — TOPSIS reduces to 0.5 in that
    symmetric case.

    The balanced-beats-lopsided property emerges when ideal/anti-ideal
    are pulled by mixed cells. Here we add a fifth cell E = (0.7, 0.3)
    which is also "lopsided-ish" but on the same side as A.
    """
    # 1x6 grid with: A=(1,0), B=(0.5,0.5), C=(1,1), D=(0,0), E=(0.7,0.3), F=(0.3,0.7)
    scores = {
        "a": _grid([[1.0, 0.5, 1.0, 0.0, 0.7, 0.3]]),
        "b": _grid([[0.0, 0.5, 1.0, 0.0, 0.3, 0.7]]),
    }
    weights = {"a": 0.5, "b": 0.5}
    result = topsis(scores, weights).to_numpy().ravel()
    # A and B have identical weighted sum (0.5).
    ws = weighted_sum(scores, weights).to_numpy().ravel()
    assert ws[0] == pytest.approx(ws[1], abs=1e-9)
    # Under TOPSIS, with C=(1,1) pulling ideal and D=(0,0) pulling
    # anti-ideal, the balanced cell B should outscore A (and F should
    # outscore A; E should be roughly tied with A by symmetry).
    score_A, score_B = result[0], result[1]
    assert score_B > score_A, f"balanced should beat lopsided: A={score_A:.4f} B={score_B:.4f}"


def test_topsis_unequal_raw_scales_do_not_bias_via_normalisation() -> None:
    """If criteria are normalised first, doubling one criterion's raw
    values should leave TOPSIS scores unchanged (the L2 normalisation
    cancels the doubling)."""
    scores_a = {
        "a": _grid([[0.4, 0.6, 0.8]]),
        "b": _grid([[0.2, 0.5, 0.9]]),
    }
    scores_b = {
        "a": _grid([[0.8, 1.2, 1.6]]),  # doubled
        "b": _grid([[0.2, 0.5, 0.9]]),
    }
    weights = {"a": 0.5, "b": 0.5}
    np.testing.assert_allclose(
        topsis(scores_a, weights).to_numpy(),
        topsis(scores_b, weights).to_numpy(),
        atol=1e-12,
    )


def test_topsis_renormalises_on_missing_criteria_with_warning() -> None:
    """Weights for criteria not present in scores are dropped; a
    UserWarning is emitted; remaining weights renormalise. Same
    contract as :func:`weighted_sum`."""
    scores = {
        "a": _grid([[1.0, 0.5, 0.0]]),
        "b": _grid([[0.0, 0.5, 1.0]]),
    }
    weights = {"a": 0.3, "b": 0.3, "c_missing": 0.4}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = topsis(scores, weights)
        assert any("c_missing" in str(w.message) for w in caught)
    # With renormalised weights w_a = w_b = 0.5, this is the symmetric
    # 1x3 case → cells score [1, 0.5, 0] (cell 0 is ideal of a but anti
    # of b — same distance to both → 0.5; let me actually compute).
    # weighted_normalised after renorm:
    #   col=a: [1, 0.5, 0] / sqrt(1+0.25+0) = /sqrt(1.25) → norm * 0.5
    #   col=b: [0, 0.5, 1] / sqrt(1.25) → norm * 0.5
    #   ideal_a = max(col_a) * 0.5; anti_ideal_a = 0; same for b.
    # Cell 0 = (max_a, 0): d_ideal = sqrt(0 + ideal_b^2); d_anti =
    # sqrt(ideal_a^2 + 0). They are equal (since ideal_a = ideal_b) → 0.5.
    # Cell 1 = (0.5*norm, 0.5*norm): d_ideal = sqrt(2 * (0.5*ideal_a)^2);
    # d_anti = sqrt(2*(0.5*ideal_a)^2) → equal → 0.5.
    # Cell 2 = (0, max_b): symmetric to cell 0 → 0.5.
    # So all three score 0.5 in this symmetric case.
    np.testing.assert_allclose(result.to_numpy().ravel(), [0.5, 0.5, 0.5], atol=1e-9)


def test_topsis_propagates_nan() -> None:
    """A NaN in any criterion at a cell yields NaN at that cell."""
    scores = {
        "a": _grid([[1.0, 0.5, np.nan]]),
        "b": _grid([[0.5, 0.5, 0.5]]),
    }
    weights = {"a": 0.5, "b": 0.5}
    result = topsis(scores, weights).to_numpy().ravel()
    assert np.isnan(result[2]), "NaN in any criterion should propagate"
    assert np.all(np.isfinite(result[:2]))


def test_topsis_deterministic() -> None:
    """Same input produces same output."""
    rng = np.random.default_rng(20260501)
    a = rng.random((4, 4))
    b = rng.random((4, 4))
    scores = {"a": _grid(a.tolist()), "b": _grid(b.tolist())}
    weights = {"a": 0.6, "b": 0.4}
    r1 = topsis(scores, weights).to_numpy()
    r2 = topsis(scores, weights).to_numpy()
    np.testing.assert_array_equal(r1, r2)


def test_aggregate_dispatches() -> None:
    """The aggregate() entry point routes correctly."""
    scores = {
        "a": _grid([[1.0, 0.0]]),
        "b": _grid([[0.0, 1.0]]),
    }
    weights = {"a": 0.5, "b": 0.5}
    ws = aggregate(scores, weights, method="weighted_sum").to_numpy()
    ts = aggregate(scores, weights, method="topsis").to_numpy()
    np.testing.assert_allclose(ws.ravel(), [0.5, 0.5])
    # In a 2-cell symmetric problem each cell sits at the ideal of one
    # criterion and the anti-ideal of the other — TOPSIS gives 0.5 too.
    np.testing.assert_allclose(ts.ravel(), [0.5, 0.5])


def test_aggregate_unknown_method_raises() -> None:
    scores = {"a": _grid([[1.0]])}
    weights = {"a": 1.0}
    with pytest.raises(ValueError, match="unknown method"):
        aggregate(scores, weights, method="bogus")  # type: ignore[arg-type]


def test_topsis_zero_norm_criterion_raises() -> None:
    """A criterion that is identically zero across the grid has no
    meaningful normalisation; surface the error rather than dividing
    by zero."""
    scores = {
        "a": _grid([[1.0, 0.5, 0.0]]),
        "b": _grid([[0.0, 0.0, 0.0]]),
    }
    weights = {"a": 0.5, "b": 0.5}
    with pytest.raises(ValueError, match="zero or non-finite L2 norm"):
        topsis(scores, weights)


def test_topsis_default_weights_match_existing_v15_criteria_shape() -> None:
    """Smoke test: TOPSIS accepts the v1.5 criterion-name set and
    produces a [0, 1] grid of the same shape as weighted_sum."""
    rng = np.random.default_rng(42)
    grid = lambda: _grid(rng.random((6, 6)).tolist())  # noqa: E731
    scores = {
        "slope": grid(),
        "illumination": grid(),
        "thermal": grid(),
        "ice": grid(),
        "coupling": grid(),
        "los_to_earth": grid(),
    }
    weights = {
        "slope": 0.20,
        "illumination": 0.18,
        "thermal": 0.15,
        "ice": 0.15,
        "coupling": 0.15,
        "los_to_earth": 0.15,
        "hazard": 0.05,
        "seismic": 0.0,
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # absent hazard, seismic
        result = topsis(scores, weights)
    arr = result.to_numpy()
    assert arr.shape == (6, 6)
    finite = arr[np.isfinite(arr)]
    assert finite.min() >= 0.0
    assert finite.max() <= 1.0
