"""Combine per-criterion score grids into a single suitability map.

Two aggregators ship as of v1.7:

- :func:`weighted_sum` — the v1.4+ default, a renormalising weighted
  linear sum.
- :func:`topsis` — the v1.7 alternative, a Technique for Order
  Preference by Similarity to Ideal Solution (Hwang & Yoon 1981) that
  rewards cells close to the per-criterion ideal *and* far from the
  per-criterion anti-ideal. TOPSIS penalises lopsided profiles (a cell
  that saturates one criterion at 1.0 and scores 0.0 elsewhere ranks
  lower than a cell with even moderate-but-balanced scores) which the
  linear sum cannot.

Both share the same input/output contract: per-criterion score grids
on a common xarray-broadcastable grid plus per-criterion weights.
:func:`aggregate` is the public entry point and dispatches via the
``method`` argument; the default stays ``weighted_sum`` so v1.5/v1.6
callers see no behaviour change.

Reference (TOPSIS): Hwang, C.-L. & Yoon, K. (1981). *Multiple Attribute
Decision Making: Methods and Applications.* Lecture Notes in Economics
and Mathematical Systems 186, Springer.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from typing import Literal

import numpy as np
import xarray as xr

AggregateMethod = Literal["weighted_sum", "topsis"]


def weighted_sum(
    scores: Mapping[str, xr.DataArray],
    weights: Mapping[str, float],
) -> xr.DataArray:
    """Combine criterion score grids into a single [0, 1] suitability map.

    Weights are renormalised across the criteria provided in ``scores``,
    so callers can drop a criterion without rebalancing the weights file.
    If ``weights`` references criteria that are absent from ``scores``
    (the common week-N rollout case where only some criteria have shipped),
    a ``UserWarning`` is emitted naming the missing criteria and the
    remaining weights are renormalised to sum to 1. Inputs must broadcast
    to a common grid; xarray handles the alignment.

    Args:
        scores: Mapping from criterion name to a [0, 1] score DataArray.
            Must be non-empty.
        weights: Mapping from criterion name to a non-negative weight.
            Must contain a key for every name in ``scores``; may carry
            additional keys for criteria not yet implemented.

    Returns:
        DataArray of [0, 1] aggregate scores on the shared grid.

    Raises:
        ValueError: If ``scores`` is empty, if any weight is negative,
            or if the sum of weights for the provided criteria is zero.
        KeyError: If any criterion in ``scores`` lacks a matching weight.
    """
    if not scores:
        raise ValueError("scores must not be empty")

    missing_weights = [name for name in scores if name not in weights]
    if missing_weights:
        raise KeyError(f"missing weights for criteria: {sorted(missing_weights)}")

    selected_weights = {name: float(weights[name]) for name in scores}
    for name, w in selected_weights.items():
        if w < 0:
            raise ValueError(f"weight for {name!r} is negative: {w}")

    absent_criteria = [name for name in weights if name not in scores]
    if absent_criteria:
        warnings.warn(
            "weighted_sum: criteria with weights but no score grid will be "
            f"ignored: {sorted(absent_criteria)}. Remaining weights "
            f"renormalised across {sorted(scores)}.",
            stacklevel=2,
        )

    total = sum(selected_weights.values())
    if total == 0:
        raise ValueError("sum of weights for the provided criteria is zero")

    result: xr.DataArray | None = None
    for name, arr in scores.items():
        contribution = arr * (selected_weights[name] / total)
        result = contribution if result is None else result + contribution

    assert result is not None  # for type narrowing; scores is non-empty
    return result


def _validate_inputs_and_weights(
    scores: Mapping[str, xr.DataArray],
    weights: Mapping[str, float],
) -> dict[str, float]:
    """Shared validation + renormalisation of weights for any aggregator.

    Mirrors the contract of :func:`weighted_sum` exactly: missing
    weights raise, negative weights raise, missing-criteria-with-weights
    emit a warning and the remaining weights are renormalised across
    the criteria present in ``scores``.
    """
    if not scores:
        raise ValueError("scores must not be empty")

    missing_weights = [name for name in scores if name not in weights]
    if missing_weights:
        raise KeyError(f"missing weights for criteria: {sorted(missing_weights)}")

    selected_weights = {name: float(weights[name]) for name in scores}
    for name, w in selected_weights.items():
        if w < 0:
            raise ValueError(f"weight for {name!r} is negative: {w}")

    absent_criteria = [name for name in weights if name not in scores]
    if absent_criteria:
        warnings.warn(
            "aggregate: criteria with weights but no score grid will be "
            f"ignored: {sorted(absent_criteria)}. Remaining weights "
            f"renormalised across {sorted(scores)}.",
            stacklevel=3,
        )

    total = sum(selected_weights.values())
    if total == 0:
        raise ValueError("sum of weights for the provided criteria is zero")
    return {name: w / total for name, w in selected_weights.items()}


def topsis(
    scores: Mapping[str, xr.DataArray],
    weights: Mapping[str, float],
) -> xr.DataArray:
    """Combine criterion score grids via TOPSIS (Hwang & Yoon 1981).

    Steps applied per cell:

    1. **Vector normalise** each criterion grid to L2 unit length
       across all finite cells. This puts every criterion on the same
       scale before the weight is applied — without this step, a
       criterion whose raw values happen to be larger would dominate
       even at equal weight.
    2. **Apply weights** (renormalised over the criteria present in
       ``scores``, mirroring :func:`weighted_sum`'s contract).
    3. **Identify ideal and anti-ideal** per criterion: the maximum
       and minimum of the weighted-normalised grid. Every criterion in
       this codebase is "higher is better", so ideal is the max; if a
       cost-direction criterion is added later, it would need to be
       inverted before reaching this function.
    4. **Euclidean distance** from each cell to ideal (``d_ideal``)
       and anti-ideal (``d_anti``) across all criteria.
    5. **Score** = ``d_anti / (d_ideal + d_anti)``. By construction
       in [0, 1]: cells at the ideal score 1.0; cells at the
       anti-ideal score 0.0; balanced "good across the board" cells
       outscore lopsided ones because their distances to the
       anti-ideal point are larger in the high-dimensional sense.

    Cells where any criterion is NaN are propagated as NaN — same as
    :func:`weighted_sum`.

    Args:
        scores: Mapping from criterion name to a [0, 1] score DataArray
            (same convention as :func:`weighted_sum`).
        weights: Mapping from criterion name to a non-negative weight.

    Returns:
        DataArray of [0, 1] TOPSIS scores on the shared grid.

    Raises:
        ValueError: On empty input, negative weights, zero weight sum,
            or zero L2 norm in any criterion (would mean a criterion
            grid is identically zero, indeterminate normalisation).
        KeyError: On missing weights.
    """
    normalised_weights = _validate_inputs_and_weights(scores, weights)

    # Stack the per-criterion grids along a new dim so we can do the
    # L2 normalisation, weighting, and ideal/anti-ideal lookup in
    # vectorised xarray ops without materialising a Python loop over
    # criteria multiple times.
    names = list(scores)
    stacked = xr.concat([scores[n] for n in names], dim="criterion")
    stacked = stacked.assign_coords(criterion=names)

    # Step 1: vector normalise per criterion.
    sqsum = (stacked**2).sum(dim=stacked.dims[1:], skipna=True)
    norm = np.sqrt(sqsum)
    if (norm <= 0).any():
        zero_crits = [
            n for n, v in zip(names, norm.to_numpy(), strict=True) if not np.isfinite(v) or v <= 0
        ]
        raise ValueError(
            f"topsis: criteria with zero or non-finite L2 norm cannot be "
            f"normalised: {sorted(zero_crits)}"
        )
    normalised = stacked / norm

    # Step 2: apply per-criterion weights.
    weight_da = xr.DataArray(
        np.asarray([normalised_weights[n] for n in names], dtype=np.float64),
        dims=("criterion",),
        coords={"criterion": names},
    )
    weighted = normalised * weight_da

    # Step 3: ideal / anti-ideal across the grid (per criterion).
    ideal = weighted.max(dim=weighted.dims[1:], skipna=True)
    anti_ideal = weighted.min(dim=weighted.dims[1:], skipna=True)

    # Step 4: per-cell Euclidean distances. Sum-of-squares over the
    # criterion dim, then sqrt. NaN in any criterion at a cell yields
    # NaN distance for that cell, and then NaN score, matching the
    # weighted_sum NaN-propagation behaviour.
    d_ideal = np.sqrt(((weighted - ideal) ** 2).sum(dim="criterion", skipna=False))
    d_anti = np.sqrt(((weighted - anti_ideal) ** 2).sum(dim="criterion", skipna=False))

    # Step 5: closeness coefficient. Where the denominator is zero
    # (ideal == anti-ideal — i.e. a cell where all criteria are
    # identical to ideal AND anti-ideal, which only happens if all
    # criteria are uniform across the grid) the result is left as NaN
    # rather than producing 0/0; that case has no meaningful ordering
    # information anyway.
    denom = d_ideal + d_anti
    score = xr.where(denom > 0, d_anti / denom, np.nan)
    return score


def aggregate(
    scores: Mapping[str, xr.DataArray],
    weights: Mapping[str, float],
    *,
    method: AggregateMethod = "weighted_sum",
) -> xr.DataArray:
    """Combine per-criterion grids via the named aggregator.

    Public entry point that dispatches to :func:`weighted_sum` (default,
    v1.4+ behaviour) or :func:`topsis` (v1.7 alternative). Both share
    the same input/output contract; the choice changes how a balanced
    profile across criteria scores against a lopsided one (TOPSIS
    rewards balance; weighted-sum is indifferent to balance).

    Args:
        scores: Per-criterion [0, 1] score grids.
        weights: Per-criterion non-negative weights.
        method: ``"weighted_sum"`` (default) or ``"topsis"``.

    Returns:
        Aggregate score DataArray on the shared grid.
    """
    if method == "weighted_sum":
        return weighted_sum(scores, weights)
    if method == "topsis":
        return topsis(scores, weights)
    raise ValueError(f"aggregate: unknown method {method!r}; expected 'weighted_sum' or 'topsis'")
