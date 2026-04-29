"""Combine per-criterion score grids into a single suitability map.

Implements a weighted linear sum (weights renormalised over the criteria
actually present in the input). TOPSIS is planned for a follow-up but
deliberately not in the v0 scaffold.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping

import xarray as xr


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
