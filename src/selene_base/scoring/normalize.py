"""Pure-numpy [0, 1] normalisation primitives.

Three families of mapping cover every criterion in the pipeline:

- :func:`min_max` — linearly rescale a value range to [0, 1].
- :func:`optimal_range` — Gaussian peak around a target value
  (e.g. "moderate temperature is best").
- :func:`inverse_threshold` — score 1 at zero, decaying linearly to 0
  at a hard cutoff (e.g. "less slope is better, anything above 10° is
  unusable").

NaNs propagate through every function unchanged so that masked pixels
stay masked downstream. These are the only real implementations in the
scaffold; the criteria modules build on top of them.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def min_max(x: ArrayLike) -> NDArray[np.float64]:
    """Linearly rescale ``x`` so its valid range maps to [0, 1].

    NaNs in the input remain NaN in the output. If every value is equal
    (or the array is all-NaN), the result is all-zeros for the finite
    cells (NaN is preserved); this avoids divide-by-zero while keeping
    the output bounded.

    Args:
        x: Array-like of numeric values.

    Returns:
        Array of the same shape as ``x`` with values in [0, 1] and NaNs
        preserved.
    """
    arr = np.asarray(x, dtype=np.float64)
    valid = ~np.isnan(arr)
    if not valid.any():
        return arr.copy()
    x_min = float(np.nanmin(arr))
    x_max = float(np.nanmax(arr))
    if x_max == x_min:
        return np.where(np.isnan(arr), np.nan, 0.0)
    return (arr - x_min) / (x_max - x_min)


def optimal_range(x: ArrayLike, target: float, sigma: float) -> NDArray[np.float64]:
    """Score ``x`` by Gaussian distance from a target value.

    The score is ``exp(-(x - target)**2 / (2 * sigma**2))``: it peaks at
    1 when ``x == target`` and decays toward 0 as ``|x - target|`` grows.

    Args:
        x: Array-like of numeric values.
        target: Value that should receive the maximum score of 1.
        sigma: Gaussian width (in the same units as ``x``); must be
            strictly positive.

    Returns:
        Array of the same shape as ``x`` with values in [0, 1] and NaNs
        preserved.

    Raises:
        ValueError: If ``sigma`` is not strictly positive.
    """
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError(f"sigma must be a positive finite number, got {sigma!r}")
    arr = np.asarray(x, dtype=np.float64)
    return np.exp(-((arr - target) ** 2) / (2.0 * sigma**2))


def inverse_threshold(x: ArrayLike, threshold: float) -> NDArray[np.float64]:
    """Score ``x`` so that zero is best and ``threshold`` is unusable.

    Returns ``1 - x / threshold`` clipped into [0, 1]: a value of 0 maps
    to 1, ``threshold`` maps to 0, and anything above ``threshold`` is
    pinned at 0. Negative values are pinned at 1 (treated as "even
    better than zero"). NaNs are preserved.

    Args:
        x: Array-like of non-negative numeric values.
        threshold: Cutoff at which the score drops to 0; must be
            strictly positive.

    Returns:
        Array of the same shape as ``x`` with values in [0, 1] and NaNs
        preserved.

    Raises:
        ValueError: If ``threshold`` is not strictly positive.
    """
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError(f"threshold must be a positive finite number, got {threshold!r}")
    arr = np.asarray(x, dtype=np.float64)
    raw = 1.0 - (arr / threshold)
    return np.clip(raw, 0.0, 1.0)
