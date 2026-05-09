"""Multi-class volatile access criterion (v2.1+).

Replaces the v1.x single ``ice`` criterion (Diviner PRP ice-depth +
PSR-proximity bonus combined into one [0, 1] score) with a three-class
volatile-access criterion that distinguishes thermal-class accessibility
inside a 2 km walking-EVA disc:

  - **H₂O cold trap**: ``temp_max < 110 K`` — the standard water-ice
    stability ceiling at the lunar surface (Paige et al. 2010, Hayne
    et al. 2015).
  - **CO₂ / NH₃ cold trap**: ``temp_max < 66 K`` — the freezing
    temperature for solid CO₂ and the upper bound for stable solid
    ammonia under polar conditions (Hayne et al. 2020 micro-cold-trap
    classification).
  - **Ultra-cold**: ``temp_max < 60 K`` — the deepest-trap regime
    where exotic volatiles (HCN, CH₃OH) survive on geologic
    timescales (Hayne et al. 2020).

The thresholds are nested (60 K ⊂ 66 K ⊂ 110 K), so by construction
``ultracold ≤ co2_nh3 ≤ h2o`` per cell. A site that accesses all
three classes within walking-EVA range is more ISRU-favorable than a
site that only accesses the H₂O class — and the criterion makes this
explicit by combining the three sub-scores rather than collapsing to
a single mask.

Score is the unweighted mean of three sub-scores, each defined as the
fraction of cells inside the 2 km EVA disc that satisfy that class's
threshold:

.. code-block:: text

    h2o_score        = fraction of disc cells with temp_max < 110 K
    co2_nh3_score    = fraction of disc cells with temp_max <  66 K
    ultracold_score  = fraction of disc cells with temp_max <  60 K
    score            = (h2o + co2_nh3 + ultracold) / 3

The equal weighting is deliberate: NASA's published priority ordering
across volatile classes is documented as a *priority list* (Lawrence
2025 FOM framework: "volatile diversity"), not as a quantitative
preference vector. Equal weights are agnostic to that ordering and
keep the criterion calibration-clean. A NASA-priority alternative
(e.g. ``(0.5, 0.3, 0.2)`` to up-weight H₂O) is exposed as
``sub_weights=`` for downstream sensitivity work.

Implementation mirrors v2.0's ``eva_psr_access`` exactly: a single
``scipy.ndimage.convolve`` per thermal-class mask against a circular
kernel built once. NaN propagation matches every other criterion in
the package — cells whose ``prp_temp_max`` is NaN are emitted as
NaN. Truncated discs at the grid boundary score against the cells
actually present (no implicit extrapolation).

References:
- Paige, D. A. et al. (2010). *Diviner Lunar Radiometer observations
  of cold traps in the Moon's south polar region.* Science 330, 479.
  ([doi:10.1126/science.1187726](https://doi.org/10.1126/science.1187726))
- Hayne, P. O. et al. (2015). *Evidence for exposed water ice in the
  Moon's south polar regions.* Icarus 255, 58.
- Hayne, P. O. et al. (2020). *Micro cold traps on the Moon.* Nature
  Astronomy 5, 169.
  ([doi:10.1038/s41550-020-1198-9](https://doi.org/10.1038/s41550-020-1198-9))
- Wueller, F. et al. (2026). JGR Planets — per-class volatile
  treatment in the peer-reviewed reference catalog.
- Lawrence, S. J. (2025). NTRS 20250008952 — site-selection FOM
  framework, including volatile-diversity FOM.

Filled in v2.1.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from scipy import ndimage

H2O_THRESHOLD_K = 110.0
CO2_NH3_THRESHOLD_K = 66.0
ULTRACOLD_THRESHOLD_K = 60.0
DEFAULT_SUB_WEIGHTS: tuple[float, float, float] = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)


def _build_disc_kernel(radius_pixels: float) -> np.ndarray:
    """Boolean disc kernel of a given radius in pixels.

    Pixels at integer offsets ``(dy, dx)`` are inside the disc when
    ``dy^2 + dx^2 <= r^2``. The kernel is square, sized to fit the
    disc tightly; the centre cell is always included.
    """
    if radius_pixels < 0:
        raise ValueError(f"radius_pixels must be non-negative, got {radius_pixels!r}")
    r_ceil = int(np.ceil(radius_pixels))
    if r_ceil == 0:
        return np.ones((1, 1), dtype=np.float32)
    yy, xx = np.ogrid[-r_ceil : r_ceil + 1, -r_ceil : r_ceil + 1]
    return (yy * yy + xx * xx <= radius_pixels * radius_pixels).astype(np.float32)


def compute_components(
    prp_temp_max: xr.DataArray,
    *,
    eva_radius_km: float = 2.0,
    h2o_threshold_k: float = H2O_THRESHOLD_K,
    co2_nh3_threshold_k: float = CO2_NH3_THRESHOLD_K,
    ultracold_threshold_k: float = ULTRACOLD_THRESHOLD_K,
    sub_weights: tuple[float, float, float] = DEFAULT_SUB_WEIGHTS,
    pixel_size_m: float = 240.0,
) -> dict[str, xr.DataArray]:
    """Per-cell three-class volatile-access score with diagnostic components.

    Returns a dict with four DataArrays aligned to ``prp_temp_max``:

    - ``"h2o_score"`` — fraction of EVA-disc cells with
      ``temp_max < h2o_threshold_k``.
    - ``"co2_nh3_score"`` — fraction with ``temp_max < co2_nh3_threshold_k``.
    - ``"ultracold_score"`` — fraction with ``temp_max < ultracold_threshold_k``.
    - ``"combined_score"`` — weighted mean of the three sub-scores,
      using ``sub_weights`` (defaults to equal 1/3 each).

    All four are in [0, 1] by construction. NaN propagates from
    ``prp_temp_max``: cells whose input is NaN emit NaN on every
    sub-score and on the combined score.

    Args:
        prp_temp_max: DataArray of Diviner PRP annual maximum surface
            temperature (K). NaN signals "no data at this cell".
        eva_radius_km: Walking-EVA radius. Defaults to 2.0 km
            (EVA-EXP-0070 Rev D + A3GT LPSC 2025).
        h2o_threshold_k: Upper bound on temp_max for the H₂O class.
            Default 110 K (water-ice stability ceiling).
        co2_nh3_threshold_k: Upper bound on temp_max for the
            CO₂ / NH₃ class. Default 66 K.
        ultracold_threshold_k: Upper bound for the ultra-cold class.
            Default 60 K.
        sub_weights: Three non-negative weights ``(w_h2o, w_co2_nh3,
            w_ultracold)`` whose sum is non-zero. Default
            ``(1/3, 1/3, 1/3)`` (equal). Renormalised internally so
            the combined score stays in [0, 1] for any weight vector
            whose sum is positive.
        pixel_size_m: Edge length of one grid cell in metres.

    Returns:
        Mapping ``{"h2o_score", "co2_nh3_score", "ultracold_score",
        "combined_score"}`` to DataArrays on the input grid.

    Raises:
        ValueError: On non-positive ``eva_radius_km`` /
            ``pixel_size_m``, on a non-monotonic threshold ordering,
            or on a non-positive sum of sub-weights.
    """
    if eva_radius_km <= 0:
        raise ValueError(f"eva_radius_km must be positive, got {eva_radius_km!r}")
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")
    if not (ultracold_threshold_k <= co2_nh3_threshold_k <= h2o_threshold_k):
        raise ValueError(
            "thresholds must satisfy ultracold ≤ co2_nh3 ≤ h2o; got "
            f"ultracold={ultracold_threshold_k!r}, co2_nh3={co2_nh3_threshold_k!r}, "
            f"h2o={h2o_threshold_k!r}"
        )
    if any(w < 0 for w in sub_weights):
        raise ValueError(f"sub_weights must be non-negative, got {sub_weights!r}")
    weight_sum = sum(sub_weights)
    if weight_sum <= 0:
        raise ValueError(f"sub_weights must have positive sum, got {sub_weights!r}")
    w_h2o, w_co2, w_uc = (w / weight_sum for w in sub_weights)

    arr = prp_temp_max.to_numpy().astype(np.float64)
    valid_mask = np.isfinite(arr)
    h2o_mask = valid_mask & (arr < h2o_threshold_k)
    co2_mask = valid_mask & (arr < co2_nh3_threshold_k)
    ultra_mask = valid_mask & (arr < ultracold_threshold_k)

    radius_pixels = eva_radius_km * 1000.0 / pixel_size_m
    kernel = _build_disc_kernel(radius_pixels)

    valid_count = ndimage.convolve(valid_mask.astype(np.float32), kernel, mode="constant", cval=0.0)
    h2o_count = ndimage.convolve(h2o_mask.astype(np.float32), kernel, mode="constant", cval=0.0)
    co2_count = ndimage.convolve(co2_mask.astype(np.float32), kernel, mode="constant", cval=0.0)
    ultra_count = ndimage.convolve(ultra_mask.astype(np.float32), kernel, mode="constant", cval=0.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        h2o_score = np.where(valid_count > 0, h2o_count / valid_count, 0.0)
        co2_score = np.where(valid_count > 0, co2_count / valid_count, 0.0)
        ultra_score = np.where(valid_count > 0, ultra_count / valid_count, 0.0)

    combined = w_h2o * h2o_score + w_co2 * co2_score + w_uc * ultra_score
    combined = np.clip(combined, 0.0, 1.0)

    nan_fill = np.where(valid_mask, 0.0, np.nan)
    h2o_score = np.where(valid_mask, h2o_score, np.nan) + nan_fill * 0.0
    co2_score = np.where(valid_mask, co2_score, np.nan) + nan_fill * 0.0
    ultra_score = np.where(valid_mask, ultra_score, np.nan) + nan_fill * 0.0
    combined = np.where(valid_mask, combined, np.nan)

    def _wrap(values: np.ndarray, name: str) -> xr.DataArray:
        out = xr.DataArray(
            values.astype(np.float64),
            coords=prp_temp_max.coords,
            dims=prp_temp_max.dims,
            name=name,
        )
        if prp_temp_max.rio.crs is not None:
            out = out.rio.write_crs(prp_temp_max.rio.crs, inplace=False)
        return out

    return {
        "h2o_score": _wrap(h2o_score, "multi_volatile_h2o_score"),
        "co2_nh3_score": _wrap(co2_score, "multi_volatile_co2_nh3_score"),
        "ultracold_score": _wrap(ultra_score, "multi_volatile_ultracold_score"),
        "combined_score": _wrap(combined, "multi_volatile_score"),
    }


def compute(
    prp_temp_max: xr.DataArray,
    *,
    eva_radius_km: float = 2.0,
    h2o_threshold_k: float = H2O_THRESHOLD_K,
    co2_nh3_threshold_k: float = CO2_NH3_THRESHOLD_K,
    ultracold_threshold_k: float = ULTRACOLD_THRESHOLD_K,
    sub_weights: tuple[float, float, float] = DEFAULT_SUB_WEIGHTS,
    pixel_size_m: float = 240.0,
) -> xr.DataArray:
    """Per-cell three-class volatile-access score (combined).

    Convenience wrapper over :func:`compute_components` that returns
    only the combined score. Use :func:`compute_components` when the
    per-class breakdown is needed for diagnostic reporting (e.g. "how
    many sites access all three thermal classes vs only H₂O?").

    See :func:`compute_components` for argument and threshold
    semantics.

    Returns:
        DataArray of [0, 1] scores aligned with ``prp_temp_max``,
        named ``multi_volatile_score``. NaN at cells where the input
        is NaN.
    """
    components = compute_components(
        prp_temp_max,
        eva_radius_km=eva_radius_km,
        h2o_threshold_k=h2o_threshold_k,
        co2_nh3_threshold_k=co2_nh3_threshold_k,
        ultracold_threshold_k=ultracold_threshold_k,
        sub_weights=sub_weights,
        pixel_size_m=pixel_size_m,
    )
    return components["combined_score"]
