"""EVA-disc PSR access criterion (v2.0+).

NASA's documented walking-EVA radius is 2 km (A3GT priorities at LPSC
2025; EVA-EXP-0070 Rev D). Selene-base v1.x used a global coupling
product (PSR-distance × ridge-distance multiplicative); v2.0 replaces
this with per-cell EVA-disc integration that directly mirrors what
crews experience at a candidate landing site: the cold-class PSR area
that is reachable on foot from the landing pad.

Score: for each candidate cell, the fraction of cells within a 2 km
disc that satisfy the cold-class PSR criterion (Diviner PRP annual
maximum temperature below 110 K). Cells already in cold-class score
1.0; cells with no cold-class neighbour inside the EVA radius score
0.0. The output is in [0, 1] by construction.

The per-cell loop in the user-facing scoping note (cKDTree
``query_ball_point`` per cell) is replaced here with a single
``scipy.ndimage.convolve`` against a circular kernel — semantically
identical on the project's regular polar-stereographic grid, with the
edge-handling property the scoping note called out (cells near the
grid boundary score against the truncated set of neighbours actually
present rather than against an idealised disc area).

Cold-class threshold rationale: the 110 K cap is the standard water-ice
stability ceiling at the lunar surface (Paige et al. 2010, Science 330,
479; Hayne et al. 2015, Icarus 255, 58). Above 110 K, water ice
sublimates on geologic timescales; below 110 K it survives indefinitely
in the absence of energetic insolation events. Cold-class therefore
identifies the cells that *can* hold a working water-ice deposit, which
is the EVA target.

References:
- Lawrence, S. J. (2025). NTRS 20250008952. Site-selection FOM
  framework for crewed lunar surface missions.
- A3GT priorities at LPSC 2025 (Artemis III Geology Team).
- EVA-EXP-0070 Rev D. NASA EVA capability spec; documents the 2 km
  walking-EVA radius.
- Paige, D. A. et al. (2010). *Diviner Lunar Radiometer observations
  of cold traps in the Moon's south polar region.* Science 330, 479.
- Hayne, P. O. et al. (2015). *Evidence for exposed water ice in the
  Moon's south polar regions from LRO ultraviolet albedo and
  temperature measurements.* Icarus 255, 58.
- Wueller, F. et al. (2026), JGR Planets. Bundled DBF carries
  ``PSR_AREA`` / ``oldest_PSR`` per-site columns retained for
  cross-check; not consumed here.

Filled in v2.0.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from scipy import ndimage


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


def compute(
    prp_temp_max: xr.DataArray,
    *,
    eva_radius_km: float = 2.0,
    cold_threshold_k: float = 110.0,
    pixel_size_m: float = 240.0,
) -> xr.DataArray:
    """Per-cell EVA-disc PSR access score.

    For each cell on the input grid:

    .. code-block:: text

        cold_count   = number of cells within eva_radius_km whose
                       prp_temp_max < cold_threshold_k
        total_count  = number of finite-temperature cells within
                       eva_radius_km
        score        = cold_count / total_count    (or 0 if total = 0)

    Cells where ``prp_temp_max`` is itself NaN are emitted as NaN —
    matching the NaN-propagation contract used by every other criterion
    in the package. Edge cells whose disc is truncated by the grid
    boundary score against the cells actually present (no implicit
    extrapolation of "missing" disc area).

    Args:
        prp_temp_max: DataArray of Diviner PRP annual maximum surface
            temperature (K). NaN signals "no data at this cell".
        eva_radius_km: Walking-EVA radius. Defaults to 2.0 km, the
            value documented in EVA-EXP-0070 Rev D and reiterated in
            the A3GT LPSC 2025 priorities.
        cold_threshold_k: Maximum-temperature cap below which a cell
            is counted as cold-class. Defaults to 110 K, the
            water-ice stability ceiling. Lower values (e.g. 90 K,
            70 K) progressively restrict the count to colder traps.
        pixel_size_m: Edge length of one grid cell in metres.

    Returns:
        DataArray of [0, 1] scores aligned with ``prp_temp_max``,
        named ``eva_psr_access_score``. NaN at cells where the input
        is NaN.

    Raises:
        ValueError: On non-positive ``eva_radius_km``,
            ``cold_threshold_k``, or ``pixel_size_m``.
    """
    if eva_radius_km <= 0:
        raise ValueError(f"eva_radius_km must be positive, got {eva_radius_km!r}")
    if cold_threshold_k <= 0:
        raise ValueError(f"cold_threshold_k must be positive, got {cold_threshold_k!r}")
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")

    arr = prp_temp_max.to_numpy().astype(np.float64)
    valid_mask = np.isfinite(arr)
    cold_mask = valid_mask & (arr < cold_threshold_k)

    radius_pixels = eva_radius_km * 1000.0 / pixel_size_m
    kernel = _build_disc_kernel(radius_pixels)

    cold_count = ndimage.convolve(cold_mask.astype(np.float32), kernel, mode="constant", cval=0.0)
    total_count = ndimage.convolve(valid_mask.astype(np.float32), kernel, mode="constant", cval=0.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        score = np.where(total_count > 0, cold_count / total_count, 0.0)
    score = np.where(valid_mask, score, np.nan).astype(np.float64)
    score = np.clip(score, 0.0, 1.0)

    out = xr.DataArray(
        score,
        coords=prp_temp_max.coords,
        dims=prp_temp_max.dims,
        name="eva_psr_access_score",
    )
    if prp_temp_max.rio.crs is not None:
        out = out.rio.write_crs(prp_temp_max.rio.crs, inplace=False)
    return out
