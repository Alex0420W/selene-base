"""Ice criterion — rewards inferred water-ice resource potential.

Two interfaces, picked by which source data is available:

- :func:`compute` (the default) reads the Diviner PRP modeled
  ice-stability depth. ``ice_depth = 0`` means water ice is stable at
  the surface (best); intermediate depths indicate ice stable a few
  centimetres down; the source uses the sentinel ``-999`` for
  "effectively no ice" which the PDS4 parser converts to NaN.
- :func:`compute_from_lend` (formerly the default; kept as a drop-in
  for when the LEND CSETN flux map shows up) inverts a min-max
  normalised neutron flux and optionally adds a near-PSR proximity
  bonus.

The PSR mask used by either path is derived from the Mazarico
illumination raster via :func:`derive_psr_mask`. The PRP signal is a
more direct habitat-relevant ice proxy than LEND epithermal flux —
LEND measures total hydrogen including non-ice forms (OH, structural
water in regolith), whereas the PRP ice-depth field is a thermal
stability calculation specifically for water ice.

Filled in week 3; PRP-driven default added week 6.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from scipy.ndimage import distance_transform_edt

from selene_base.scoring.normalize import min_max

PRP_ICE_DEPTH_MAX_M = 2.87  # PRP cap; anything deeper was sentinel -999 -> NaN


def derive_psr_mask(
    illumination: xr.DataArray,
    *,
    threshold: float = 0.001,
) -> xr.DataArray:
    """Boolean mask of permanently shadowed regions.

    A PSR is a pixel that, over the 18.6-year illumination cycle,
    receives essentially no direct sunlight. We define it as
    ``illumination < threshold`` (the small positive cutoff lets the
    rasterisation noise floor pass without producing false positives).

    Args:
        illumination: DataArray of illumination fraction in [0, 1].
        threshold: Fraction below which a pixel counts as a PSR;
            must be in (0, 1).

    Returns:
        Boolean DataArray (True at PSR pixels) aligned with
        ``illumination``. NaN illumination becomes False.

    Raises:
        ValueError: If ``threshold`` is outside ``(0, 1)``.
    """
    if not (0 < threshold < 1):
        raise ValueError(f"threshold must be in (0, 1), got {threshold!r}")

    arr = illumination.to_numpy()
    mask = (arr < threshold) & np.isfinite(arr)

    out = xr.DataArray(
        mask,
        coords=illumination.coords,
        dims=illumination.dims,
        name="psr_mask",
    )
    if illumination.rio.crs is not None:
        out = out.rio.write_crs(illumination.rio.crs, inplace=False)
    return out


def compute(
    ice_depth_m: xr.DataArray,
    psr_mask: xr.DataArray | None = None,
    *,
    surface_ice_bonus: float = 0.5,
    near_psr_bonus: float = 0.2,
    near_psr_radius_km: float = 5.0,
    pixel_size_m: float = 240.0,
    max_depth_m: float = PRP_ICE_DEPTH_MAX_M,
) -> xr.DataArray:
    """Map Diviner PRP ice-stability depth to a [0, 1] resource score.

    Base score:

    - ``ice_depth == 0`` (stable at surface) → 1.0
    - ``0 < ice_depth ≤ max_depth_m`` → ``1 - ice_depth / max_depth_m``
    - ``ice_depth`` NaN (was sentinel ``-999``: effectively no ice) → 0.0

    Optional bonuses, summed and clipped to [0, 1]:

    - ``+ surface_ice_bonus`` where ``ice_depth == 0``.
    - ``+ near_psr_bonus`` where the cell is within
      ``near_psr_radius_km`` of any PSR pixel (when ``psr_mask`` is
      supplied).

    Args:
        ice_depth_m: DataArray of modeled water-ice stability depth
            (metres). NaN signals "no ice".
        psr_mask: Optional boolean DataArray on the same grid; True at
            PSR pixels (see :func:`derive_psr_mask`).
        surface_ice_bonus: Additive bonus applied where
            ``ice_depth == 0``. Must be in ``[0, 1]``.
        near_psr_bonus: Additive bonus applied within
            ``near_psr_radius_km`` of any PSR. Must be in ``[0, 1]``.
        near_psr_radius_km: Bonus-application radius around PSR pixels.
        pixel_size_m: Pixel edge length in metres; controls the
            distance-transform unit conversion.
        max_depth_m: Depth at which the linear ramp reaches zero.

    Returns:
        DataArray of [0, 1] scores aligned with ``ice_depth_m``.

    Raises:
        ValueError: On out-of-range parameters or shape mismatches.
    """
    if not (0 <= surface_ice_bonus <= 1):
        raise ValueError(f"surface_ice_bonus must be in [0, 1], got {surface_ice_bonus!r}")
    if not (0 <= near_psr_bonus <= 1):
        raise ValueError(f"near_psr_bonus must be in [0, 1], got {near_psr_bonus!r}")
    if near_psr_radius_km <= 0:
        raise ValueError(f"near_psr_radius_km must be positive, got {near_psr_radius_km!r}")
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")
    if max_depth_m <= 0:
        raise ValueError(f"max_depth_m must be positive, got {max_depth_m!r}")

    arr = ice_depth_m.to_numpy().astype(np.float64)
    no_ice = np.isnan(arr)
    safe_depth = np.where(no_ice, max_depth_m, arr)
    base = 1.0 - np.clip(safe_depth / max_depth_m, 0.0, 1.0)
    base = np.where(no_ice, 0.0, base)

    surface_mask = (~no_ice) & (np.abs(arr) < 1e-9)
    base = base + np.where(surface_mask, surface_ice_bonus, 0.0)

    if psr_mask is not None:
        mask_arr = np.asarray(psr_mask.to_numpy(), dtype=bool)
        if mask_arr.shape != base.shape:
            raise ValueError(
                f"psr_mask shape {mask_arr.shape!r} does not match ice_depth shape {base.shape!r}"
            )
        if mask_arr.any():
            radius_pixels = (near_psr_radius_km * 1000.0) / pixel_size_m
            distance_pixels = distance_transform_edt(~mask_arr)
            near = distance_pixels <= radius_pixels
            base = base + np.where(near, near_psr_bonus, 0.0)

    score = np.clip(base, 0.0, 1.0)

    out = xr.DataArray(
        score,
        coords=ice_depth_m.coords,
        dims=ice_depth_m.dims,
        name="ice_score",
    )
    if ice_depth_m.rio.crs is not None:
        out = out.rio.write_crs(ice_depth_m.rio.crs, inplace=False)
    return out


def compute_from_lend(
    epithermal_flux: xr.DataArray,
    psr_mask: xr.DataArray | None = None,
    *,
    near_psr_radius_km: float = 5.0,
    near_psr_bonus: float = 0.3,
    pixel_size_m: float = 240.0,
) -> xr.DataArray:
    """Map LEND epithermal-neutron flux to a [0, 1] ice-resource score.

    Base score: ``1 - min_max(flux)`` — lower flux is interpreted as
    more hydrogen and scored higher. When ``psr_mask`` is supplied, an
    additive bonus of ``near_psr_bonus`` is granted to pixels within
    ``near_psr_radius_km`` of any PSR pixel; the result is clipped to
    [0, 1].

    Kept as a drop-in for the day a south-polar LEND CSETN map shows
    up; the default :func:`compute` (PRP-based) is the better proxy
    for habitat-relevant ice.

    Args:
        epithermal_flux: DataArray of LEND neutron flux / count rate.
        psr_mask: Optional boolean DataArray on the same grid; True at
            PSR pixels.
        near_psr_radius_km: Radius (km) within which a pixel inherits
            the proximity bonus. Must be strictly positive.
        near_psr_bonus: Additive bonus applied inside the near-PSR
            neighbourhood. Must be in ``[0, 1]``.
        pixel_size_m: Pixel edge length in metres.

    Returns:
        DataArray of [0, 1] scores; NaN where flux is NaN.

    Raises:
        ValueError: On out-of-range parameters or shape mismatches.
    """
    if near_psr_radius_km <= 0:
        raise ValueError(f"near_psr_radius_km must be positive, got {near_psr_radius_km!r}")
    if not (0 <= near_psr_bonus <= 1):
        raise ValueError(f"near_psr_bonus must be in [0, 1], got {near_psr_bonus!r}")
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")

    flux = epithermal_flux.to_numpy().astype(np.float64)
    base = 1.0 - min_max(flux)

    if psr_mask is not None:
        mask_arr = np.asarray(psr_mask.to_numpy(), dtype=bool)
        if mask_arr.shape != base.shape:
            raise ValueError(
                f"psr_mask shape {mask_arr.shape!r} does not match flux shape {base.shape!r}"
            )
        if mask_arr.any():
            radius_pixels = (near_psr_radius_km * 1000.0) / pixel_size_m
            distance_pixels = distance_transform_edt(~mask_arr)
            near = distance_pixels <= radius_pixels
            base = base + np.where(near, near_psr_bonus, 0.0)
            base = np.clip(base, 0.0, 1.0)

    score = np.where(np.isnan(flux), np.nan, base)
    out = xr.DataArray(
        score,
        coords=epithermal_flux.coords,
        dims=epithermal_flux.dims,
        name="ice_score",
    )
    if epithermal_flux.rio.crs is not None:
        out = out.rio.write_crs(epithermal_flux.rio.crs, inplace=False)
    return out
