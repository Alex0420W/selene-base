"""Ice criterion — rewards inferred water-ice resource potential.

LEND epithermal-neutron suppression maps trace hydrogen abundance,
which under polar conditions is widely interpreted as buried water ice.
The base score is one minus the min-max-normalised neutron flux (lower
flux → more H → higher score). Optionally, sites within
``near_psr_radius_km`` of any permanently-shadowed-region (PSR) pixel
get an additive bonus; PSRs themselves are derived from the Mazarico
illumination map via :func:`derive_psr_mask`.

Filled in week 3.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from scipy.ndimage import distance_transform_edt

from selene_base.scoring.normalize import min_max


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
        ``illumination``. NaN illumination becomes False (unknown is
        not a PSR).

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
    epithermal_flux: xr.DataArray,
    psr_mask: xr.DataArray | None = None,
    *,
    near_psr_radius_km: float = 5.0,
    near_psr_bonus: float = 0.3,
    pixel_size_m: float = 240.0,
) -> xr.DataArray:
    """Map epithermal neutron flux (and optional PSR proximity) to a [0, 1] score.

    Base score: ``1 - min_max(flux)`` — lower flux is interpreted as
    more hydrogen and scored higher. When ``psr_mask`` is supplied, an
    additive bonus of ``near_psr_bonus`` is granted to pixels within
    ``near_psr_radius_km`` of any PSR pixel; the result is clipped to
    [0, 1]. With ``psr_mask=None`` only the base score is returned.

    Args:
        epithermal_flux: DataArray of LEND neutron flux / count rate.
        psr_mask: Optional boolean DataArray on the same grid; True at
            PSR pixels (see :func:`derive_psr_mask`).
        near_psr_radius_km: Radius (km) within which a pixel inherits
            the proximity bonus. Must be strictly positive.
        near_psr_bonus: Additive bonus (in score units) applied inside
            the near-PSR neighbourhood. Must be in [0, 1].
        pixel_size_m: Pixel edge length in metres; used for the
            distance transform that grows the PSR mask. Defaults to
            the project's 240 m grid.

    Returns:
        DataArray of [0, 1] scores; NaN where flux is NaN.

    Raises:
        ValueError: On bad parameter values (radius non-positive, bonus
            outside [0, 1], pixel_size non-positive, or shape mismatch
            between flux and PSR mask).
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
