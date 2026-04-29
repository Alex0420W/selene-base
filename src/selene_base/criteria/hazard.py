"""Hazard criterion — penalises proximity to fresh craters.

Operates on a precomputed crater-density grid produced by
:func:`selene_base.data.rasterize.rasterize_crater_density`. Score
saturates at zero once the local crater count crosses a user-set
threshold, modelling impact ejecta, rim instability, and trafficability
in a single composite.

Filled in week 3.
"""

from __future__ import annotations

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr


def compute(
    crater_density: xr.DataArray,
    *,
    saturation_count: float = 50.0,
) -> xr.DataArray:
    """Map crater density to a [0, 1] safety score.

    ``score = clip(1 - density / saturation_count, 0, 1)``: zero
    craters within the radius the density was computed over scores 1.0,
    ``saturation_count`` craters scores 0.0, and anything above
    saturation is clipped to 0.0.

    Args:
        crater_density: DataArray of crater counts per pixel (float
            values OK; produced by ``rasterize_crater_density``).
        saturation_count: Density (in the same units as the input)
            mapping to score 0.0; must be strictly positive.

    Returns:
        DataArray of [0, 1] scores aligned with ``crater_density``;
        NaN where the input is NaN.

    Raises:
        ValueError: If ``saturation_count`` is non-positive.
    """
    if saturation_count <= 0:
        raise ValueError(f"saturation_count must be positive, got {saturation_count!r}")

    arr = crater_density.to_numpy().astype(np.float64)
    raw = 1.0 - (arr / saturation_count)
    score = np.clip(raw, 0.0, 1.0)

    out = xr.DataArray(
        score,
        coords=crater_density.coords,
        dims=crater_density.dims,
        name="hazard_score",
    )
    if crater_density.rio.crs is not None:
        out = out.rio.write_crs(crater_density.rio.crs, inplace=False)
    return out
