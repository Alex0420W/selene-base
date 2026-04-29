"""Smoke test for the Diviner Polar Resource Product loader.

Loads the cached PRP rasters (after ``selene preprocess``) and checks
that the temperature and ice-depth layers are physically sensible.
Skipped automatically when the cache files are not present so CI stays
green without 600 MB of source data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr

PROCESSED = Path("data/processed")
LAYER_PATHS = {
    name: PROCESSED / f"diviner_{name}_southpole_240m.tif"
    for name in ("temp_avg", "temp_max", "ice_depth")
}

pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in LAYER_PATHS.values()),
    reason=f"PRP cache not present: {sorted(p.as_posix() for p in LAYER_PATHS.values())}",
)


def _load(name: str) -> xr.DataArray:
    return rioxarray.open_rasterio(LAYER_PATHS[name], masked=True).squeeze("band", drop=True)


def test_temp_avg_in_plausible_range() -> None:
    da = _load("temp_avg")
    finite = da.to_numpy()[np.isfinite(da.to_numpy())]
    assert finite.size > 0
    # PSR cold-trap floors are ~25 K; sunlit rims peak below ~250 K for
    # an annual mean. Allow a generous range.
    assert 10.0 < finite.min() < 100.0
    assert 100.0 < finite.max() < 300.0


def test_temp_max_at_least_temp_avg() -> None:
    avg = _load("temp_avg").to_numpy()
    mx = _load("temp_max").to_numpy()
    mask = np.isfinite(avg) & np.isfinite(mx)
    assert mask.any()
    # Tmax should be >= Tavg almost everywhere; allow tiny interpolation
    # slack at mesh boundaries.
    assert (mx[mask] >= avg[mask] - 5.0).all()


def test_ice_depth_finite_or_nan() -> None:
    ice = _load("ice_depth").to_numpy()
    finite = ice[np.isfinite(ice)]
    if finite.size:
        # Triangle-mesh nearest-neighbour interpolation can pick up small
        # negative depths near mesh boundaries; accept everything from
        # ~-1 to the PRP cap of 2.87 m.
        assert finite.min() >= -2.0
        assert finite.max() <= 3.0


def test_layers_share_grid_shape() -> None:
    shapes = {name: _load(name).shape for name in LAYER_PATHS}
    assert len({tuple(s) for s in shapes.values()}) == 1, shapes
