"""Smoke test for the LOLA LDEM south-polar DEM."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from selene_base.data.load import load_lola_ldem

LOLA_PATH = Path("data/raw/lola/ldem_80s_80m.lbl")

pytestmark = pytest.mark.skipif(
    not LOLA_PATH.exists(),
    reason=f"data not downloaded: {LOLA_PATH}",
)


def test_returns_dataarray() -> None:
    da = load_lola_ldem(LOLA_PATH)
    assert isinstance(da, xr.DataArray)


def test_two_d_grid_with_yx_dims() -> None:
    da = load_lola_ldem(LOLA_PATH)
    assert set(da.dims) >= {"y", "x"}
    assert da.sizes["y"] > 0
    assert da.sizes["x"] > 0


def test_some_finite_pixels() -> None:
    da = load_lola_ldem(LOLA_PATH)
    sample = da.isel(y=slice(0, 256), x=slice(0, 256)).values
    assert np.isfinite(sample).any()


def test_elevation_in_plausible_range() -> None:
    da = load_lola_ldem(LOLA_PATH)
    sample = da.isel(y=slice(0, 256), x=slice(0, 256)).values
    finite = sample[np.isfinite(sample)]
    assert finite.size > 0
    assert finite.min() > -20_000
    assert finite.max() < 20_000


def test_crs_is_set() -> None:
    da = load_lola_ldem(LOLA_PATH)
    assert da.rio.crs is not None
