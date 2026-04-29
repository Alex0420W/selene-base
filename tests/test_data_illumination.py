"""Smoke test for the Mazarico illumination raster."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from selene_base.data.load import load_illumination

ILLUM_PATH = Path("data/raw/illumination/avgvisib_65s_240m_201608.img")

pytestmark = pytest.mark.skipif(
    not ILLUM_PATH.exists(),
    reason=f"data not downloaded: {ILLUM_PATH}",
)


def test_returns_dataarray() -> None:
    da = load_illumination(ILLUM_PATH)
    assert isinstance(da, xr.DataArray)


def test_two_d_grid_with_yx_dims() -> None:
    da = load_illumination(ILLUM_PATH)
    assert set(da.dims) >= {"y", "x"}
    assert da.sizes["y"] > 0
    assert da.sizes["x"] > 0


def test_some_finite_pixels() -> None:
    da = load_illumination(ILLUM_PATH)
    sample = da.isel(y=slice(0, 256), x=slice(0, 256)).values
    assert np.isfinite(sample).any()


def test_values_are_fractions() -> None:
    da = load_illumination(ILLUM_PATH)
    sample = da.isel(y=slice(0, 256), x=slice(0, 256)).values
    finite = sample[np.isfinite(sample)]
    assert finite.size > 0
    assert finite.min() >= 0.0
    assert finite.max() <= 1.0 + 1e-6


def test_crs_is_set() -> None:
    da = load_illumination(ILLUM_PATH)
    assert da.rio.crs is not None
