"""Smoke test for the LEND south-polar epithermal-neutron map."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from selene_base.data.load import load_lend

LEND_PATH = Path("data/raw/lend/lend_csetn_sp.img")

pytestmark = pytest.mark.skipif(
    not LEND_PATH.exists(),
    reason=f"data not downloaded: {LEND_PATH}",
)


def test_returns_dataarray() -> None:
    da = load_lend(LEND_PATH)
    assert isinstance(da, xr.DataArray)


def test_two_d_grid_with_yx_dims() -> None:
    da = load_lend(LEND_PATH)
    assert set(da.dims) >= {"y", "x"}
    assert da.sizes["y"] > 0
    assert da.sizes["x"] > 0


def test_some_finite_pixels() -> None:
    da = load_lend(LEND_PATH)
    sample = da.isel(y=slice(0, 64), x=slice(0, 64)).values
    assert np.isfinite(sample).any()


def test_values_non_negative() -> None:
    da = load_lend(LEND_PATH)
    sample = da.isel(y=slice(0, 64), x=slice(0, 64)).values
    finite = sample[np.isfinite(sample)]
    assert finite.size > 0
    assert (finite >= 0).all(), "neutron counts/flux should be non-negative"
