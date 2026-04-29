"""Smoke test for the Diviner annual Tbol max/min mosaics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from selene_base.data.load import load_diviner

TMAX_PATH = Path("data/raw/diviner/diviner_tbol_max_sp.tif")
TMIN_PATH = Path("data/raw/diviner/diviner_tbol_min_sp.tif")

pytestmark = pytest.mark.skipif(
    not (TMAX_PATH.exists() and TMIN_PATH.exists()),
    reason=f"data not downloaded: {TMAX_PATH}, {TMIN_PATH}",
)


def test_returns_dataset_with_both_vars() -> None:
    ds = load_diviner(TMAX_PATH, TMIN_PATH)
    assert isinstance(ds, xr.Dataset)
    assert "tbol_max" in ds.data_vars
    assert "tbol_min" in ds.data_vars


def test_temperatures_in_plausible_range() -> None:
    ds = load_diviner(TMAX_PATH, TMIN_PATH)
    for var in ("tbol_max", "tbol_min"):
        sample = ds[var].isel(y=slice(0, 256), x=slice(0, 256)).values
        finite = sample[np.isfinite(sample)]
        assert finite.size > 0, f"{var}: all pixels are NaN in sample"
        assert finite.min() > 10.0, f"{var}: min temp {finite.min()} K is implausible"
        assert finite.max() < 500.0, f"{var}: max temp {finite.max()} K is implausible"


def test_tmax_geq_tmin_on_overlap() -> None:
    ds = load_diviner(TMAX_PATH, TMIN_PATH)
    if ds["tbol_max"].shape != ds["tbol_min"].shape:
        pytest.skip("Tmax/Tmin grids differ in shape; defer comparison to reproject step")
    delta = (ds["tbol_max"] - ds["tbol_min"]).isel(y=slice(0, 256), x=slice(0, 256)).values
    finite = delta[np.isfinite(delta)]
    assert finite.size > 0
    assert (finite >= -1.0).all(), "Tmax should be >= Tmin everywhere"
