"""GPU equivalence test for ``derive_horizon_profile``.

Skipped automatically when CuPy is not importable or no CUDA device is
visible, so the suite stays green on CPU-only CI. On GPU hosts the test
asserts that the CuPy-backed ray-march agrees with the NumPy/SciPy
reference to float32 precision (max |diff| <= 1e-3 deg) and that the
NaN mask is identical.

Locked to a small synthetic 64x64 grid so the test is seconds, not
minutes; the full 240 m equivalence (max |diff| ~= 2.3e-4 deg measured
on the v1.5 GPU port) is exercised by the v1.5 preprocess run, not here.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from selene_base.criteria.los_to_earth import derive_horizon_profile

cupy_unavailable = False
try:
    import cupy as cp

    if cp.cuda.runtime.getDeviceCount() == 0:
        cupy_unavailable = True
except Exception:  # ImportError, CUDARuntimeError on hosts without a driver
    cupy_unavailable = True

requires_gpu = pytest.mark.skipif(
    cupy_unavailable, reason="CuPy not installed or no CUDA device"
)


def _synthetic_elevation(height: int = 64, width: int = 64) -> xr.DataArray:
    rng = np.random.default_rng(20260430)
    z = rng.normal(0.0, 50.0, size=(height, width)).astype(np.float32)
    # A ridge so the horizon profile has structure to compare on.
    yy, xx = np.mgrid[0:height, 0:width]
    z += (300.0 * np.exp(-((xx - width / 2) ** 2) / (2 * 8**2))).astype(np.float32)
    return xr.DataArray(
        z,
        dims=("y", "x"),
        coords={"y": np.arange(height), "x": np.arange(width)},
        name="elevation_m",
    )


@requires_gpu
def test_gpu_matches_cpu_horizon_profile_to_float_precision() -> None:
    elevation = _synthetic_elevation()
    cpu = derive_horizon_profile(elevation, pixel_size_m=240.0).to_numpy()
    gpu = derive_horizon_profile(elevation, pixel_size_m=240.0, use_gpu=True).to_numpy()

    assert cpu.shape == gpu.shape
    assert np.array_equal(np.isnan(cpu), np.isnan(gpu))
    finite = np.isfinite(cpu) & np.isfinite(gpu)
    assert finite.any()
    diff = np.abs(cpu[finite] - gpu[finite])
    assert diff.max() < 1e-3, f"max|diff|={diff.max():.3e} deg exceeds 1e-3"


@requires_gpu
def test_gpu_horizon_profile_is_deterministic() -> None:
    elevation = _synthetic_elevation()
    a = derive_horizon_profile(elevation, pixel_size_m=240.0, use_gpu=True).to_numpy()
    b = derive_horizon_profile(elevation, pixel_size_m=240.0, use_gpu=True).to_numpy()
    assert np.array_equal(np.isnan(a), np.isnan(b))
    finite = np.isfinite(a) & np.isfinite(b)
    assert (a[finite] == b[finite]).all()
