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

from selene_base.criteria.los_to_earth import (
    compute_earth_visibility_fraction,
    derive_horizon_profile,
)

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


def _synthetic_visibility_inputs() -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """Build a small horizon profile + lat/lon/gamma grid for the GPU equivalence test.

    The synthetic horizon has structure (some azimuths blocked, some
    open) so the libration sweep produces non-degenerate visibility
    fractions to compare on.
    """
    n_az, h, w = 36, 32, 32
    rng = np.random.default_rng(20260430)
    horizon_arr = (
        np.degrees(rng.normal(0.0, 0.05, size=(n_az, h, w)))
        - np.linspace(0.0, 1.5, n_az)[:, None, None]
    ).astype(np.float32)
    horizon = xr.DataArray(
        horizon_arr,
        dims=("azimuth", "y", "x"),
        coords={
            "azimuth": np.degrees(2.0 * np.pi * np.arange(n_az) / n_az),
            "y": np.arange(h),
            "x": np.arange(w),
        },
    )
    # Lat 80° S to 89° S, lon spread around the pole, gamma derived from synthetic projected coords.
    lat = np.linspace(-80.0, -89.0, h)[:, None] * np.ones(w)
    lon = np.linspace(0.0, 90.0, w)[None, :] * np.ones(h)[:, None]
    xs = np.linspace(-200_000.0, 200_000.0, w)
    ys = np.linspace(200_000.0, -200_000.0, h)
    xx, yy = np.meshgrid(xs, ys)
    gamma = np.arctan2(xx, yy)
    pixel_lat = xr.DataArray(lat, dims=("y", "x"), coords={"y": ys, "x": xs})
    pixel_lon = xr.DataArray(lon, dims=("y", "x"), coords={"y": ys, "x": xs})
    grid_conv = xr.DataArray(gamma, dims=("y", "x"), coords={"y": ys, "x": xs})
    return horizon, pixel_lat, pixel_lon, grid_conv


@requires_gpu
def test_gpu_matches_cpu_visibility_fraction() -> None:
    horizon, pixel_lat, pixel_lon, gamma = _synthetic_visibility_inputs()
    cpu = compute_earth_visibility_fraction(horizon, pixel_lat, pixel_lon, gamma).to_numpy()
    gpu = compute_earth_visibility_fraction(
        horizon, pixel_lat, pixel_lon, gamma, use_gpu=True
    ).to_numpy()
    assert cpu.shape == gpu.shape
    # Visibility fractions are exact rationals (count / 24); CPU and
    # GPU paths must agree to bitwise equality on finite cells.
    finite = np.isfinite(cpu) & np.isfinite(gpu)
    assert finite.any()
    assert np.array_equal(cpu[finite], gpu[finite])
