"""Tests for the v1.9 per-azimuth chunked + mmap-backed horizon-profile path.

The correctness contract is: chunked output must match single-pass output
within float32 precision. The horizon profile is the input to the LOS-to-
Earth criterion, which feeds the aggregate score and the per-region
ranking — any off-by-epsilon drift would silently bias the catalog.

These tests run on small synthetic elevation grids (no real LOLA data
needed). The 10 m end-to-end run on the actual GB10 + PGDA mosaic is
covered by the v1.9 release-time diagnostic on HW + SP + G2 (one-off,
~5–6 hours, results recorded in the v1.9 commit body and Roadmap entry).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr

from selene_base.criteria import los_to_earth
from selene_base.data.load import LUNAR_SOUTH_POLAR_CRS
from selene_base.pipeline import preprocess_tiled

POLAR_CRS = str(LUNAR_SOUTH_POLAR_CRS)


def _synthetic_tile(size: int = 64, resolution_m: float = 240.0) -> xr.DataArray:
    """Smooth bowl + central peak elevation tile suitable for ray-march tests.

    Small enough to run instantly, structured enough that the horizon
    profile has non-trivial variation per azimuth.
    """
    half = (size * resolution_m) / 2.0
    xs = np.linspace(-half + resolution_m / 2, half - resolution_m / 2, size)
    ys = np.linspace(half - resolution_m / 2, -half + resolution_m / 2, size)
    xx, yy = np.meshgrid(xs, ys)
    z = (xx**2 + yy**2) / 1e7
    z += 800.0 * np.exp(-((xx**2 + yy**2) / (2 * 5_000.0**2)))
    da = xr.DataArray(
        z.astype(np.float32),
        dims=("y", "x"),
        coords={"y": ys, "x": xs},
        name="elevation_m",
    )
    return da.rio.write_crs(POLAR_CRS, inplace=False)


def _run(elevation: xr.DataArray, **kwargs):
    return los_to_earth.derive_horizon_profile(
        elevation,
        n_azimuths=12,
        max_horizon_km=5.0,
        pixel_size_m=240.0,
        n_distance_samples=20,
        use_gpu=False,
        **kwargs,
    )


# ---- correctness: chunked == single-pass within float32 precision ---------


def test_chunked_output_matches_single_pass_chunk_size_3() -> None:
    elevation = _synthetic_tile()
    single = _run(elevation).to_numpy()
    chunked = _run(elevation, azimuth_chunk_size=3).to_numpy()
    np.testing.assert_array_equal(single, chunked)


def test_chunked_output_matches_single_pass_chunk_size_1() -> None:
    """Edge case: one azimuth per chunk — every iteration is its own chunk."""
    elevation = _synthetic_tile()
    single = _run(elevation).to_numpy()
    chunked = _run(elevation, azimuth_chunk_size=1).to_numpy()
    np.testing.assert_array_equal(single, chunked)


def test_chunked_output_matches_single_pass_chunk_size_equals_n_azimuths() -> None:
    """``chunk_size == n_azimuths`` is the unchunked path; result must be bit-identical."""
    elevation = _synthetic_tile()
    single = _run(elevation).to_numpy()
    chunked = _run(elevation, azimuth_chunk_size=12).to_numpy()
    np.testing.assert_array_equal(single, chunked)


def test_chunked_output_matches_single_pass_chunk_size_larger_than_n_azimuths() -> None:
    """``chunk_size > n_azimuths`` is clamped to ``n_azimuths``; same result."""
    elevation = _synthetic_tile()
    single = _run(elevation).to_numpy()
    chunked = _run(elevation, azimuth_chunk_size=99).to_numpy()
    np.testing.assert_array_equal(single, chunked)


# ---- mmap-backed accumulator -----------------------------------------------


def test_mmap_backed_out_param_writes_to_disk(tmp_path: Path) -> None:
    """``out=`` accepting a memmap writes results to the backing file.

    The load-bearing claim for the v1.9 architecture is: the on-disk file
    matches a single-pass in-memory run bit-for-bit. The DataArray's
    backing-array identity is implementation detail; only the disk content
    is contractual.
    """
    elevation = _synthetic_tile()
    n_az = 12
    h, w = elevation.shape
    npy_path = tmp_path / "horizon.npy"
    mmap = np.lib.format.open_memmap(npy_path, mode="w+", dtype=np.float32, shape=(n_az, h, w))

    _run(elevation, azimuth_chunk_size=3, out=mmap)
    mmap.flush()
    del mmap

    on_disk = np.load(npy_path, mmap_mode="r")
    single = _run(elevation).to_numpy()
    np.testing.assert_array_equal(np.asarray(on_disk), single)


def test_out_param_dtype_validated() -> None:
    elevation = _synthetic_tile()
    h, w = elevation.shape
    bad_out = np.zeros((12, h, w), dtype=np.float64)
    with pytest.raises(ValueError, match="out dtype must be float32"):
        _run(elevation, out=bad_out)


def test_out_param_shape_validated() -> None:
    elevation = _synthetic_tile()
    h, w = elevation.shape
    bad_out = np.zeros((11, h, w), dtype=np.float32)  # n_azimuths off by 1
    with pytest.raises(ValueError, match="out shape"):
        _run(elevation, out=bad_out)


def test_negative_azimuth_chunk_size_rejected() -> None:
    elevation = _synthetic_tile()
    with pytest.raises(ValueError, match="azimuth_chunk_size must be positive"):
        _run(elevation, azimuth_chunk_size=0)


# ---- default_azimuth_chunk_size heuristic ----------------------------------


def test_default_chunk_size_v1_5_resolutions_no_chunking() -> None:
    """At v1.5 tile sizes (≤ ~12000 px) the full buffer fits the budget."""
    # 240m: tile ≈ 2533×2533, full buffer ≈ 0.9 GB → no chunking
    assert (
        preprocess_tiled.default_azimuth_chunk_size(n_azimuths=36, height=2533, width=2533) is None
    )
    # 20m: tile ≈ 11000×11000, full buffer ≈ 17 GB → no chunking under 20 GB cap
    assert (
        preprocess_tiled.default_azimuth_chunk_size(n_azimuths=36, height=11000, width=11000)
        is None
    )


def test_default_chunk_size_10m_returns_chunked() -> None:
    """10m tile (~22000×22000) exceeds the 20 GB cap; expect a chunk size of 9."""
    chunk = preprocess_tiled.default_azimuth_chunk_size(n_azimuths=36, height=22000, width=22000)
    assert chunk is not None
    assert chunk == 9


def test_default_chunk_size_5m_returns_smaller_chunk() -> None:
    """5m tile (~44000×44000) needs further chunking; expect chunk == 2."""
    chunk = preprocess_tiled.default_azimuth_chunk_size(n_azimuths=36, height=44000, width=44000)
    assert chunk is not None
    assert chunk == 2


def test_default_chunk_size_respects_explicit_budget() -> None:
    """A 1 GB budget at v1.5's 20m tile size triggers chunking even there."""
    chunk = preprocess_tiled.default_azimuth_chunk_size(
        n_azimuths=36, height=11000, width=11000, budget_bytes=1024**3
    )
    assert chunk is not None
    assert 0 < chunk < 36


# ---- end-to-end: tiled driver chunked path produces identical NPZ ----------


def test_run_tiled_per_region_chunked_matches_single_pass(tmp_path: Path) -> None:
    """The driver's mmap-backed chunked path emits a ``.npy`` + ``.meta.npz``
    pair whose horizon array is identical to the single-pass ``.npz`` cache.

    The two paths intentionally use different on-disk formats — compressed
    single-file NPZ (v1.5) vs. mmap-loadable plain NPY plus a tiny metadata
    sidecar (v1.9, required at fine resolutions on unified-memory hardware
    so the rank stage can stream-load without OOM-ing). The bulk float32
    horizon array must still match bit-for-bit.
    """
    from tests.test_preprocess_tiled import _synthetic_lola_source

    src = _synthetic_lola_source(resolution_m=240.0)
    src_tif = tmp_path / "raw" / "lola" / "synthetic.tif"
    src_tif.parent.mkdir(parents=True, exist_ok=True)
    src.rio.to_raster(src_tif, driver="GTiff", compress="DEFLATE")

    single_dir = tmp_path / "processed_single"
    chunked_dir = tmp_path / "processed_chunked"

    common = dict(
        resolution_m=240.0,
        region_codes=["SP"],
        source_path=src_tif,
        use_gpu=False,
        echo=lambda _msg: None,
    )

    preprocess_tiled.run_tiled_per_region(processed_dir=single_dir, **common)
    preprocess_tiled.run_tiled_per_region(processed_dir=chunked_dir, azimuth_chunk_size=9, **common)

    single_npz = next(single_dir.glob("horizon_profile_southpole_*_sp.npz"))
    chunked_npy = next(chunked_dir.glob("horizon_profile_southpole_*_sp.npy"))
    chunked_meta = chunked_npy.with_suffix(".meta.npz")

    s = np.load(single_npz)
    c_horizon = np.load(chunked_npy, mmap_mode="r")
    c_meta = np.load(chunked_meta)
    np.testing.assert_array_equal(s["horizon_profile_deg"], np.asarray(c_horizon))
    np.testing.assert_array_equal(s["azimuth_deg"], c_meta["azimuth_deg"])
    # Chunked path emits .npy + .meta.npz only — no sibling single-file .npz.
    assert not list(chunked_dir.glob("horizon_profile_southpole_*_sp.npz"))
    # No leftover temp .npy files.
    assert not list(chunked_dir.glob("*.tmp.npy"))


# ---- compute_earth_visibility_fraction row-chunked path -------------------


def _libration_inputs(elevation: xr.DataArray):
    """Build minimal lat/lon/gamma DataArrays aligned with ``elevation``.

    The visibility math depends on these spatially-varying fields; their
    actual numerical values aren't load-bearing for the chunked-vs-unchunked
    equivalence check, only that they vary across the grid.
    """
    h, w = elevation.shape
    lat = xr.DataArray(
        np.linspace(-89.5, -88.0, h * w).reshape(h, w),
        dims=("y", "x"),
        coords={"y": elevation.coords["y"], "x": elevation.coords["x"]},
        name="lat",
    )
    lon = xr.DataArray(
        np.linspace(0.0, 359.0, h * w).reshape(h, w),
        dims=("y", "x"),
        coords={"y": elevation.coords["y"], "x": elevation.coords["x"]},
        name="lon",
    )
    gamma = xr.DataArray(
        np.linspace(-1.0, 1.0, h * w).reshape(h, w),
        dims=("y", "x"),
        coords={"y": elevation.coords["y"], "x": elevation.coords["x"]},
        name="gamma",
    )
    return lat, lon, gamma


def test_visibility_row_chunked_matches_unchunked() -> None:
    elevation = _synthetic_tile(size=64)
    horizon = _run(elevation)
    lat, lon, gamma = _libration_inputs(elevation)

    unchunked = los_to_earth.compute_earth_visibility_fraction(
        horizon, lat, lon, gamma, n_libration_samples=8
    )
    chunked = los_to_earth.compute_earth_visibility_fraction(
        horizon, lat, lon, gamma, n_libration_samples=8, row_chunk_size=16
    )
    np.testing.assert_array_equal(unchunked.to_numpy(), chunked.to_numpy())


def test_visibility_row_chunked_size_1_matches_unchunked() -> None:
    """Edge case: one row per chunk."""
    elevation = _synthetic_tile(size=32)
    horizon = _run(elevation)
    lat, lon, gamma = _libration_inputs(elevation)

    unchunked = los_to_earth.compute_earth_visibility_fraction(
        horizon, lat, lon, gamma, n_libration_samples=4
    )
    chunked = los_to_earth.compute_earth_visibility_fraction(
        horizon, lat, lon, gamma, n_libration_samples=4, row_chunk_size=1
    )
    np.testing.assert_array_equal(unchunked.to_numpy(), chunked.to_numpy())


def test_visibility_negative_row_chunk_rejected() -> None:
    elevation = _synthetic_tile(size=16)
    horizon = _run(elevation)
    lat, lon, gamma = _libration_inputs(elevation)

    with pytest.raises(ValueError, match="row_chunk_size must be positive"):
        los_to_earth.compute_earth_visibility_fraction(
            horizon, lat, lon, gamma, n_libration_samples=2, row_chunk_size=0
        )


def test_visibility_with_mmap_horizon_matches_in_memory(tmp_path: Path) -> None:
    """The mmap-loaded horizon (v1.9 cache shape) produces the same visibility
    fraction as a fully in-memory horizon — the rank-side load-bearing claim."""
    elevation = _synthetic_tile(size=32)
    horizon = _run(elevation)
    lat, lon, gamma = _libration_inputs(elevation)

    npy_path = tmp_path / "h.npy"
    np.save(npy_path, horizon.to_numpy())

    in_memory = los_to_earth.compute_earth_visibility_fraction(
        horizon, lat, lon, gamma, n_libration_samples=4
    )
    mmap_arr = np.load(npy_path, mmap_mode="r")
    horizon_mmap = xr.DataArray(
        mmap_arr,
        dims=horizon.dims,
        coords={d: horizon.coords[d] for d in horizon.dims if d in horizon.coords},
    )
    from_mmap = los_to_earth.compute_earth_visibility_fraction(
        horizon_mmap, lat, lon, gamma, n_libration_samples=4, row_chunk_size=8
    )
    np.testing.assert_array_equal(in_memory.to_numpy(), from_mmap.to_numpy())


# ---- preprocess + rank end-to-end via the new mmap cache ------------------


def test_existing_horizon_cache_path_prefers_npy(tmp_path: Path) -> None:
    """The cache resolver picks the v1.9 .npy + .meta.npz pair when both
    formats happen to be on disk — important on machines where rank can't
    safely materialise the legacy compressed .npz at high resolutions."""
    npz_path, npy_path, meta_path = preprocess_tiled.horizon_cache_paths(
        tmp_path, resolution_m=10.0, region_code="HW"
    )
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    npz_path.write_bytes(b"")
    npy_path.write_bytes(b"")
    meta_path.write_bytes(b"")

    chosen = preprocess_tiled.existing_horizon_cache_path(
        tmp_path, resolution_m=10.0, region_code="HW"
    )
    assert chosen == npy_path


def test_existing_horizon_cache_path_falls_back_to_npz(tmp_path: Path) -> None:
    """Legacy v1.5 caches (single ``.npz``) still resolve correctly."""
    npz_path, _, _ = preprocess_tiled.horizon_cache_paths(
        tmp_path, resolution_m=20.0, region_code="HW"
    )
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    npz_path.write_bytes(b"")

    chosen = preprocess_tiled.existing_horizon_cache_path(
        tmp_path, resolution_m=20.0, region_code="HW"
    )
    assert chosen == npz_path


def test_existing_horizon_cache_path_returns_none_when_missing(tmp_path: Path) -> None:
    chosen = preprocess_tiled.existing_horizon_cache_path(
        tmp_path, resolution_m=20.0, region_code="HW"
    )
    assert chosen is None


def test_existing_horizon_cache_path_orphan_npy_falls_through(tmp_path: Path) -> None:
    """A ``.npy`` without its sibling ``.meta.npz`` is incomplete; the
    resolver must not return it (would break the rank-side loader)."""
    _, npy_path, _ = preprocess_tiled.horizon_cache_paths(
        tmp_path, resolution_m=10.0, region_code="HW"
    )
    npy_path.parent.mkdir(parents=True, exist_ok=True)
    npy_path.write_bytes(b"")  # but no .meta.npz

    chosen = preprocess_tiled.existing_horizon_cache_path(
        tmp_path, resolution_m=10.0, region_code="HW"
    )
    assert chosen is None
