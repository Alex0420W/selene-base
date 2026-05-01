"""Tests for the v1.5 per-region tiled preprocess driver.

The tiled driver is exercised end-to-end against a small synthetic
elevation field on a local polar-stereographic grid. The synthetic
source covers all nine USGS polygon bounding boxes plus a 100 km
buffer, so :func:`run_tiled_per_region` can window each tile without
hitting the source edge. The horizon profile is computed on the CPU
(``use_gpu=False``) so the test runs on any host.

What the test does *not* try to verify:

- Numerical correctness of the GPU path (covered by
  ``tests/test_los_to_earth_gpu.py``).
- Agreement against the v1.4 240 m cached product (different grid
  origin → different pixel locations).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr

from selene_base.data.load import LUNAR_SOUTH_POLAR_CRS
from selene_base.pipeline import preprocess_tiled
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

POLAR_CRS = str(LUNAR_SOUTH_POLAR_CRS)


def _synthetic_lola_source(resolution_m: float = 240.0) -> xr.DataArray:
    """Build a synthetic LOLA-like raster covering all 9 USGS polygons + 100 km buffer.

    The output is a smooth bowl plus a Gaussian ridge so the horizon
    profile has structure to compare against. CRS is set to the v1.4
    south-polar stereographic so the tiled driver can reproject from it
    without a transform.
    """
    polygons = regions_polygons_to_geodataframe(target_crs=POLAR_CRS)
    bx_min, by_min, bx_max, by_max = polygons.total_bounds
    pad_m = 110_000.0  # > 100 km buffer used by the driver
    xmin = float(np.floor((bx_min - pad_m) / resolution_m) * resolution_m)
    ymin = float(np.floor((by_min - pad_m) / resolution_m) * resolution_m)
    xmax = float(np.ceil((bx_max + pad_m) / resolution_m) * resolution_m)
    ymax = float(np.ceil((by_max + pad_m) / resolution_m) * resolution_m)

    width = int(round((xmax - xmin) / resolution_m))
    height = int(round((ymax - ymin) / resolution_m))
    xs = np.linspace(xmin + resolution_m / 2, xmax - resolution_m / 2, width)
    ys = np.linspace(ymax - resolution_m / 2, ymin + resolution_m / 2, height)
    xx, yy = np.meshgrid(xs, ys)

    z = (xx**2 + yy**2) / 1e9  # gentle bowl, max ~hundreds of m
    z += 500.0 * np.exp(-((xx**2 + yy**2) / (2 * 50_000.0**2)))  # central peak

    da = xr.DataArray(
        z.astype(np.float32),
        dims=("y", "x"),
        coords={"y": ys, "x": xs},
        name="elevation_m",
    )
    return da.rio.write_crs(POLAR_CRS, inplace=False)


def test_compute_tile_specs_covers_nine_regions() -> None:
    specs = preprocess_tiled.compute_tile_specs(
        target_crs=POLAR_CRS, buffer_m=100_000.0, resolution_m=20.0
    )
    assert len(specs) == 9
    codes = {s.region_code for s in specs}
    assert codes == {"n2", "mm", "ma", "g2", "mp", "sp", "cb", "n1", "hw"}
    for s in specs:
        # Bounds snapped to 20 m grid
        for v in (s.xmin, s.ymin, s.xmax, s.ymax):
            assert (v % 20.0) == 0.0
        # Tile is at least 200 km wide once the 100 km buffer is added on each side.
        assert s.width_m >= 200_000.0 - 1.0
        assert s.height_m >= 200_000.0 - 1.0


def test_compute_tile_specs_filters_by_region_code() -> None:
    specs = preprocess_tiled.compute_tile_specs(
        target_crs=POLAR_CRS, buffer_m=100_000.0, resolution_m=20.0, region_codes=["SP"]
    )
    assert [s.region_code for s in specs] == ["sp"]
    assert specs[0].region_name == "Slater Plain"


def test_run_tiled_per_region_produces_valid_npz(tmp_path: Path) -> None:
    src = _synthetic_lola_source(resolution_m=240.0)
    src_tif = tmp_path / "raw" / "lola" / "synthetic.tif"
    src_tif.parent.mkdir(parents=True, exist_ok=True)
    src.rio.to_raster(src_tif, driver="GTiff", compress="DEFLATE")

    processed = tmp_path / "processed"
    results = preprocess_tiled.run_tiled_per_region(
        resolution_m=240.0,
        region_codes=["SP"],
        processed_dir=processed,
        source_path=src_tif,
        use_gpu=False,
        echo=lambda _msg: None,
    )

    assert len(results) == 1
    [r] = results
    assert r.region_code == "sp"
    assert r.status == "cached"
    assert r.output_path is not None and r.output_path.exists()

    npz = np.load(r.output_path, allow_pickle=False)
    horizon = npz["horizon_profile_deg"]
    assert horizon.ndim == 3
    assert horizon.shape[0] == 36  # default n_azimuths
    # Tile shape should match the spec.
    spec = preprocess_tiled.compute_tile_specs(
        target_crs=POLAR_CRS, buffer_m=100_000.0, resolution_m=240.0, region_codes=["SP"]
    )[0]
    expected_h, expected_w = spec.shape(240.0)
    assert horizon.shape[1:] == (expected_h, expected_w)
    assert np.isfinite(horizon).all(), "tiled horizon contains NaN"

    # The 100 km buffer guarantees that pixels at least 100 km from any
    # tile edge can ray-march the full default max_horizon_km without
    # leaving the tile, so the -90° initial sentinel must not survive
    # there. Pixels in the buffer ring may still carry -90° (rays going
    # outward escape the source) — that is documented behavior.
    interior_px = int(np.ceil(100_000.0 / 240.0))
    interior = horizon[:, interior_px:-interior_px, interior_px:-interior_px]
    assert interior.size > 0
    assert interior.max() <= 90.0
    assert interior.min() > -90.0, (
        "every interior pixel should observe at least one finite ray; "
        f"got min={interior.min():.3f}°"
    )
    # Curvature drop alone caps the horizon angle; on a smooth bowl +
    # central peak the apparent elevation should not exceed a few
    # degrees — well within a generous range.
    assert interior.min() >= -10.0
    assert interior.max() < 30.0

    # Cached metadata round-trips.
    assert int(npz["resolution_m"]) == 240
    assert npz["region_code"].item() == "sp"


def test_run_tiled_per_region_is_idempotent(tmp_path: Path) -> None:
    src = _synthetic_lola_source(resolution_m=240.0)
    src_tif = tmp_path / "raw" / "lola" / "synthetic.tif"
    src_tif.parent.mkdir(parents=True, exist_ok=True)
    src.rio.to_raster(src_tif, driver="GTiff", compress="DEFLATE")

    processed = tmp_path / "processed"
    kwargs = {
        "resolution_m": 240.0,
        "region_codes": ["SP"],
        "processed_dir": processed,
        "source_path": src_tif,
        "use_gpu": False,
        "echo": lambda _msg: None,
    }
    first = preprocess_tiled.run_tiled_per_region(**kwargs)
    second = preprocess_tiled.run_tiled_per_region(**kwargs)
    assert first[0].status == "cached"
    assert second[0].status == "cached"
    # Second run should not have re-derived (n_pixels is 0 for skip path).
    assert second[0].n_pixels == 0
    assert second[0].elapsed_s == pytest.approx(0.0)
