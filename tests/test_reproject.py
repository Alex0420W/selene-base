"""Tests for :mod:`selene_base.data.reproject`.

All tests are synthetic — no downloaded data required, so they run in
CI on every push.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS

from selene_base.data.reproject import (
    cache_processed,
    is_cog,
    reproject_to_grid,
)

LUNAR_GEOGRAPHIC = CRS.from_proj4("+proj=longlat +R=1737400 +no_defs +type=crs")
LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _make_lonlat_raster(
    values: np.ndarray, lon_extent: tuple[float, float], lat_extent: tuple[float, float]
) -> xr.DataArray:
    """Build a CRS-tagged DataArray on a regular lon/lat grid."""
    h, w = values.shape
    lon = np.linspace(lon_extent[0], lon_extent[1], w)
    lat = np.linspace(lat_extent[1], lat_extent[0], h)  # decreasing y
    da = xr.DataArray(values, dims=("y", "x"), coords={"y": lat, "x": lon})
    return da.rio.write_crs(LUNAR_GEOGRAPHIC, inplace=False)


class TestReprojectToGrid:
    def test_missing_crs_raises(self) -> None:
        da = xr.DataArray(np.zeros((4, 4)), dims=("y", "x"))
        with pytest.raises(ValueError, match="no CRS"):
            reproject_to_grid(
                da,
                target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
                bounds_m=(-1000, -1000, 1000, 1000),
                resolution_m=100,
            )

    def test_unknown_resampling_raises(self) -> None:
        da = _make_lonlat_raster(np.ones((10, 10)), (0, 10), (-90, -80))
        with pytest.raises(ValueError, match="unknown resampling"):
            reproject_to_grid(
                da,
                target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
                bounds_m=(-100_000, -100_000, 100_000, 100_000),
                resolution_m=10_000,
                resampling="not-a-method",
            )

    def test_non_positive_resolution_raises(self) -> None:
        da = _make_lonlat_raster(np.ones((4, 4)), (0, 1), (-90, -89))
        with pytest.raises(ValueError, match="resolution_m must be positive"):
            reproject_to_grid(
                da,
                target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
                bounds_m=(-1000, -1000, 1000, 1000),
                resolution_m=0,
            )

    def test_degenerate_bounds_raise(self) -> None:
        da = _make_lonlat_raster(np.ones((4, 4)), (0, 1), (-90, -89))
        with pytest.raises(ValueError, match="bounds_m"):
            reproject_to_grid(
                da,
                target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
                bounds_m=(1000, 1000, -1000, -1000),
                resolution_m=100,
            )

    def test_constant_field_remains_constant(self) -> None:
        # A uniform field reprojected to any compatible target stays uniform.
        da = _make_lonlat_raster(np.full((100, 100), 0.5), (0, 360), (-90, -80))
        out = reproject_to_grid(
            da,
            target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
            bounds_m=(-200_000, -200_000, 200_000, 200_000),
            resolution_m=10_000,
        )
        assert out.rio.crs == LUNAR_SOUTH_POLAR
        assert set(out.dims) == {"y", "x"}
        assert out.sizes["x"] == 40 and out.sizes["y"] == 40
        finite = out.to_numpy()[np.isfinite(out.to_numpy())]
        assert finite.size > 0
        np.testing.assert_allclose(finite, 0.5, atol=1e-6)

    def test_output_shape_matches_bounds_and_resolution(self) -> None:
        da = _make_lonlat_raster(np.zeros((50, 50)), (0, 360), (-90, -80))
        out = reproject_to_grid(
            da,
            target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
            bounds_m=(-50_000, -50_000, 50_000, 50_000),
            resolution_m=2_500,
        )
        assert out.sizes["x"] == 40
        assert out.sizes["y"] == 40

    def test_transform_is_valid(self) -> None:
        da = _make_lonlat_raster(np.zeros((20, 20)), (0, 360), (-90, -80))
        out = reproject_to_grid(
            da,
            target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
            bounds_m=(-100_000, -100_000, 100_000, 100_000),
            resolution_m=5_000,
        )
        t = out.rio.transform()
        assert t.a == 5_000  # pixel width in metres
        assert t.e == -5_000  # pixel height (negative because y decreases)


class TestCacheProcessed:
    def test_writes_cog_skips_when_present(self, tmp_path: Path) -> None:
        rng = np.random.default_rng(42)
        values = rng.uniform(size=(1024, 1024)).astype(np.float32)
        da = _make_lonlat_raster(values, (0, 360), (-90, -80))
        warped = reproject_to_grid(
            da,
            target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
            bounds_m=(-200_000, -200_000, 200_000, 200_000),
            resolution_m=400,
        )
        out = cache_processed(warped, "test_dataset", tmp_path)
        assert out.exists()
        assert out.name == "test_dataset_southpole_240m.tif"

        first_mtime = out.stat().st_mtime_ns
        again = cache_processed(warped, "test_dataset", tmp_path)
        assert again == out
        assert again.stat().st_mtime_ns == first_mtime  # skipped

    def test_overwrite_actually_overwrites(self, tmp_path: Path) -> None:
        da = _make_lonlat_raster(np.zeros((512, 512)), (0, 360), (-90, -80))
        warped = reproject_to_grid(
            da,
            target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
            bounds_m=(-50_000, -50_000, 50_000, 50_000),
            resolution_m=200,
        )
        out = cache_processed(warped, "ds", tmp_path)
        first_mtime = out.stat().st_mtime_ns
        again = cache_processed(warped, "ds", tmp_path, overwrite=True)
        assert again.stat().st_mtime_ns >= first_mtime  # rewritten

    def test_output_is_valid_cog_with_overviews(self, tmp_path: Path) -> None:
        # 1024×1024 large enough for COG driver to add overviews.
        rng = np.random.default_rng(0)
        values = rng.uniform(size=(1024, 1024)).astype(np.float32)
        da = _make_lonlat_raster(values, (0, 360), (-90, -80))
        warped = reproject_to_grid(
            da,
            target_crs=str(LUNAR_SOUTH_POLAR.to_proj4()),
            bounds_m=(-200_000, -200_000, 200_000, 200_000),
            resolution_m=400,
        )
        out = cache_processed(warped, "cog_test", tmp_path)

        with rasterio.open(out) as ds:
            assert ds.driver in {"COG", "GTiff"}
            block_h, block_w = ds.block_shapes[0]
            assert block_h >= 256 and block_w >= 256, (
                f"expected internal tiling >=256, got {ds.block_shapes[0]}"
            )
            assert ds.overviews(1), "expected at least one overview level"
        assert is_cog(out)

    def test_missing_crs_raises(self, tmp_path: Path) -> None:
        da = xr.DataArray(np.zeros((10, 10)), dims=("y", "x"))
        with pytest.raises(ValueError, match="no CRS"):
            cache_processed(da, "broken", tmp_path)
