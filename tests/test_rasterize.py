"""Tests for :mod:`selene_base.data.rasterize`."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import Point

from selene_base.data.rasterize import rasterize_crater_density

LUNAR_GEOGRAPHIC = CRS.from_proj4("+proj=longlat +R=1737400 +no_defs +type=crs")
LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _polar_grid(width: int, height: int, pixel_m: float) -> xr.DataArray:
    """A blank target grid in lunar south polar stereographic, projected metres.

    Origin at (0, 0); +x east, +y north (so pixel (0, 0) is at xmin/ymax).
    """
    half_x = (width / 2.0) * pixel_m
    half_y = (height / 2.0) * pixel_m
    da = xr.DataArray(
        np.zeros((height, width), dtype=np.float32),
        dims=("y", "x"),
        coords={
            "y": np.linspace(half_y - pixel_m / 2, -half_y + pixel_m / 2, height),
            "x": np.linspace(-half_x + pixel_m / 2, half_x - pixel_m / 2, width),
        },
    )
    return da.rio.write_crs(LUNAR_SOUTH_POLAR, inplace=False)


def _craters_at_polar_xy(coords_m: list[tuple[float, float]]) -> gpd.GeoDataFrame:
    """Build a GeoDataFrame from polar-xy coords by reprojecting to lon/lat."""
    from pyproj import Transformer

    t = Transformer.from_crs(LUNAR_SOUTH_POLAR, LUNAR_GEOGRAPHIC, always_xy=True)
    pts: list[Point] = []
    for x, y in coords_m:
        lon, lat = t.transform(x, y)
        pts.append(Point(lon, lat))
    return gpd.GeoDataFrame({"id": range(len(pts))}, geometry=pts, crs=LUNAR_GEOGRAPHIC)


class TestRasterizeCraterDensity:
    def test_empty_catalog_returns_zero_grid(self) -> None:
        grid = _polar_grid(8, 8, pixel_m=1000.0)
        empty = gpd.GeoDataFrame({"id": []}, geometry=[], crs=LUNAR_GEOGRAPHIC)
        out = rasterize_crater_density(empty, grid, radius_km=10.0)
        assert out.shape == (8, 8)
        np.testing.assert_array_equal(out.to_numpy(), 0.0)

    def test_single_crater_centre_pixel_count_is_one(self) -> None:
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        # One crater at the polar origin (0, 0). Search radius 0.1 km should
        # only catch the centre pixel.
        craters = _craters_at_polar_xy([(0.0, 0.0)])
        out = rasterize_crater_density(craters, grid, radius_km=0.1, chunk_rows=4).to_numpy()
        assert out.sum() == 1, f"expected exactly one positive pixel, got {out.sum()}"
        # The grid is 11x11 so the centre is row=5, col=5.
        assert out[5, 5] == 1

    def test_radius_grows_count(self) -> None:
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        craters = _craters_at_polar_xy([(0.0, 0.0)])
        small = rasterize_crater_density(craters, grid, radius_km=0.5).to_numpy()
        large = rasterize_crater_density(craters, grid, radius_km=3.0).to_numpy()
        assert large.sum() > small.sum()
        # Every pixel should have count 0 or 1 since there's only one crater.
        assert set(np.unique(small).tolist()).issubset({0, 1})
        assert set(np.unique(large).tolist()).issubset({0, 1})

    def test_overlapping_craters_count_each(self) -> None:
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        # Three craters within 0.5 km of (0, 0): centre pixel should see all three
        craters = _craters_at_polar_xy([(0.0, 0.0), (200.0, 0.0), (-200.0, 0.0)])
        out = rasterize_crater_density(craters, grid, radius_km=1.0).to_numpy()
        assert out[5, 5] == 3

    def test_diameter_filter_drops_small_craters(self) -> None:
        grid = _polar_grid(5, 5, pixel_m=1000.0)
        from pyproj import Transformer

        t = Transformer.from_crs(LUNAR_SOUTH_POLAR, LUNAR_GEOGRAPHIC, always_xy=True)
        rows = []
        for d in [0.5, 1.5, 3.0]:
            lon, lat = t.transform(0.0, 0.0)
            rows.append({"diam_km": d, "geometry": Point(lon, lat)})
        gdf = gpd.GeoDataFrame(rows, crs=LUNAR_GEOGRAPHIC)
        out = rasterize_crater_density(gdf, grid, radius_km=10.0, diameter_col="diam_km").to_numpy()
        # Only the two craters with diameter ≥ 1 km contribute.
        np.testing.assert_array_equal(out, 2.0)

    def test_non_positive_radius_rejected(self) -> None:
        grid = _polar_grid(4, 4, pixel_m=1000.0)
        craters = _craters_at_polar_xy([(0.0, 0.0)])
        with pytest.raises(ValueError, match="radius_km"):
            rasterize_crater_density(craters, grid, radius_km=0.0)

    def test_missing_grid_crs_rejected(self) -> None:
        bad = xr.DataArray(np.zeros((4, 4)), dims=("y", "x"))
        craters = _craters_at_polar_xy([(0.0, 0.0)])
        with pytest.raises(ValueError, match="target_grid has no CRS"):
            rasterize_crater_density(craters, bad)

    def test_missing_crater_crs_rejected(self) -> None:
        grid = _polar_grid(4, 4, pixel_m=1000.0)
        bad = gpd.GeoDataFrame({"id": [0]}, geometry=[Point(0, 0)])
        with pytest.raises(ValueError, match="craters has no CRS"):
            rasterize_crater_density(bad, grid)
