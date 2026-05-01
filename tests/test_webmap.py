"""Tests for :mod:`selene_base.viz.webmap`."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.nasa_regions import regions_to_geodataframe
from selene_base.viz.webmap import build_map

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _polar_score_cog(tmp_path: Path) -> Path:
    rng = np.random.default_rng(0)
    h = w = 256
    pixel_m = 1000.0
    half_x = (w / 2.0) * pixel_m
    half_y = (h / 2.0) * pixel_m
    da = xr.DataArray(
        rng.uniform(0, 1, size=(h, w)).astype(np.float32),
        dims=("y", "x"),
        coords={
            "y": np.linspace(half_y - pixel_m / 2, -half_y + pixel_m / 2, h),
            "x": np.linspace(-half_x + pixel_m / 2, half_x - pixel_m / 2, w),
        },
    ).rio.write_crs(LUNAR_SOUTH_POLAR, inplace=False)
    out = tmp_path / "score.tif"
    da.rio.to_raster(
        out,
        driver="COG",
        compress="DEFLATE",
        BLOCKSIZE=256,
        OVERVIEWS="AUTO",
        BIGTIFF="IF_SAFER",
    )
    return out


def _three_top_sites() -> gpd.GeoDataFrame:
    rows = [
        {
            "site_id": f"site_{i + 1:02d}",
            "rank": i + 1,
            "score": 0.9 - 0.05 * i,
            "lat": -85.0 - i * 0.5,
            "lon": -10.0 + i * 5.0,
            "x_m": 0.0,
            "y_m": 0.0,
            "score_slope": 0.7,
            "score_illumination": 0.8,
            "score_thermal": float("nan"),
            "score_ice": float("nan"),
            "score_hazard": 0.95,
            "score_seismic": float("nan"),
        }
        for i in range(3)
    ]
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


def test_build_map_writes_html(tmp_path: Path) -> None:
    score_cog = _polar_score_cog(tmp_path)
    sites = _three_top_sites()
    nasa = regions_to_geodataframe()
    out = tmp_path / "webmap.html"
    result_path = build_map(score_cog, sites, nasa, out)
    assert result_path == out
    assert out.exists()
    # Non-trivial HTML — large enough to have markers, polygons, embedded
    # PNG overlay, and per-site popup tables. The 256×256 score COG yields
    # roughly 80–120 KB on this fixture; assert generously.
    assert out.stat().st_size > 50_000


def test_html_contains_layer_names_and_per_site_markers(tmp_path: Path) -> None:
    score_cog = _polar_score_cog(tmp_path)
    sites = _three_top_sites()
    nasa = regions_to_geodataframe()
    out = tmp_path / "webmap.html"
    build_map(score_cog, sites, nasa, out)
    html = out.read_text(encoding="utf-8")
    assert "Aggregate suitability score" in html
    assert "NASA Artemis IV" in html
    assert "selene-base top sites" in html
    for site_id in sites["site_id"]:
        assert site_id in html
    for region_name in nasa["name"]:
        assert region_name in html


def test_skip_per_criterion_layers_when_processed_dir_missing(tmp_path: Path) -> None:
    score_cog = _polar_score_cog(tmp_path)
    sites = _three_top_sites()
    nasa = regions_to_geodataframe()
    out = tmp_path / "webmap.html"
    build_map(score_cog, sites, nasa, out, processed_dir=tmp_path / "nonexistent")
    html = out.read_text(encoding="utf-8")
    assert "Aggregate suitability score" in html
    # Per-criterion overlays would announce "slope score" / "hazard score"; absent here.
    assert "slope score" not in html.lower()


def test_missing_score_cog_raises(tmp_path: Path) -> None:
    sites = _three_top_sites()
    nasa = regions_to_geodataframe()
    with pytest.raises((FileNotFoundError, OSError)):
        build_map(tmp_path / "nope.tif", sites, nasa, tmp_path / "webmap.html")
