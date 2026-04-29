"""Tests for :mod:`selene_base.viz.site_report`."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.nasa_regions import regions_to_geodataframe
from selene_base.viz.site_report import (
    generate_site_index,
    generate_site_report,
)

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _score_cog(tmp_path: Path) -> Path:
    rng = np.random.default_rng(1)
    h = w = 200
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
    da.rio.to_raster(out, driver="GTiff", compress="DEFLATE")
    return out


def _site_row() -> gpd.GeoSeries:
    record = {
        "site_id": "site_07",
        "rank": 7,
        "score": 0.84,
        "lat": -86.5,
        "lon": -5.0,
        "x_m": 0.0,
        "y_m": 0.0,
        "score_slope": 0.92,
        "score_illumination": 0.78,
        "score_thermal": float("nan"),
        "score_ice": float("nan"),
        "score_hazard": 0.80,
        "score_seismic": float("nan"),
    }
    geom = Point(record["lon"], record["lat"])
    gdf = gpd.GeoDataFrame([record], geometry=[geom], crs=LUNAR_GEOGRAPHIC_CRS)
    return gdf.iloc[0]


def test_report_writes_html_with_site_id_in_title(tmp_path: Path) -> None:
    score_cog = _score_cog(tmp_path)
    nasa = regions_to_geodataframe()
    site = _site_row()
    out_dir = tmp_path / "sites"
    path = generate_site_report(site, score_cog, out_dir, nasa_regions=nasa)
    assert path.exists()
    assert path.name == "site_07.html"
    html = path.read_text(encoding="utf-8")
    assert "site_07" in html
    assert "rank 7" in html
    assert "0.840" in html  # aggregate score formatted
    # Embedded base64 payloads (bar chart + mini map) should be present.
    assert html.count("data:image/png;base64,") >= 2


def test_report_lists_dominant_criterion_and_nearest_nasa(tmp_path: Path) -> None:
    score_cog = _score_cog(tmp_path)
    nasa = regions_to_geodataframe()
    site = _site_row()
    path = generate_site_report(site, score_cog, tmp_path / "sites", nasa_regions=nasa)
    html = path.read_text(encoding="utf-8")
    assert "slope" in html
    assert "Haworth" in html  # site is at Haworth's centroid


def test_report_runs_without_nasa_regions(tmp_path: Path) -> None:
    score_cog = _score_cog(tmp_path)
    site = _site_row()
    path = generate_site_report(site, score_cog, tmp_path / "sites", nasa_regions=None)
    html = path.read_text(encoding="utf-8")
    assert "Nearest NASA" not in html
    assert "site_07" in html


def test_index_links_every_site(tmp_path: Path) -> None:
    sites_df = pd.DataFrame(
        {
            "site_id": [f"site_{i:02d}" for i in range(1, 4)],
            "rank": [1, 2, 3],
            "score": [0.95, 0.92, 0.90],
            "lat": [-85.0, -86.0, -87.0],
            "lon": [-10.0, 0.0, 10.0],
        }
    )
    sites = gpd.GeoDataFrame(
        sites_df,
        geometry=[Point(row["lon"], row["lat"]) for _, row in sites_df.iterrows()],
        crs=LUNAR_GEOGRAPHIC_CRS,
    )
    out_dir = tmp_path / "sites"
    index_path = generate_site_index(sites, out_dir)
    assert index_path.exists()
    html = index_path.read_text(encoding="utf-8")
    for site_id in sites["site_id"]:
        assert f"{site_id}.html" in html
