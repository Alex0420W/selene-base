"""Tests for :mod:`selene_base.validation.wueller_scoring` (v1.5).

Verifies that ``score_wueller_sites`` samples the right pixel and
returns the correct columns. Uses a tiny in-memory raster stack
written to a temporary directory so the test does not depend on the
240 m COGs being on disk.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr
from pyproj import Transformer
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.wueller_comparison import POLAR_PROJ
from selene_base.validation.wueller_scoring import (
    AGGREGATE_FILENAME,
    CRITERION_RASTER_FILENAMES,
    DEFAULT_SCORED_SUBDIR,
    RAW_RASTER_FILENAMES,
    score_wueller_sites,
)


def _make_raster(value_grid: np.ndarray, out_path: Path) -> None:
    """Write a tiny 2D raster in polar-stereo CRS to ``out_path``."""
    height, width = value_grid.shape
    pixel_size = 240.0
    xs = np.linspace(-pixel_size * (width - 1) / 2, pixel_size * (width - 1) / 2, width)
    ys = np.linspace(pixel_size * (height - 1) / 2, -pixel_size * (height - 1) / 2, height)
    da = xr.DataArray(
        value_grid.astype(np.float32),
        dims=("y", "x"),
        coords={"y": ys, "x": xs},
        name=out_path.stem,
    ).rio.write_crs(POLAR_PROJ, inplace=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    da.rio.to_raster(out_path)


def _wueller_two_sites() -> gpd.GeoDataFrame:
    """Two Wueller-style sites near the south pole, in lunar geographic."""
    transformer = Transformer.from_crs(POLAR_PROJ, LUNAR_GEOGRAPHIC_CRS, always_xy=True)
    lon0, lat0 = transformer.transform(0.0, 0.0)  # pole-adjacent pixel center
    lon1, lat1 = transformer.transform(240.0, 0.0)  # one pixel east
    rows = [
        {
            "wueller_site_id": "W01",
            "region": "Mons Mouton",
            "lat": float(lat0),
            "lon": float(lon0),
            "in_usgs_scope": True,
        },
        {
            "wueller_site_id": "W02",
            "region": "Haworth",
            "lat": float(lat1),
            "lon": float(lon1),
            "in_usgs_scope": True,
        },
    ]
    geoms = [Point(r["lon"], r["lat"]) for r in rows]
    return gpd.GeoDataFrame(rows, geometry=geoms, crs=LUNAR_GEOGRAPHIC_CRS)


def test_score_wueller_sites_samples_correct_pixel(tmp_path: Path) -> None:
    """Two Wueller sites land on two distinct pixels and pick up the right values."""
    processed_dir = tmp_path / "processed"
    outputs_dir = tmp_path / "outputs"
    scored_dir = processed_dir / DEFAULT_SCORED_SUBDIR

    # 3x3 grid; pole is at the centre pixel (1,1). Site W01 sits at (0, 0)
    # in polar-stereo metres → centre cell. Site W02 sits at (240, 0) → one
    # cell east, i.e. column 2 (the +x side).
    centre = 1
    east = 2
    grids: dict[str, np.ndarray] = {}
    for crit in CRITERION_RASTER_FILENAMES:
        g = np.zeros((3, 3), dtype=np.float32)
        g[centre, centre] = 0.5  # centre pixel score
        g[centre, east] = 0.9  # eastern pixel score
        grids[crit] = g
        _make_raster(g, scored_dir / CRITERION_RASTER_FILENAMES[crit])
    # Aggregate
    agg = np.zeros((3, 3), dtype=np.float32)
    agg[centre, centre] = 0.6
    agg[centre, east] = 0.7
    _make_raster(agg, outputs_dir / AGGREGATE_FILENAME)
    # Raw inputs — slope 1° (compliant) at centre, 20° (steep) at east;
    # illumination 0.5 (compliant) at both; los_visibility 0.6 (compliant)
    # at centre, 0.1 (fail) at east.
    slope = np.zeros((3, 3), dtype=np.float32)
    slope[centre, centre] = 1.0
    slope[centre, east] = 20.0
    _make_raster(slope, processed_dir / RAW_RASTER_FILENAMES["slope_deg"])
    illum = np.full((3, 3), 0.5, dtype=np.float32)
    _make_raster(illum, processed_dir / RAW_RASTER_FILENAMES["illumination"])
    los = np.full((3, 3), 0.1, dtype=np.float32)
    los[centre, centre] = 0.6
    _make_raster(los, processed_dir / RAW_RASTER_FILENAMES["los_visibility"])

    sites = _wueller_two_sites()
    df = score_wueller_sites(
        wueller_sites=sites,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
    )

    assert list(df["wueller_site_id"]) == ["W01", "W02"]
    # W01 (centre): score_slope = 0.5; W02 (east): score_slope = 0.9
    assert df.loc[0, "score_slope"] == pytest.approx(0.5, rel=1e-6)
    assert df.loc[1, "score_slope"] == pytest.approx(0.9, rel=1e-6)
    assert df.loc[0, "aggregate_score"] == pytest.approx(0.6, rel=1e-6)
    assert df.loc[1, "aggregate_score"] == pytest.approx(0.7, rel=1e-6)
    # Raw values
    assert df.loc[0, "slope_deg"] == pytest.approx(1.0, rel=1e-6)
    assert df.loc[1, "slope_deg"] == pytest.approx(20.0, rel=1e-6)
    # HLS: centre passes (slope ok, illum ok, los ok); east fails (slope > 8°,
    # los < 0.5).
    assert bool(df.loc[0, "hls_compliant"]) is True
    assert bool(df.loc[1, "hls_compliant"]) is False


def test_score_wueller_sites_handles_missing_rasters(tmp_path: Path) -> None:
    """Missing per-criterion rasters yield NaN, not a crash."""
    processed_dir = tmp_path / "processed"
    outputs_dir = tmp_path / "outputs"
    scored_dir = processed_dir / DEFAULT_SCORED_SUBDIR
    scored_dir.mkdir(parents=True)
    outputs_dir.mkdir(parents=True)

    sites = _wueller_two_sites()
    df = score_wueller_sites(
        wueller_sites=sites,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
    )
    for crit in CRITERION_RASTER_FILENAMES:
        assert df[f"score_{crit}"].isna().all(), f"expected NaN column for {crit}"
    assert df["aggregate_score"].isna().all()
    # HLS evaluates to False when all raw inputs are NaN (np.isfinite is False).
    assert (df["hls_compliant"] == False).all()  # noqa: E712 — explicit numpy compare
