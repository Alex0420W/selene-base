"""Tests for :mod:`selene_base.scoring.ranking`."""

from __future__ import annotations

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS

from selene_base.scoring.ranking import top_n_sites

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _score_grid(values: np.ndarray, pixel_m: float = 1000.0) -> xr.DataArray:
    h, w = values.shape
    half_x = (w / 2.0) * pixel_m
    half_y = (h / 2.0) * pixel_m
    da = xr.DataArray(
        values.astype(np.float64),
        dims=("y", "x"),
        coords={
            "y": np.linspace(half_y - pixel_m / 2, -half_y + pixel_m / 2, h),
            "x": np.linspace(-half_x + pixel_m / 2, half_x - pixel_m / 2, w),
        },
    )
    return da.rio.write_crs(LUNAR_SOUTH_POLAR, inplace=False)


class TestTopNSites:
    def test_picks_highest_first(self) -> None:
        arr = np.zeros((10, 10))
        arr[5, 5] = 0.95
        arr[2, 2] = 0.85
        arr[7, 7] = 0.75
        score = _score_grid(arr)
        sites = top_n_sites(score, n=3, min_distance_m=1.0, min_score=0.5)
        assert list(sites["score"]) == [0.95, 0.85, 0.75]
        assert list(sites["site_id"]) == ["site_01", "site_02", "site_03"]
        assert list(sites["rank"]) == [1, 2, 3]

    def test_min_distance_suppresses_neighbours(self) -> None:
        # Two peaks, 2 pixels (=2 km) apart. With min_distance 5 km, only one wins.
        arr = np.zeros((10, 10))
        arr[5, 4] = 0.9
        arr[5, 6] = 0.8
        score = _score_grid(arr)
        sites = top_n_sites(score, n=5, min_distance_m=5_000.0, min_score=0.5)
        assert len(sites) == 1
        assert sites["score"].iloc[0] == pytest.approx(0.9)

    def test_below_min_score_excluded(self) -> None:
        arr = np.zeros((6, 6))
        arr[1, 1] = 0.4  # below default min_score=0.5
        arr[4, 4] = 0.7
        score = _score_grid(arr)
        sites = top_n_sites(score, n=5, min_distance_m=1.0)
        assert len(sites) == 1
        assert sites["score"].iloc[0] == pytest.approx(0.7)

    def test_empty_when_nothing_meets_threshold(self) -> None:
        score = _score_grid(np.zeros((5, 5)))
        sites = top_n_sites(score, n=5, min_distance_m=1.0, min_score=0.5)
        assert len(sites) == 0
        # Schema must still include the per-criterion columns (v2.2 list).
        for crit in (
            "slope",
            "illumination",
            "thermal",
            "multi_volatile",
            "hazard",
            "pgv_seismic",
            "eva_psr_access",
            "los_to_earth",
        ):
            assert f"score_{crit}" in sites.columns

    def test_sub_scores_attached(self) -> None:
        arr = np.zeros((5, 5))
        arr[2, 2] = 0.9
        score = _score_grid(arr)
        slope = _score_grid(np.full((5, 5), 0.5))
        hazard = _score_grid(np.full((5, 5), 0.8))
        sites = top_n_sites(
            score,
            n=1,
            min_distance_m=1.0,
            min_score=0.5,
            sub_scores={"slope": slope, "hazard": hazard},
        )
        assert sites["score_slope"].iloc[0] == pytest.approx(0.5)
        assert sites["score_hazard"].iloc[0] == pytest.approx(0.8)
        # Untouched criteria stay NaN (v2.2 list).
        for crit in ("illumination", "thermal", "multi_volatile", "pgv_seismic"):
            assert np.isnan(sites[f"score_{crit}"].iloc[0])

    def test_lat_lon_present_and_finite(self) -> None:
        arr = np.zeros((5, 5))
        arr[2, 2] = 0.9  # centre pixel sits at polar (0, 0) → lat -90
        score = _score_grid(arr)
        sites = top_n_sites(score, n=1, min_distance_m=1.0, min_score=0.5)
        assert -90.0 <= sites["lat"].iloc[0] <= 90.0
        assert -360.0 <= sites["lon"].iloc[0] <= 360.0

    def test_n_must_be_positive(self) -> None:
        score = _score_grid(np.full((3, 3), 0.9))
        with pytest.raises(ValueError, match="n must be"):
            top_n_sites(score, n=0)

    def test_min_distance_must_be_positive(self) -> None:
        score = _score_grid(np.full((3, 3), 0.9))
        with pytest.raises(ValueError, match="min_distance_m"):
            top_n_sites(score, min_distance_m=0.0)

    def test_missing_crs_raises(self) -> None:
        bad = xr.DataArray(np.full((3, 3), 0.9), dims=("y", "x"))
        with pytest.raises(ValueError, match="no CRS"):
            top_n_sites(bad)
