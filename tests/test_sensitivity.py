"""Tests for :mod:`selene_base.validation.sensitivity`."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.sensitivity import (
    best_weights,
    latin_hypercube_weights,
    sweep_weights,
)

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


class TestLatinHypercubeWeights:
    def test_rows_sum_to_one(self) -> None:
        weights = latin_hypercube_weights(50, ["a", "b", "c", "d"])
        assert weights.shape == (50, 4)
        np.testing.assert_allclose(weights.sum(axis=1), 1.0, atol=1e-9)

    def test_all_non_negative(self) -> None:
        weights = latin_hypercube_weights(50, ["a", "b", "c"])
        assert (weights >= 0).all()
        assert (weights <= 1).all()

    def test_one_dim_returns_ones(self) -> None:
        weights = latin_hypercube_weights(8, ["only"])
        assert weights.shape == (8, 1)
        np.testing.assert_allclose(weights, 1.0)

    def test_seed_reproducibility(self) -> None:
        first = latin_hypercube_weights(10, ["a", "b", "c"], seed=7)
        second = latin_hypercube_weights(10, ["a", "b", "c"], seed=7)
        np.testing.assert_allclose(first, second)
        third = latin_hypercube_weights(10, ["a", "b", "c"], seed=11)
        assert not np.allclose(first, third)

    def test_n_samples_validated(self) -> None:
        with pytest.raises(ValueError, match="n_samples"):
            latin_hypercube_weights(0, ["a", "b"])

    def test_empty_criterion_names_rejected(self) -> None:
        with pytest.raises(ValueError, match="criterion_names"):
            latin_hypercube_weights(5, [])


def _polar_score_grid(values: np.ndarray, pixel_m: float = 5000.0) -> xr.DataArray:
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


class TestSweepWeights:
    def test_shape_and_columns(self) -> None:
        # Two criteria, A and B. A peaks in the top-left quadrant, B in
        # bottom-right; weight regimes will accordingly favour different
        # corners of the synthetic grid.
        h = w = 41
        ramp = np.linspace(0.0, 1.0, w)
        a = np.broadcast_to(ramp, (h, w)).copy()
        b = np.broadcast_to(ramp[::-1], (h, w)).copy()
        score_maps = {"a": _polar_score_grid(a), "b": _polar_score_grid(b)}

        # Single fake "NASA region" at the center, where score balance
        # has the most leverage on the per-region distance number.
        nasa = gpd.GeoDataFrame(
            {"name": ["fake"], "lat": [-90.0], "lon": [0.0], "radius_km": [10.0]},
            geometry=[Point(0.0, -90.0)],
            crs=LUNAR_GEOGRAPHIC_CRS,
        )

        samples = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]])
        results = sweep_weights(
            score_maps,
            samples,
            nasa,
            top_n=5,
            min_distance_km=1.0,
            proximity_threshold_km=10.0,
            criterion_order=["a", "b"],
        )
        assert len(results) == 3
        for col in (
            "w_a",
            "w_b",
            "n_inside_region",
            "n_within_proximity_km",
            "n_regions_with_top_site",
            "n_regions_within_proximity_km",
            "mean_score_at_nasa_centroids",
            "mean_top_n_score",
        ):
            assert col in results.columns
        # Top-N mean score under the all-A weighting equals the all-A
        # max score (1.0); under all-B also 1.0; under 50/50 is the
        # equal-mix peak (0.5 + 0.5 = 1.0). That should cap at 1.
        assert results["mean_top_n_score"].max() <= 1.0 + 1e-9

    def test_shape_mismatch_rejected(self) -> None:
        score_maps = {
            "a": _polar_score_grid(np.zeros((5, 5))),
            "b": _polar_score_grid(np.zeros((5, 5))),
        }
        nasa = gpd.GeoDataFrame(
            {"name": ["fake"], "lat": [-89.0], "lon": [0.0], "radius_km": [10.0]},
            geometry=[Point(0.0, -89.0)],
            crs=LUNAR_GEOGRAPHIC_CRS,
        )
        bad = np.array([[0.5, 0.5, 0.0]])  # 3 columns for 2 criteria
        with pytest.raises(ValueError, match="weight_samples"):
            sweep_weights(
                score_maps,
                bad,
                nasa,
                criterion_order=["a", "b"],
            )

    def test_empty_score_maps_rejected(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            sweep_weights({}, np.zeros((1, 1)), gpd.GeoDataFrame())


def test_best_weights_picks_max_match_count() -> None:
    import pandas as pd

    df = pd.DataFrame(
        {
            "w_a": [0.1, 0.5, 0.9],
            "w_b": [0.9, 0.5, 0.1],
            "n_regions_within_proximity_km": [3, 7, 5],
            "n_regions_with_top_site": [1, 2, 3],
            "n_inside_region": [0, 1, 0],
            "mean_score_at_nasa_centroids": [0.4, 0.6, 0.5],
        }
    )
    best = best_weights(df, ["a", "b"])
    # Sample with 7 region matches wins — that's the second row.
    assert best == {"a": 0.5, "b": 0.5}


def test_best_weights_empty_raises() -> None:
    import pandas as pd

    with pytest.raises(ValueError, match="empty"):
        best_weights(pd.DataFrame(), ["a"])
