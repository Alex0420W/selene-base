"""Tests for :mod:`selene_base.criteria.pgv_seismic` (v2.2+)."""

from __future__ import annotations

import math

import geopandas as gpd
import numpy as np
import pytest
import xarray as xr
from shapely.geometry import LineString, Point

from selene_base.criteria.pgv_seismic import (
    DEFAULT_ATTENUATION_L_KM,
    DEFAULT_CUTOFF_RADIUS_KM,
    DEFAULT_SIGMOID_MIDPOINT,
    DEFAULT_SIGMOID_SCALE,
    _densify_with_scarp_ids,
    compute,
    compute_components,
)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _grid(extent_m: float = 100_000, n: int = 5) -> xr.DataArray:
    """Square grid centred on origin in polar-stereographic metres."""
    arr = np.zeros((n, n), dtype=np.float32)
    half = extent_m / 2.0
    da = xr.DataArray(
        arr,
        dims=("y", "x"),
        coords={
            "y": np.linspace(half, -half, n),
            "x": np.linspace(-half, half, n),
        },
        name="dummy",
    )
    return da.rio.write_crs(POLAR_PROJ, inplace=False)


def _scarps_gdf(geoms: list) -> gpd.GeoDataFrame:
    if not geoms:
        return gpd.GeoDataFrame(geometry=[], crs=POLAR_PROJ)
    return gpd.GeoDataFrame(geometry=geoms, crs=POLAR_PROJ)


# ----- defaults -----


class TestDefaults:
    def test_attenuation_l_km_is_50(self) -> None:
        assert DEFAULT_ATTENUATION_L_KM == 50.0

    def test_cutoff_is_5l(self) -> None:
        assert DEFAULT_CUTOFF_RADIUS_KM == 5.0 * DEFAULT_ATTENUATION_L_KM

    def test_sigmoid_midpoint_anchored_at_strong_shaking(self) -> None:
        # exp(-40 km / 50 km) = exp(-0.8) ≈ 0.4493 — the Watters 2024
        # strong-shaking distance maps to score 0.5.
        anchor = math.exp(-40.0 / 50.0)
        assert abs(DEFAULT_SIGMOID_MIDPOINT - anchor) < 0.01


# ----- score range / basic behaviour -----


class TestScoreRange:
    def test_score_in_range_zero_to_one(self) -> None:
        # Random scarps; score must always land in [0, 1].
        rng = np.random.default_rng(0)
        coords = rng.uniform(-50_000, 50_000, size=(10, 4))
        scarps = _scarps_gdf([LineString([(c[0], c[1]), (c[2], c[3])]) for c in coords])
        out = compute(scarps, _grid())
        arr = out.to_numpy()
        assert np.all(arr >= 0.0)
        assert np.all(arr <= 1.0)

    def test_score_one_no_scarps_in_range(self) -> None:
        # Scarp far outside the cutoff radius. With grid extent ±50 km,
        # placing the scarp at x = 400 km means the closest cell (x =
        # +50 km) is 350 km from the scarp, comfortably beyond the
        # default 250 km cutoff.
        scarps = _scarps_gdf([LineString([(400_000, -500_000), (400_000, 500_000)])])
        out = compute_components(scarps, _grid(extent_m=100_000, n=5))
        # cum_pgv must be 0 for every cell; score sigmoid(0) = ~0.91.
        assert np.allclose(out["cum_pgv"].to_numpy(), 0.0)
        # Score saturates at the sigmoid value at cum_pgv=0:
        # 1 / (1 + exp((0 - 0.45)/0.20)) = 1 / (1 + exp(-2.25)) ≈ 0.905.
        np.testing.assert_allclose(
            out["score"].to_numpy(),
            1.0 / (1.0 + math.exp(-DEFAULT_SIGMOID_MIDPOINT / DEFAULT_SIGMOID_SCALE)),
        )

    def test_score_low_at_scarp_overlap(self) -> None:
        # Cell sitting on a scarp: cum_pgv approaches 1 (exp(-0/L) = 1),
        # score saturates near 0 (well below 0.5 midpoint).
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        out = compute(scarps, _grid())
        # Centre column (x=0) is on the scarp; corner (x=±50 km) is 50 km away.
        arr = out.to_numpy()
        assert arr[2, 2] < 0.1  # on scarp
        assert arr[2, 0] > arr[2, 2]  # far cell scores higher (safer)


# ----- attenuation kernel math -----


class TestAttenuationKernel:
    def test_decay_correct_at_known_distances(self) -> None:
        # Single straight scarp at x=0, dense enough for accurate
        # per-scarp-min distance (1km densification by default).
        scarps = _scarps_gdf([LineString([(0, -200_000), (0, 200_000)])])
        out = compute_components(scarps, _grid(extent_m=100_000, n=5))
        # Expected cum_pgv at x = ±25 km, ±50 km along the centre row.
        cum = out["cum_pgv"].to_numpy()
        # n=5 grid spans -50 km .. +50 km; cells at indices 0,1,2,3,4
        # correspond to x = -50, -25, 0, 25, 50 km (exact, since coords
        # are linspace(-half, half, n)).
        assert cum[2, 2] == pytest.approx(math.exp(-0.0 / 50.0), abs=1e-3)  # x=0
        assert cum[2, 1] == pytest.approx(math.exp(-25.0 / 50.0), abs=1e-3)  # x=-25 km
        assert cum[2, 3] == pytest.approx(math.exp(-25.0 / 50.0), abs=1e-3)  # x=+25 km
        assert cum[2, 0] == pytest.approx(math.exp(-50.0 / 50.0), abs=1e-3)  # x=-50 km
        assert cum[2, 4] == pytest.approx(math.exp(-50.0 / 50.0), abs=1e-3)  # x=+50 km

    def test_cumulative_summation_two_scarps(self) -> None:
        # Two parallel scarps at x = -25 km and x = +25 km. The cell at
        # x = 0 sees BOTH scarps at 25 km, so cum_pgv ≈ 2 * exp(-25/50).
        scarps = _scarps_gdf(
            [
                LineString([(-25_000, -200_000), (-25_000, 200_000)]),
                LineString([(25_000, -200_000), (25_000, 200_000)]),
            ]
        )
        out = compute_components(scarps, _grid(extent_m=100_000, n=5))
        cum = out["cum_pgv"].to_numpy()
        n_per = out["n_contributing_scarps"].to_numpy()
        expected = 2.0 * math.exp(-25_000 / 50_000)
        assert cum[2, 2] == pytest.approx(expected, abs=1e-3)
        # Both scarps inside cutoff at the centre cell.
        assert n_per[2, 2] == 2

    def test_cutoff_radius_excludes_distant_scarps(self) -> None:
        # Two scarps: one at x=0 (in range), one at x=300 km (outside
        # default 250 km cutoff). Cum_pgv at x=0 should reflect ONLY
        # the in-range scarp.
        scarps = _scarps_gdf(
            [
                LineString([(0, -200_000), (0, 200_000)]),
                LineString([(300_000, -200_000), (300_000, 200_000)]),
            ]
        )
        out = compute_components(scarps, _grid(extent_m=100_000, n=5))
        cum = out["cum_pgv"].to_numpy()
        n_per = out["n_contributing_scarps"].to_numpy()
        # cell at (2, 2) is at x=0, on the in-range scarp; only that one contributes.
        assert cum[2, 2] == pytest.approx(1.0, abs=1e-3)
        assert n_per[2, 2] == 1

    def test_kernel_parameters_respected(self) -> None:
        # Smaller L → faster decay, so cum_pgv at 25 km drops more.
        scarps = _scarps_gdf([LineString([(0, -200_000), (0, 200_000)])])
        out_l50 = compute_components(scarps, _grid(), attenuation_l_km=50.0)
        out_l25 = compute_components(scarps, _grid(), attenuation_l_km=25.0)
        # At x = ±25 km: L=50 gives exp(-0.5); L=25 gives exp(-1.0).
        assert out_l50["cum_pgv"].to_numpy()[2, 1] > out_l25["cum_pgv"].to_numpy()[2, 1]


# ----- sigmoid normalization -----


class TestSigmoid:
    def test_midpoint_at_strong_shaking(self) -> None:
        # cum_pgv == 0.45 → score 0.5 (default sigmoid).
        # Construct a synthetic single-scarp scene where the centre cell of
        # a 3x3 grid is exactly 40 km from a scarp at x = -40 km, so
        # exp(-40/50) ≈ 0.4493 ≈ midpoint.
        arr = np.zeros((3, 3), dtype=np.float32)
        da = xr.DataArray(
            arr,
            dims=("y", "x"),
            coords={"y": np.array([10_000, 0, -10_000]), "x": np.array([-10_000, 0, 10_000])},
            name="dummy",
        ).rio.write_crs(POLAR_PROJ, inplace=False)
        # Scarp at x = -40 km.
        s = _scarps_gdf([LineString([(-40_000, -200_000), (-40_000, 200_000)])])
        out = compute_components(s, da)
        # Centre cell distance ≈ 40 km.
        assert out["nearest_scarp_distance_km"].to_numpy()[1, 1] == pytest.approx(40.0, abs=0.05)
        # cum_pgv ≈ 0.4493; score ≈ 0.5.
        score = out["score"].to_numpy()[1, 1]
        assert score == pytest.approx(0.5, abs=0.02)

    def test_higher_cum_pgv_lowers_score(self) -> None:
        # Stack 5 parallel scarps at the same location → 5× cum_pgv at
        # nearby cells → much lower score than the single-scarp case.
        scarps_one = _scarps_gdf([LineString([(0, -200_000), (0, 200_000)])])
        scarps_five = _scarps_gdf(
            [
                LineString([(i, -200_000), (i, 200_000)])
                for i in [-200, -100, 0, 100, 200]  # 5 scarps within 1 km of x=0
            ]
        )
        out_one = compute(scarps_one, _grid()).to_numpy()
        out_five = compute(scarps_five, _grid()).to_numpy()
        # At x=0 column, the 5-scarp scene scores lower (riskier).
        assert out_five[2, 2] < out_one[2, 2]


# ----- edges + components -----


class TestEdgesAndComponents:
    def test_components_returns_four_keys(self) -> None:
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        out = compute_components(scarps, _grid())
        assert set(out.keys()) == {
            "score",
            "cum_pgv",
            "nearest_scarp_distance_km",
            "n_contributing_scarps",
        }

    def test_handles_empty_scarp_data(self) -> None:
        scarps = _scarps_gdf([])
        out = compute_components(scarps, _grid())
        # Empty catalog: cum_pgv = 0 everywhere, n_contributing = 0,
        # nearest = inf, score saturates at sigmoid(0).
        assert np.allclose(out["cum_pgv"].to_numpy(), 0.0)
        assert np.all(out["n_contributing_scarps"].to_numpy() == 0)
        assert np.all(np.isinf(out["nearest_scarp_distance_km"].to_numpy()))
        # Score is the sigmoid(0) saturation value.
        expected = 1.0 / (1.0 + math.exp(-DEFAULT_SIGMOID_MIDPOINT / DEFAULT_SIGMOID_SCALE))
        np.testing.assert_allclose(out["score"].to_numpy(), expected)

    def test_densified_points_at_correct_spacing(self) -> None:
        # 100 km straight line, 5 km densify spacing → ~21 points.
        scarps = _scarps_gdf([LineString([(0, 0), (100_000, 0)])])
        pts, ids = _densify_with_scarp_ids(scarps, POLAR_PROJ, spacing_m=5_000)
        assert pts.shape[1] == 2
        # Point count: ceil(100/5)+1 = 21
        assert pts.shape[0] >= 20
        # All points belong to the same scarp (id 0).
        assert np.all(ids == 0)

    def test_score_smoothness_versus_distance(self) -> None:
        # Single scarp, sweep distance: score should monotonically
        # increase with distance (further = safer).
        scarps = _scarps_gdf([LineString([(0, -200_000), (0, 200_000)])])
        # Build a row of cells at x = 0, 10, 20, 30, 40, 50, 60 km.
        xs_km = np.array([0, 10, 20, 30, 40, 50, 60])
        arr = np.zeros((1, len(xs_km)), dtype=np.float32)
        da = xr.DataArray(
            arr,
            dims=("y", "x"),
            coords={"y": np.array([0]), "x": xs_km * 1000.0},
            name="dummy",
        ).rio.write_crs(POLAR_PROJ, inplace=False)
        out = compute(scarps, da).to_numpy()[0]
        # Monotonic increase: score(0 km) ≤ score(10 km) ≤ ...
        diffs = np.diff(out)
        assert np.all(diffs >= -1e-6)


# ----- parameter validation -----


class TestParameterValidation:
    def test_invalid_attenuation_l_raises(self) -> None:
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        with pytest.raises(ValueError, match="attenuation_l_km"):
            compute(scarps, _grid(), attenuation_l_km=0.0)

    def test_invalid_cutoff_raises(self) -> None:
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        with pytest.raises(ValueError, match="cutoff_radius_km"):
            compute(scarps, _grid(), cutoff_radius_km=-10.0)

    def test_invalid_densify_spacing_raises(self) -> None:
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        with pytest.raises(ValueError, match="densify_spacing_m"):
            compute(scarps, _grid(), densify_spacing_m=0.0)

    def test_invalid_sigmoid_scale_raises(self) -> None:
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        with pytest.raises(ValueError, match="sigmoid_scale"):
            compute(scarps, _grid(), sigmoid_scale=0.0)

    def test_missing_target_grid_crs_raises(self) -> None:
        bad = xr.DataArray(
            np.zeros((3, 3), dtype=np.float32),
            dims=("y", "x"),
            coords={"y": np.array([0, 1, 2]), "x": np.array([0, 1, 2])},
        )
        scarps = _scarps_gdf([LineString([(0, -100_000), (0, 100_000)])])
        with pytest.raises(ValueError, match="no CRS"):
            compute(scarps, bad)

    def test_missing_scarp_crs_raises(self) -> None:
        scarps = gpd.GeoDataFrame(
            geometry=[LineString([(0, -100_000), (0, 100_000)])],
            crs=None,
        )
        with pytest.raises(ValueError, match="no CRS"):
            compute(scarps, _grid())


# ----- point-geometry handling -----


class TestPointGeometries:
    def test_point_scarps_treated_as_centroids(self) -> None:
        # The Watters 2015 polar bundle is Point geometries; the criterion
        # should still produce a sensible score (each point treated as a
        # zero-length scarp).
        scarps = _scarps_gdf([Point(0, 0), Point(100_000, 0)])
        out = compute_components(scarps, _grid())
        # cum_pgv at the centre (0, 0) is at least exp(0) = 1 from the
        # point at origin.
        assert out["cum_pgv"].to_numpy()[2, 2] >= 1.0
        # Both points are inside the 250 km cutoff → n_contributing = 2.
        assert out["n_contributing_scarps"].to_numpy()[2, 2] == 2
