"""Tests for :mod:`selene_base.criteria.seismic` (v1.8 — eighth criterion live).

Two suites:

- :class:`TestDistanceToScarps` — KDTree-based distance from each grid
  cell to the nearest densified scarp vertex. Synthetic scarps in the
  lunar south polar CRS.
- :class:`TestCompute` — logistic mapping from distance (km) to a
  ``[0, 1]`` safety score, with the v1.8 defaults centred on Civilini
  et al. (2023)'s shallow-moonquake-to-scarp clustering distance.

Plus a smoke test that loads the bundled Mishra & Kumar 2022 catalog
end-to-end through ``distance_to_scarps`` + ``compute`` on a small
polar-stereo grid (the "criterion works on the polar grid" test from
the v1.8 spec).
"""

from __future__ import annotations

import math

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS
from shapely.geometry import LineString, Point

from selene_base.criteria.seismic import (
    BUNDLED_MISHRA_KUMAR_2022,
    compute,
    distance_to_scarps,
    load_bundled_catalog,
)

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _polar_grid(width: int, height: int, pixel_m: float) -> xr.DataArray:
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


# --- Bundled catalog ---------------------------------------------------


def test_load_bundled_catalog_returns_polylines() -> None:
    """The bundled Mishra & Kumar 2022 file loads as polylines, not points."""
    cat = load_bundled_catalog()
    assert len(cat) == 704, "v1.8 ships exactly 704 main-segment features"
    assert cat.crs is not None, "bundled catalog must carry an explicit CRS"
    geom_types = set(cat.geometry.type.unique())
    assert geom_types <= {"LineString", "MultiLineString"}, (
        f"expected line geometries, got {sorted(geom_types)}"
    )


def test_bundled_catalog_path_constant_resolves() -> None:
    assert BUNDLED_MISHRA_KUMAR_2022.exists(), (
        f"bundled shapefile path not found: {BUNDLED_MISHRA_KUMAR_2022}"
    )


# --- Distance -----------------------------------------------------------


class TestDistanceToScarps:
    def test_empty_catalog_returns_inf(self) -> None:
        grid = _polar_grid(5, 5, pixel_m=1000.0)
        empty = gpd.GeoDataFrame({"id": []}, geometry=[], crs=LUNAR_SOUTH_POLAR)
        out = distance_to_scarps(empty, grid).to_numpy()
        assert np.all(np.isinf(out))

    def test_point_at_origin(self) -> None:
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        scarps = gpd.GeoDataFrame({"id": [0]}, geometry=[Point(0.0, 0.0)], crs=LUNAR_SOUTH_POLAR)
        out = distance_to_scarps(scarps, grid).to_numpy()
        # Centre pixel (5, 5) is at (0, 0) — distance ~ 0.
        assert out[5, 5] == pytest.approx(0.0, abs=0.05)
        # Diagonal corner at (~5 km, 5 km) → ~7.07 km
        diag = out[0, 0]
        assert 6.5 < diag < 8.0

    def test_line_distance(self) -> None:
        # Vertical scarp along x=0; pixel (row=5, col=8) is 3 pixels east.
        grid = _polar_grid(11, 11, pixel_m=1000.0)
        line = LineString([(0.0, -10_000.0), (0.0, 10_000.0)])
        scarps = gpd.GeoDataFrame({"id": [0]}, geometry=[line], crs=LUNAR_SOUTH_POLAR)
        out = distance_to_scarps(scarps, grid).to_numpy()
        assert out[5, 5] == pytest.approx(0.0, abs=0.05)
        assert out[5, 8] == pytest.approx(3.0, abs=0.1)
        assert out[5, 2] == pytest.approx(3.0, abs=0.1)


# --- Scoring ------------------------------------------------------------


class TestCompute:
    """Logistic mapping ``score(d) = 1 / (1 + exp(-(d - midpoint)/steepness))``.

    With v1.8 defaults (midpoint 25 km, steepness 8 km):
      d=0  → 0.04   d=5  → 0.08   d=15 → 0.22   d=25 → 0.50
      d=35 → 0.78   d=50 → 0.96   d→∞  → 1.00
    """

    @pytest.mark.parametrize(
        "d, expected",
        [
            (0.0, 1.0 / (1.0 + math.exp(25.0 / 8.0))),
            (5.0, 1.0 / (1.0 + math.exp(20.0 / 8.0))),
            (25.0, 0.5),
            (50.0, 1.0 / (1.0 + math.exp(-25.0 / 8.0))),
            (1000.0, 1.0),  # asymptote
        ],
    )
    def test_logistic_at_known_distances(self, d: float, expected: float) -> None:
        out = compute(xr.DataArray(np.array([[d]]), dims=("y", "x"))).to_numpy()
        assert out[0, 0] == pytest.approx(expected, abs=1e-3)

    def test_close_cell_low_score(self) -> None:
        """A cell within Civilini's ~5 km moonquake-clustering radius must
        score in the unsafe band."""
        d = xr.DataArray(np.array([[3.0]]), dims=("y", "x"))
        out = compute(d).to_numpy()
        assert out[0, 0] < 0.10

    def test_far_cell_high_score(self) -> None:
        """A cell well beyond 50 km must score effectively safe."""
        d = xr.DataArray(np.array([[80.0]]), dims=("y", "x"))
        out = compute(d).to_numpy()
        assert out[0, 0] > 0.99

    def test_scoring_function_monotonic(self) -> None:
        """Closer distance → strictly lower score, across a sweep."""
        d = xr.DataArray(
            np.array([[0.0, 5.0, 10.0, 25.0, 50.0, 100.0]]),
            dims=("y", "x"),
        )
        out = compute(d).to_numpy().ravel()
        assert np.all(np.diff(out) > 0), f"expected strictly increasing, got {out}"

    def test_scoring_function_in_range_0_1(self) -> None:
        """Across a wide distance sweep, scores stay in [0, 1] strictly."""
        rng = np.random.default_rng(20260501)
        d = xr.DataArray(rng.uniform(0.0, 200.0, size=(10, 10)), dims=("y", "x"))
        out = compute(d).to_numpy()
        finite = out[np.isfinite(out)]
        assert finite.min() >= 0.0
        assert finite.max() <= 1.0

    def test_inf_distance_treated_as_safe(self) -> None:
        d = xr.DataArray(np.array([[np.inf]]), dims=("y", "x"))
        out = compute(d).to_numpy()
        assert out[0, 0] == pytest.approx(1.0)

    def test_nan_propagates(self) -> None:
        d = xr.DataArray(np.array([[np.nan, 25.0]]), dims=("y", "x"))
        out = compute(d).to_numpy()
        assert math.isnan(out[0, 0])
        assert out[0, 1] == pytest.approx(0.5, abs=1e-9)

    def test_midpoint_must_be_positive(self) -> None:
        d = xr.DataArray(np.array([[10.0]]), dims=("y", "x"))
        with pytest.raises(ValueError, match="midpoint_km"):
            compute(d, midpoint_km=0.0)

    def test_steepness_must_be_positive(self) -> None:
        d = xr.DataArray(np.array([[10.0]]), dims=("y", "x"))
        with pytest.raises(ValueError, match="steepness_km"):
            compute(d, steepness_km=-1.0)


# --- End-to-end smoke test ---------------------------------------------


def test_seismic_criterion_works_for_polar_grid() -> None:
    """Bundled Mishra & Kumar catalog → distance grid → score grid on a
    coarse south-polar test grid. Verifies the full criterion chain runs
    without CRS or memory issues against real data."""
    catalog = load_bundled_catalog()
    # 25-km pixels over a ±300 km square — small enough for a unit-test
    # KDTree query, large enough that some catalog segments fall inside.
    grid = _polar_grid(25, 25, pixel_m=25_000.0)
    distance = distance_to_scarps(catalog, grid)
    score = compute(distance)

    arr = score.to_numpy()
    assert arr.shape == (25, 25)
    finite = arr[np.isfinite(arr)]
    assert finite.size == arr.size, "no NaN should appear with a non-empty catalog"
    assert finite.min() >= 0.0
    assert finite.max() <= 1.0
    # The 704-segment catalog covers the polar region densely; at least
    # some pixels in a ±300 km square should be in the unsafe band.
    assert finite.min() < 0.5, "some grid cells should sit inside the unsafe band"
    # And at least some should be safe (the catalog isn't dense everywhere).
    assert finite.max() > 0.5, "some grid cells should sit in the safe band"
