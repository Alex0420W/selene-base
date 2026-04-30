"""Tests for :func:`selene_base.scoring.ranking.top_n_sites_per_region` (week 11).

The function searches *within* each NASA region polygon, applies the
NASA HLS hard-constraint filters (slope, slope-buffer, illumination,
DTE visibility), then NMS-ranks the survivors by score. These tests
exercise:

- the all-compliant happy path (a polygon with cells passing every
  filter returns top-N),
- the no-compliant path (a polygon where every cell fails at least one
  filter — the function records nothing for that region but doesn't
  error),
- each individual hard filter (slope, buffer, illumination, DTE),
- the per-region NMS (close pairs are suppressed at min_distance_m),
- shape mismatches and parameter validation.

Synthetic grids are used throughout — no real LRO data required.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from shapely.geometry import Polygon

from selene_base.scoring.ranking import (
    HLS_BUFFER_M,
    HLS_DTE_VISIBILITY_MIN,
    HLS_ILLUMINATION_MIN,
    HLS_SLOPE_MAX_DEG,
    top_n_sites_per_region,
)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _make_grid(
    arr: np.ndarray,
    *,
    pixel_size_m: float = 240.0,
    name: str = "data",
) -> xr.DataArray:
    """Build a 2D xarray-with-CRS test grid centred at (0, 0) projected metres."""
    height, width = arr.shape
    half_w = width * pixel_size_m / 2.0
    half_h = height * pixel_size_m / 2.0
    x_coords = np.linspace(-half_w + pixel_size_m / 2.0, half_w - pixel_size_m / 2.0, width)
    # Standard COG y orientation: row 0 sits at the top (largest y).
    y_coords = np.linspace(half_h - pixel_size_m / 2.0, -half_h + pixel_size_m / 2.0, height)
    da = xr.DataArray(
        arr.astype(np.float64),
        dims=("y", "x"),
        coords={"y": y_coords, "x": x_coords},
        name=name,
    ).rio.write_crs(POLAR_PROJ, inplace=False)
    return da


def _square_polygon(
    half_size_m: float, *, region: str = "Test", code: str = "TT"
) -> gpd.GeoDataFrame:
    """A single square polygon centred on the origin, in POLAR_PROJ metres."""
    poly = Polygon(
        [
            (-half_size_m, -half_size_m),
            (half_size_m, -half_size_m),
            (half_size_m, half_size_m),
            (-half_size_m, half_size_m),
            (-half_size_m, -half_size_m),
        ]
    )
    return gpd.GeoDataFrame(
        {"Region": [region], "RegionCode": [code], "Area_km2": [400.0]},
        geometry=[poly],
        crs=POLAR_PROJ,
    )


class TestHappyPath:
    def test_compliant_polygon_returns_top_n(self) -> None:
        # 21x21 grid of pixels (240m each → 5040 m extent), polygon
        # covers the centre 9x9 region. Score increases toward centre,
        # all four filters trivially passed everywhere.
        grid_size = 21
        score = np.zeros((grid_size, grid_size), dtype=np.float64)
        cy, cx = grid_size // 2, grid_size // 2
        for r in range(grid_size):
            for c in range(grid_size):
                d = ((r - cy) ** 2 + (c - cx) ** 2) ** 0.5
                score[r, c] = 1.0 / (1.0 + d)

        slope = np.full_like(score, 2.0)  # 2° everywhere — passes 8° filter
        illum = np.full_like(score, 0.7)
        los = np.full_like(score, 0.8)
        score_da = _make_grid(score, name="score")
        slope_da = _make_grid(slope, name="slope_deg")
        illum_da = _make_grid(illum, name="illumination")
        los_da = _make_grid(los, name="los_visibility")
        polygon = _square_polygon(half_size_m=1080.0)  # ±1080 m → 9x9 cells

        result = top_n_sites_per_region(
            score_da,
            polygon,
            slope_deg=slope_da,
            illumination=illum_da,
            los_visibility=los_da,
            n_per_region=3,
            min_distance_m=240.0,
        )
        assert len(result) == 3
        assert (result["region_name"] == "Test").all()
        assert (result["region_code"] == "TT").all()
        assert list(result["rank_in_region"]) == [1, 2, 3]
        # Best-ranked site should be the highest-scoring one.
        assert result["score"].iloc[0] >= result["score"].iloc[1]
        assert result["score"].iloc[1] >= result["score"].iloc[2]
        assert result["hls_compliant"].all()


class TestHLSFilters:
    """Each of the four HLS filters disqualifies the right cells."""

    def _baseline_compliant_grid(
        self,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
        score = np.full((21, 21), 0.8, dtype=np.float64)
        slope = np.full((21, 21), 2.0, dtype=np.float64)
        illum = np.full((21, 21), 0.7, dtype=np.float64)
        los = np.full((21, 21), 0.8, dtype=np.float64)
        return (
            _make_grid(score, name="score"),
            _make_grid(slope, name="slope_deg"),
            _make_grid(illum, name="illumination"),
            _make_grid(los, name="los_visibility"),
        )

    def test_slope_too_steep_disqualifies(self) -> None:
        score, slope, illum, los = self._baseline_compliant_grid()
        slope = slope * 0.0 + (HLS_SLOPE_MAX_DEG + 5.0)  # all steeper than 8°
        polygon = _square_polygon(half_size_m=1080.0)
        result = top_n_sites_per_region(
            score,
            polygon,
            slope_deg=slope,
            illumination=illum,
            los_visibility=los,
        )
        assert len(result) == 0

    def test_illumination_too_low_disqualifies(self) -> None:
        score, slope, illum, los = self._baseline_compliant_grid()
        illum = illum * 0.0 + (HLS_ILLUMINATION_MIN - 0.05)
        polygon = _square_polygon(half_size_m=1080.0)
        result = top_n_sites_per_region(
            score,
            polygon,
            slope_deg=slope,
            illumination=illum,
            los_visibility=los,
        )
        assert len(result) == 0

    def test_dte_visibility_too_low_disqualifies(self) -> None:
        score, slope, illum, los = self._baseline_compliant_grid()
        los = los * 0.0 + (HLS_DTE_VISIBILITY_MIN - 0.05)
        polygon = _square_polygon(half_size_m=1080.0)
        result = top_n_sites_per_region(
            score,
            polygon,
            slope_deg=slope,
            illumination=illum,
            los_visibility=los,
        )
        assert len(result) == 0

    def test_buffer_filter_excludes_near_steep_cells(self) -> None:
        # Build a 41×41 grid where all cells are buildable (slope=2°)
        # except for a strip along col=20 which is 15° (steep). That
        # strip should propagate a 100 m exclusion zone around it; cells
        # within 100 m of the strip should be disqualified by the
        # buffer filter even though they themselves are flat.
        size = 41
        score = np.full((size, size), 0.8, dtype=np.float64)
        slope = np.full((size, size), 2.0, dtype=np.float64)
        slope[:, 20] = 15.0
        illum = np.full((size, size), 0.7, dtype=np.float64)
        los = np.full((size, size), 0.8, dtype=np.float64)
        score_da = _make_grid(score, name="score")
        slope_da = _make_grid(slope, name="slope_deg")
        illum_da = _make_grid(illum, name="illumination")
        los_da = _make_grid(los, name="los_visibility")
        polygon = _square_polygon(half_size_m=4920.0)  # 41 cells × 240m

        # 100m buffer at 240 m/cell ≈ exclude any cell within ~0.42 of
        # the steep cell, which means immediate neighbours (1 pixel
        # away = 240 m, > 100 m) are NOT excluded. Switch to a 300m
        # buffer to make the test bite: any cell within 300 m of the
        # steep strip should be disqualified.
        result = top_n_sites_per_region(
            score_da,
            polygon,
            slope_deg=slope_da,
            illumination=illum_da,
            los_visibility=los_da,
            hls_buffer_m=300.0,
            n_per_region=10,
            min_distance_m=240.0,
        )
        # No accepted site should be on or directly adjacent to the
        # steep strip at col 20: the strip's pixel-x extent is x ≈ 0,
        # so any accepted (x_m) within ±300 m of 0 should be empty.
        for x_m in result["x_m"]:
            assert abs(x_m) > 240.0, f"site at x={x_m} was within buffer zone"

    def test_buffer_filter_at_default_threshold_keeps_close_neighbours(self) -> None:
        # At the 100m default buffer, neighbours one pixel (240m) away
        # from a steep strip should pass: their distance to the steep
        # cell is 240m > 100m. Same setup as above.
        size = 41
        score = np.full((size, size), 0.8, dtype=np.float64)
        slope = np.full((size, size), 2.0, dtype=np.float64)
        slope[:, 20] = 15.0
        illum = np.full((size, size), 0.7, dtype=np.float64)
        los = np.full((size, size), 0.8, dtype=np.float64)
        score_da = _make_grid(score, name="score")
        slope_da = _make_grid(slope, name="slope_deg")
        illum_da = _make_grid(illum, name="illumination")
        los_da = _make_grid(los, name="los_visibility")
        polygon = _square_polygon(half_size_m=4920.0)

        result = top_n_sites_per_region(
            score_da,
            polygon,
            slope_deg=slope_da,
            illumination=illum_da,
            los_visibility=los_da,
            hls_buffer_m=HLS_BUFFER_M,
            n_per_region=20,
            min_distance_m=240.0,
        )
        # We should get up to 20 sites; the steep strip itself is
        # excluded but everything outside it should be eligible.
        assert len(result) > 0


class TestNMSWithinRegion:
    def test_min_distance_enforced_within_region(self) -> None:
        # Two adjacent peaks separated by 1 pixel (240 m); set
        # min_distance_m = 480 m so only one survives.
        size = 11
        score = np.full((size, size), 0.5, dtype=np.float64)
        score[5, 4] = 1.0
        score[5, 6] = 0.99  # just slightly lower so order is deterministic
        slope = np.full((size, size), 2.0, dtype=np.float64)
        illum = np.full((size, size), 0.7, dtype=np.float64)
        los = np.full((size, size), 0.8, dtype=np.float64)
        score_da = _make_grid(score, name="score")
        slope_da = _make_grid(slope, name="slope_deg")
        illum_da = _make_grid(illum, name="illumination")
        los_da = _make_grid(los, name="los_visibility")
        polygon = _square_polygon(half_size_m=1320.0)

        result = top_n_sites_per_region(
            score_da,
            polygon,
            slope_deg=slope_da,
            illumination=illum_da,
            los_visibility=los_da,
            n_per_region=2,
            min_distance_m=600.0,
        )
        # The two near-twin peaks were 480 m apart (2 pixels at 240 m
        # each). min_distance_m = 600 m strictly exceeds 480 m, so the
        # second peak is suppressed and the next accepted site falls on
        # the 0.5 plateau elsewhere in the polygon.
        assert len(result) == 2
        assert result["score"].iloc[0] == pytest.approx(1.0)
        assert result["score"].iloc[1] < 0.95


class TestMultiRegionBookkeeping:
    def test_site_ids_globally_unique_across_regions(self) -> None:
        size = 11
        score = np.full((size, size), 0.8, dtype=np.float64)
        slope = np.full((size, size), 2.0, dtype=np.float64)
        illum = np.full((size, size), 0.7, dtype=np.float64)
        los = np.full((size, size), 0.8, dtype=np.float64)
        score_da = _make_grid(score, name="score")
        slope_da = _make_grid(slope, name="slope_deg")
        illum_da = _make_grid(illum, name="illumination")
        los_da = _make_grid(los, name="los_visibility")

        # Two non-overlapping square polygons within the grid.
        poly_a = Polygon(
            [(-1000, -1000), (-200, -1000), (-200, -200), (-1000, -200), (-1000, -1000)]
        )
        poly_b = Polygon([(200, 200), (1000, 200), (1000, 1000), (200, 1000), (200, 200)])
        polygons = gpd.GeoDataFrame(
            {"Region": ["A", "B"], "RegionCode": ["AA", "BB"], "Area_km2": [0.6, 0.6]},
            geometry=[poly_a, poly_b],
            crs=POLAR_PROJ,
        )

        result = top_n_sites_per_region(
            score_da,
            polygons,
            slope_deg=slope_da,
            illumination=illum_da,
            los_visibility=los_da,
            n_per_region=2,
            min_distance_m=240.0,
        )
        assert len(result) == 4
        assert sorted(result["site_id"]) == [1, 2, 3, 4]
        assert set(result["region_name"]) == {"A", "B"}
        # rank_in_region resets per-region.
        for region in ("A", "B"):
            ranks = sorted(result.loc[result["region_name"] == region, "rank_in_region"])
            assert ranks == [1, 2]


class TestParameterValidation:
    def test_non_positive_n_per_region(self) -> None:
        score, slope, illum, los = [_make_grid(np.zeros((5, 5))) for _ in range(4)]
        with pytest.raises(ValueError, match="n_per_region"):
            top_n_sites_per_region(
                score,
                _square_polygon(half_size_m=240.0),
                slope_deg=slope,
                illumination=illum,
                los_visibility=los,
                n_per_region=0,
            )

    def test_non_positive_min_distance(self) -> None:
        score, slope, illum, los = [_make_grid(np.zeros((5, 5))) for _ in range(4)]
        with pytest.raises(ValueError, match="min_distance_m"):
            top_n_sites_per_region(
                score,
                _square_polygon(half_size_m=240.0),
                slope_deg=slope,
                illumination=illum,
                los_visibility=los,
                min_distance_m=0.0,
            )

    def test_shape_mismatch(self) -> None:
        score = _make_grid(np.zeros((5, 5)))
        slope = _make_grid(np.zeros((6, 5)))
        illum = _make_grid(np.zeros((5, 5)))
        los = _make_grid(np.zeros((5, 5)))
        with pytest.raises(ValueError, match="shape mismatch"):
            top_n_sites_per_region(
                score,
                _square_polygon(half_size_m=240.0),
                slope_deg=slope,
                illumination=illum,
                los_visibility=los,
            )

    def test_no_compliant_cells_returns_empty_geodataframe(self) -> None:
        # Polygon exists, all cells fail the slope filter — function
        # must return an empty GeoDataFrame with the right schema, not
        # raise.
        size = 11
        score = np.full((size, size), 0.8, dtype=np.float64)
        slope = np.full((size, size), 30.0, dtype=np.float64)
        illum = np.full((size, size), 0.7, dtype=np.float64)
        los = np.full((size, size), 0.8, dtype=np.float64)
        result = top_n_sites_per_region(
            _make_grid(score, name="score"),
            _square_polygon(half_size_m=1320.0),
            slope_deg=_make_grid(slope, name="slope_deg"),
            illumination=_make_grid(illum, name="illumination"),
            los_visibility=_make_grid(los, name="los_visibility"),
        )
        assert len(result) == 0
        for col in (
            "site_id",
            "region_name",
            "region_code",
            "rank_in_region",
            "score",
            "lat",
            "lon",
            "hls_compliant",
        ):
            assert col in result.columns
