"""Extract geographically-distinct top-N sites from an aggregate score map.

Implements non-maximum suppression in projected metres: pick the highest
remaining cell, emit it, blank a circular ``min_distance_m`` neighbourhood,
repeat until ``n`` sites are collected or the map is exhausted. Then
attach lat/lon and per-criterion sub-scores so the GeoDataFrame is the
single source of truth for the rank step.

Filled in week 3.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import Transformer
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

DEFAULT_CRITERIA = (
    "slope",
    "illumination",
    "coupling",
    "thermal",
    "ice",
    "hazard",
    "los_to_earth",
    "seismic",
)


def _index_to_world(
    score_map: xr.DataArray,
    rows: np.ndarray,
    cols: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert pixel (row, col) indices to projected (x, y) coordinates."""
    transform = score_map.rio.transform()
    xs = transform.c + transform.a * (cols + 0.5) + transform.b * (rows + 0.5)
    ys = transform.f + transform.d * (cols + 0.5) + transform.e * (rows + 0.5)
    return xs, ys


def top_n_sites(
    score_map: xr.DataArray,
    n: int = 20,
    *,
    min_distance_m: float = 5000.0,
    min_score: float = 0.5,
    sub_scores: Mapping[str, xr.DataArray] | None = None,
) -> gpd.GeoDataFrame:
    """Extract the top-N highest-scoring sites separated by ``min_distance_m``.

    Algorithm:
    1. Drop NaN and below-``min_score`` pixels.
    2. Sort the survivors by score, descending.
    3. Walk the sorted list; for each candidate, accept it if it sits
       outside ``min_distance_m`` of every already-accepted site.
    4. Stop at ``n`` accepted sites or when the list is exhausted.

    For each accepted site, attach lat/lon (geographic CRS, R=1737400)
    and per-criterion sub-scores read from the optional ``sub_scores``
    mapping (criteria not present in the mapping become NaN columns).

    Args:
        score_map: Aggregated [0, 1] suitability scores in a projected
            CRS with ``rio.transform()`` set.
        n: Maximum number of sites to return; must be positive.
        min_distance_m: Minimum pairwise distance between returned
            sites, in metres; must be positive.
        min_score: Floor on candidate score; sites below this never
            enter the running.
        sub_scores: Optional mapping from criterion name to the
            criterion's score grid (same shape as ``score_map``);
            sampled at each accepted site to populate the
            ``score_<criterion>`` columns.

    Returns:
        GeoDataFrame in :data:`selene_base.data.load.LUNAR_GEOGRAPHIC_CRS`
        with columns:

        ``site_id``, ``geometry`` (Point), ``score``, ``lat``, ``lon``,
        ``x_m``, ``y_m``, ``score_slope``, ``score_illumination``,
        ``score_thermal``, ``score_ice``, ``score_hazard``,
        ``score_seismic``.

    Raises:
        ValueError: On non-positive ``n`` or ``min_distance_m``, or if
            the score map has no projected CRS.
    """
    if n <= 0:
        raise ValueError(f"n must be positive, got {n!r}")
    if min_distance_m <= 0:
        raise ValueError(f"min_distance_m must be positive, got {min_distance_m!r}")
    if score_map.rio.crs is None:
        raise ValueError("score_map has no CRS; set one before ranking")

    arr = score_map.to_numpy().astype(np.float64)
    valid = np.isfinite(arr) & (arr >= min_score)
    if not valid.any():
        return _empty_geodataframe()

    rows_all, cols_all = np.where(valid)
    scores_all = arr[rows_all, cols_all]
    order = np.argsort(-scores_all, kind="stable")
    rows_all = rows_all[order]
    cols_all = cols_all[order]
    scores_all = scores_all[order]

    transform = score_map.rio.transform()
    pixel_size = float(abs(transform.a))
    min_distance_pix = (min_distance_m / pixel_size) ** 2  # squared, for comparisons

    selected_rows: list[int] = []
    selected_cols: list[int] = []
    selected_scores: list[float] = []

    for r, c, s in zip(rows_all, cols_all, scores_all, strict=True):
        keep = True
        for sr, sc in zip(selected_rows, selected_cols, strict=True):
            dr = float(r - sr)
            dc = float(c - sc)
            if dr * dr + dc * dc < min_distance_pix:
                keep = False
                break
        if keep:
            selected_rows.append(int(r))
            selected_cols.append(int(c))
            selected_scores.append(float(s))
            if len(selected_rows) == n:
                break

    if not selected_rows:
        return _empty_geodataframe()

    rows_arr = np.asarray(selected_rows)
    cols_arr = np.asarray(selected_cols)
    xs, ys = _index_to_world(score_map, rows_arr, cols_arr)

    transformer = Transformer.from_crs(score_map.rio.crs, LUNAR_GEOGRAPHIC_CRS, always_xy=True)
    lons, lats = transformer.transform(xs, ys)

    table: dict[str, list[object]] = {
        "site_id": [f"site_{i + 1:02d}" for i in range(len(selected_rows))],
        "rank": list(range(1, len(selected_rows) + 1)),
        "score": list(selected_scores),
        "lat": list(lats),
        "lon": list(lons),
        "x_m": [float(v) for v in xs],
        "y_m": [float(v) for v in ys],
    }

    sub_scores = sub_scores or {}
    for crit in DEFAULT_CRITERIA:
        col = f"score_{crit}"
        if crit in sub_scores:
            sub = sub_scores[crit].to_numpy()
            try:
                values = sub[rows_arr, cols_arr]
            except IndexError:
                values = np.full(len(selected_rows), np.nan)
            table[col] = [float(v) if np.isfinite(v) else np.nan for v in values]
        else:
            table[col] = [float("nan")] * len(selected_rows)

    geometries = [Point(lon, lat) for lon, lat in zip(lons, lats, strict=True)]
    gdf = gpd.GeoDataFrame(table, geometry=geometries, crs=LUNAR_GEOGRAPHIC_CRS)
    return gdf


def _empty_geodataframe() -> gpd.GeoDataFrame:
    cols = ["site_id", "rank", "score", "lat", "lon", "x_m", "y_m"]
    cols += [f"score_{c}" for c in DEFAULT_CRITERIA]
    df = pd.DataFrame(columns=cols)
    return gpd.GeoDataFrame(df, geometry=[], crs=LUNAR_GEOGRAPHIC_CRS)


def load_sub_scores(scored_dir: Path) -> dict[str, xr.DataArray]:
    """Load whichever per-criterion score COGs are present in ``scored_dir``."""
    out: dict[str, xr.DataArray] = {}
    for crit in DEFAULT_CRITERIA:
        path = scored_dir / f"{crit}_score_southpole_240m.tif"
        if path.exists():
            out[crit] = rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)
    return out


# ---------------------------------------------------------------------------
# Per-region ranking with NASA HLS hard-constraint filters (week 11).
# ---------------------------------------------------------------------------

# NASA Human Landing System hard-constraint thresholds.
# Sources:
#   - NASA HLS specification (NASA 2019)
#   - Gracy & Lee 2024, "Update on the Artemis III Reference Mission",
#     LPSC Abstract #1695.
#   - Wueller, L., et al. 2026. "A Catalog of 130 Candidate Landing Sites
#     for Artemis III", JGR Planets, doi:10.1029/2025JE009434.
HLS_SLOPE_MAX_DEG = 8.0
HLS_BUFFER_M = 100.0
HLS_ILLUMINATION_MIN = 0.33
HLS_DTE_VISIBILITY_MIN = 0.50


def top_n_sites_per_region(
    score_map: xr.DataArray,
    regions_polygons: gpd.GeoDataFrame,
    *,
    slope_deg: xr.DataArray,
    illumination: xr.DataArray,
    los_visibility: xr.DataArray,
    sub_scores: Mapping[str, xr.DataArray] | None = None,
    n_per_region: int = 3,
    min_distance_m: float = 2000.0,
    hls_slope_max_deg: float = HLS_SLOPE_MAX_DEG,
    hls_buffer_m: float = HLS_BUFFER_M,
    hls_illumination_min: float = HLS_ILLUMINATION_MIN,
    hls_dte_visibility_min: float = HLS_DTE_VISIBILITY_MIN,
) -> gpd.GeoDataFrame:
    """For each NASA region polygon, find top-N HLS-compliant landing sites.

    Hard filters applied within each polygon (multiplicative AND):

    1. ``slope_deg <= hls_slope_max_deg`` — the landing pad itself is
       buildable.
    2. Distance from the cell to the nearest cell with
       ``slope_deg > hls_slope_max_deg`` is at least ``hls_buffer_m`` —
       the lander has buffer before any tip-over hazard. Computed via
       :func:`scipy.ndimage.distance_transform_edt` on the
       ``slope <= max`` mask.
    3. ``illumination >= hls_illumination_min`` — direct solar
       illumination ≥ 33 % over the representative period.
    4. ``los_visibility >= hls_dte_visibility_min`` — direct-to-Earth
       visibility ≥ 50 % over the libration cycle.

    These thresholds are NASA's published values (HLS spec; Gracy & Lee
    2024 LPSC #1695); they are not tuneable from the validation result.

    After hard filtering, the surviving cells inside each polygon are
    ranked by ``score_map``. Up to ``n_per_region`` sites are emitted
    via greedy NMS at ``min_distance_m`` separation. If a polygon has
    no compliant cells, the region produces zero sites — the result
    GeoDataFrame simply has no rows for that region.

    The buffer distance transform is computed on the cells inside the
    polygon's bounding box (with a small margin) rather than globally,
    so each region's buffer is local to its terrain. Cells near the
    grid edge may have artificially small buffer values; this is not a
    correctness issue for the south polar polygons, all of which sit
    interior to the analysis grid.

    Args:
        score_map: Aggregate ``[0, 1]`` suitability map in a projected
            CRS with ``rio.transform()`` set.
        regions_polygons: GeoDataFrame of NASA region polygons with at
            minimum a ``Region`` column (and ``RegionCode`` if present).
        slope_deg: Slope-in-degrees raster aligned with ``score_map``.
        illumination: Illumination-fraction raster aligned with
            ``score_map``.
        los_visibility: Earth visibility-fraction raster aligned with
            ``score_map``.
        sub_scores: Optional per-criterion score grids; sampled at each
            accepted site to populate the ``score_<crit>`` columns.
        n_per_region: Maximum sites to return per polygon. Default 3
            matches the Wueller et al. 2026 / NASA Artemis-III selection
            cadence (a small handful of candidate landing pads per
            region).
        min_distance_m: Minimum pairwise distance between sites within
            the same region, in metres. Default 2 km — tighter than the
            global rank's 5 km because each polygon is small (~400 km²
            for most regions); 2 km still gives non-overlapping
            ~100 m landing pads.
        hls_slope_max_deg, hls_buffer_m, hls_illumination_min,
            hls_dte_visibility_min: HLS thresholds; documented above.

    Returns:
        GeoDataFrame with one row per accepted site, columns:
        ``site_id`` (1-indexed across all regions), ``region_name``,
        ``region_code``, ``rank_in_region``, ``score``, ``lat``,
        ``lon``, ``x_m``, ``y_m``, ``hls_compliant`` (always True),
        plus a ``score_<crit>`` column per criterion in
        ``DEFAULT_CRITERIA``. The frame is in
        :data:`selene_base.data.load.LUNAR_GEOGRAPHIC_CRS`.

    Raises:
        ValueError: On non-positive parameters or a missing CRS on
            ``score_map``.
    """
    if n_per_region <= 0:
        raise ValueError(f"n_per_region must be positive, got {n_per_region!r}")
    if min_distance_m <= 0:
        raise ValueError(f"min_distance_m must be positive, got {min_distance_m!r}")
    if hls_slope_max_deg <= 0:
        raise ValueError(f"hls_slope_max_deg must be positive, got {hls_slope_max_deg!r}")
    if hls_buffer_m < 0:
        raise ValueError(f"hls_buffer_m must be non-negative, got {hls_buffer_m!r}")
    if score_map.rio.crs is None:
        raise ValueError("score_map has no CRS; set one before ranking")

    from rasterio.features import geometry_mask
    from scipy.ndimage import distance_transform_edt

    raster_crs = score_map.rio.crs
    polygons = regions_polygons.to_crs(raster_crs)

    score_arr = score_map.to_numpy().astype(np.float64)
    slope_arr = slope_deg.to_numpy().astype(np.float64)
    illum_arr = illumination.to_numpy().astype(np.float64)
    los_arr = los_visibility.to_numpy().astype(np.float64)
    if not (slope_arr.shape == illum_arr.shape == los_arr.shape == score_arr.shape):
        raise ValueError(
            f"shape mismatch: score={score_arr.shape!r} slope={slope_arr.shape!r} "
            f"illum={illum_arr.shape!r} los={los_arr.shape!r}"
        )

    transform = score_map.rio.transform()
    pixel_size_m = float(abs(transform.a))
    height, width = score_arr.shape

    # Slope-buffer mask is the same global field for every region — a
    # cell's buffer to the nearest steep cell doesn't depend on which
    # polygon we are looking at — so compute it once. distance_transform_edt
    # returns distance in pixels to the nearest *zero*; we want distance
    # to the nearest cell where slope > threshold, i.e. invert that mask.
    safe_slope_mask = np.isfinite(slope_arr) & (slope_arr <= hls_slope_max_deg)
    if safe_slope_mask.any():
        distance_to_steep_pix = distance_transform_edt(safe_slope_mask)
        distance_to_steep_m = distance_to_steep_pix.astype(np.float64) * pixel_size_m
    else:
        distance_to_steep_m = np.zeros_like(slope_arr)

    rows: list[dict[str, object]] = []
    next_site_id = 1
    transformer = Transformer.from_crs(raster_crs, LUNAR_GEOGRAPHIC_CRS, always_xy=True)
    sub_scores = sub_scores or {}

    for _, region in polygons.iterrows():
        polygon = region.geometry
        region_name = str(region.get("Region", region.get("name", "?")))
        region_code = str(region.get("RegionCode", region.get("code", "")))

        polygon_mask = geometry_mask(
            [polygon],
            out_shape=(height, width),
            transform=transform,
            invert=True,  # True inside the polygon
            all_touched=True,
        )
        if not polygon_mask.any():
            continue

        compliant_mask = (
            polygon_mask
            & np.isfinite(score_arr)
            & np.isfinite(slope_arr)
            & np.isfinite(illum_arr)
            & np.isfinite(los_arr)
            & (slope_arr <= hls_slope_max_deg)
            & (distance_to_steep_m >= hls_buffer_m)
            & (illum_arr >= hls_illumination_min)
            & (los_arr >= hls_dte_visibility_min)
        )
        if not compliant_mask.any():
            continue

        candidate_rows, candidate_cols = np.where(compliant_mask)
        candidate_scores = score_arr[candidate_rows, candidate_cols]
        order = np.argsort(-candidate_scores, kind="stable")
        candidate_rows = candidate_rows[order]
        candidate_cols = candidate_cols[order]
        candidate_scores = candidate_scores[order]

        # NMS within the region.
        min_distance_pix_sq = (min_distance_m / pixel_size_m) ** 2
        accepted_rows: list[int] = []
        accepted_cols: list[int] = []
        accepted_scores: list[float] = []
        for r, c, s in zip(candidate_rows, candidate_cols, candidate_scores, strict=True):
            keep = True
            for ar, ac in zip(accepted_rows, accepted_cols, strict=True):
                dr = float(r - ar)
                dc = float(c - ac)
                if dr * dr + dc * dc < min_distance_pix_sq:
                    keep = False
                    break
            if keep:
                accepted_rows.append(int(r))
                accepted_cols.append(int(c))
                accepted_scores.append(float(s))
                if len(accepted_rows) == n_per_region:
                    break

        if not accepted_rows:
            continue

        rows_arr = np.asarray(accepted_rows)
        cols_arr = np.asarray(accepted_cols)
        xs, ys = _index_to_world(score_map, rows_arr, cols_arr)
        lons, lats = transformer.transform(xs, ys)

        for k, (r, c, s, x, y, lon, lat) in enumerate(
            zip(
                accepted_rows,
                accepted_cols,
                accepted_scores,
                xs,
                ys,
                lons,
                lats,
                strict=True,
            )
        ):
            row: dict[str, object] = {
                "site_id": next_site_id,
                "region_name": region_name,
                "region_code": region_code,
                "rank_in_region": k + 1,
                "score": float(s),
                "lat": float(lat),
                "lon": float(lon),
                "x_m": float(x),
                "y_m": float(y),
                "hls_compliant": True,
            }
            for crit in DEFAULT_CRITERIA:
                col = f"score_{crit}"
                if crit in sub_scores:
                    sub = sub_scores[crit].to_numpy()
                    try:
                        v = float(sub[r, c])
                    except IndexError:
                        v = float("nan")
                    row[col] = v if np.isfinite(v) else float("nan")
                else:
                    row[col] = float("nan")
            row["geometry"] = Point(lon, lat)
            rows.append(row)
            next_site_id += 1

    if not rows:
        return _empty_per_region_geodataframe()

    df = pd.DataFrame(rows)
    return gpd.GeoDataFrame(df, geometry="geometry", crs=LUNAR_GEOGRAPHIC_CRS)


def _empty_per_region_geodataframe() -> gpd.GeoDataFrame:
    cols = [
        "site_id",
        "region_name",
        "region_code",
        "rank_in_region",
        "score",
        "lat",
        "lon",
        "x_m",
        "y_m",
        "hls_compliant",
    ]
    cols += [f"score_{c}" for c in DEFAULT_CRITERIA]
    df = pd.DataFrame(columns=cols)
    return gpd.GeoDataFrame(df, geometry=[], crs=LUNAR_GEOGRAPHIC_CRS)
