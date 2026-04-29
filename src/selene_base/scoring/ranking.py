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
    "thermal",
    "ice",
    "hazard",
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
