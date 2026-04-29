"""Seismic criterion — penalises proximity to active lobate scarps.

Civilini et al. (2023) re-located several Apollo-era shallow moonquakes
to within tens of kilometres of young lobate scarps from the Watters
catalog, suggesting active thrust faulting. The score is the per-pixel
distance to the nearest scarp, capped at ``safe_distance_km``.

Two helpers:

- :func:`distance_to_scarps` rasterises the scarp catalog into a
  per-pixel "km to nearest scarp" grid via a KDTree on densified scarp
  vertices.
- :func:`compute` maps that distance grid to a [0, 1] score.

The scarp catalog itself is not yet wired (see
``selene_base.data.download.download_scarps``), but the criterion runs
against any GeoDataFrame the user supplies.

Filled in week 3.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString


def _densify_to_points(
    scarps: gpd.GeoDataFrame,
    target_crs: str,
    *,
    spacing_m: float = 1000.0,
) -> np.ndarray:
    """Reproject scarp geometries and sample dense points along them.

    Returns an (N, 2) array of (x, y) in ``target_crs`` units. Points
    on lines/multilines are sampled every ``spacing_m`` metres so that
    a KDTree query gives accurate distance-to-line.
    """
    transformer = Transformer.from_crs(scarps.crs, target_crs, always_xy=True)
    pts: list[tuple[float, float]] = []
    for geom in scarps.geometry:
        if geom is None or geom.is_empty:
            continue
        if isinstance(geom, MultiLineString):
            lines: list[LineString] = list(geom.geoms)
        elif isinstance(geom, LineString):
            lines = [geom]
        else:  # point or polygon — fall back to centroid
            x, y = transformer.transform(geom.centroid.x, geom.centroid.y)
            pts.append((x, y))
            continue
        for line in lines:
            length = line.length
            if length <= 0:
                continue
            n_steps = max(2, int(np.ceil(length / max(spacing_m, 1e-3))) + 1)
            for s in np.linspace(0.0, length, n_steps):
                px, py = transformer.transform(*line.interpolate(s).coords[0])
                pts.append((px, py))
    if not pts:
        return np.empty((0, 2), dtype=np.float64)
    return np.asarray(pts, dtype=np.float64)


def distance_to_scarps(
    scarps: gpd.GeoDataFrame,
    target_grid: xr.DataArray,
    *,
    chunk_rows: int = 256,
) -> xr.DataArray:
    """Compute per-pixel distance (km) to the nearest scarp feature.

    Reprojects scarp geometries into ``target_grid``'s CRS, densifies
    them into points every 1 km, and runs a KDTree nearest-neighbour
    query for every pixel centre.

    Args:
        scarps: GeoDataFrame with a CRS set; supports point, line, and
            multi-line geometries.
        target_grid: DataArray on the projected analysis grid; values
            are unused.
        chunk_rows: Pixel rows queried at once (memory knob).

    Returns:
        DataArray of distances in kilometres on the same grid as
        ``target_grid``. When the catalog is empty, every pixel
        receives ``+inf``.

    Raises:
        ValueError: If either CRS is missing or the grid lacks ``y``
            and ``x`` dims.
    """
    if target_grid.rio.crs is None:
        raise ValueError("target_grid has no CRS; set one before computing distances")
    if scarps.crs is None:
        raise ValueError("scarps has no CRS; set one before computing distances")
    if "y" not in target_grid.dims or "x" not in target_grid.dims:
        raise ValueError(f"target_grid must have ('y', 'x') dims, got {target_grid.dims!r}")

    target_crs = str(target_grid.rio.crs.to_proj4())
    pts = _densify_to_points(scarps, target_crs)
    height, width = target_grid.sizes["y"], target_grid.sizes["x"]

    if pts.size == 0:
        out_arr = np.full((height, width), np.inf, dtype=np.float32)
    else:
        tree = cKDTree(pts)
        transform = target_grid.rio.transform()
        cols = np.arange(width)
        x_row = transform.c + transform.a * (cols + 0.5)

        out_arr = np.empty((height, width), dtype=np.float32)
        for r0 in range(0, height, chunk_rows):
            r1 = min(r0 + chunk_rows, height)
            rows = np.arange(r0, r1)
            ys = transform.f + transform.e * (rows + 0.5)
            xs = np.tile(x_row, len(rows))
            ys_full = np.repeat(ys, width)
            d, _ = tree.query(np.column_stack([xs, ys_full]), k=1)
            out_arr[r0:r1, :] = (d / 1000.0).reshape(r1 - r0, width)

    out = xr.DataArray(
        out_arr,
        coords=target_grid.transpose("y", "x").coords,
        dims=("y", "x"),
        name="scarp_distance_km",
    )
    return out.rio.write_crs(target_grid.rio.crs, inplace=False)


def compute(
    distance_km: xr.DataArray,
    *,
    safe_distance_km: float = 50.0,
) -> xr.DataArray:
    """Map distance-to-scarp (km) to a [0, 1] safety score.

    ``score = clip(distance_km / safe_distance_km, 0, 1)``: a site
    farther than ``safe_distance_km`` from the nearest scarp scores
    1.0; closer sites lose score linearly.

    Args:
        distance_km: DataArray of distances in km from
            :func:`distance_to_scarps`.
        safe_distance_km: Distance at and above which the score is
            1.0. Must be strictly positive.

    Returns:
        DataArray of [0, 1] scores aligned with ``distance_km``.

    Raises:
        ValueError: If ``safe_distance_km`` is non-positive.
    """
    if safe_distance_km <= 0:
        raise ValueError(f"safe_distance_km must be positive, got {safe_distance_km!r}")

    arr = distance_km.to_numpy().astype(np.float64)
    score = np.clip(arr / safe_distance_km, 0.0, 1.0)
    score = np.where(np.isfinite(arr), score, 1.0)  # +inf (empty catalog) -> safe
    score = np.where(np.isnan(arr), np.nan, score)

    out = xr.DataArray(
        score,
        coords=distance_km.coords,
        dims=distance_km.dims,
        name="seismic_score",
    )
    if distance_km.rio.crs is not None:
        out = out.rio.write_crs(distance_km.rio.crs, inplace=False)
    return out
