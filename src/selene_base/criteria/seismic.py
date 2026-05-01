"""Seismic criterion — penalises proximity to active lobate scarps.

Civilini et al. (2023) re-located several Apollo-era shallow moonquakes
and showed several of them cluster within tens of kilometres of young
lobate scarps from the Watters et al. (2015) catalog, suggesting
active thrust faulting. The score is a per-pixel logistic of distance
to the nearest mapped scarp segment, with the transition tuned to
that documented distance regime (cells inside ~5 km of a mapped scarp
score near 0, cells beyond ~50 km score near 1).

Three helpers:

- :func:`load_bundled_catalog` loads the Mishra & Kumar (2022) primary
  scarp catalog from the in-repo bundle at
  :mod:`selene_base.criteria.data.scarps_mishra_kumar_2022`.
- :func:`distance_to_scarps` rasterises a scarp catalog into a
  per-pixel "km to nearest scarp" grid via a KDTree on densified
  vertices (~1 km spacing along each line).
- :func:`compute` maps that distance grid to a ``[0, 1]`` score via a
  logistic sigmoid (default midpoint 25 km, steepness scale 8 km).

Activated as a contributing criterion in v1.8 (eighth of eight); the
scaffold landed in week 3 against synthetic catalogs.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString

# In-repo scarp catalog bundle paths. The Mishra & Kumar 2022 file is the
# primary data source consumed by `selene score`'s seismic criterion;
# the LROC SOC Watters 2015 polar shapefile sits alongside as an
# attribution anchor (see the bundle READMEs for full provenance).
_BUNDLED_DATA_DIR = Path(__file__).parent / "data"
BUNDLED_MISHRA_KUMAR_2022 = _BUNDLED_DATA_DIR / "scarps_mishra_kumar_2022" / "main_segments.shp"
BUNDLED_WATTERS_2015_POLAR = (
    _BUNDLED_DATA_DIR / "scarps_watters_2015_polar" / "polar_scarp_locations.shp"
)


def load_bundled_catalog(path: Path | None = None) -> gpd.GeoDataFrame:
    """Load the in-repo scarp catalog used by the seismic criterion.

    Defaults to the Mishra & Kumar (2022) 704-segment south-polar
    shapefile bundled at
    :data:`BUNDLED_MISHRA_KUMAR_2022`. Pass an explicit path to use a
    different catalog (e.g. the bundled Watters 2015 polar locations
    for a sanity comparison).

    Args:
        path: Override path. ``None`` uses the bundled Mishra & Kumar
            primary file.

    Returns:
        GeoDataFrame with the catalog's native CRS preserved. The
        seismic criterion's :func:`distance_to_scarps` reprojects to
        the analysis grid's CRS at consumption time.

    Raises:
        FileNotFoundError: If the requested shapefile is missing.
    """
    src = Path(path) if path is not None else BUNDLED_MISHRA_KUMAR_2022
    if not src.exists():
        raise FileNotFoundError(
            f"bundled scarp catalog not found at {src}; "
            "the v1.8 release ships this file in-repo, so a missing "
            "file usually indicates a broken installation."
        )
    return gpd.read_file(src)


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
    midpoint_km: float = 25.0,
    steepness_km: float = 8.0,
) -> xr.DataArray:
    """Map distance-to-scarp (km) to a ``[0, 1]`` safety score via logistic.

    ``score(d) = 1 / (1 + exp(-(d - midpoint_km) / steepness_km))``.

    With the v1.8 defaults ``midpoint_km = 25``, ``steepness_km = 8``:

    - ``d = 5 km``  → 0.08 (well inside Civilini et al. 2023's
      shallow-moonquake-to-scarp clustering distance; high seismic risk)
    - ``d = 15 km`` → 0.22
    - ``d = 25 km`` → 0.50 (transition; matches the typical 10–50 km
      epicentral-uncertainty regime documented for the relocated
      Apollo-era moonquakes)
    - ``d = 35 km`` → 0.78
    - ``d = 50 km`` → 0.96 (effectively safe under the Civilini regime)

    The logistic shape is preferred over a hard linear ramp because the
    seismic-distance evidence is itself probabilistic — Civilini's
    moonquake-to-scarp distances are noisy estimates with epicentral
    uncertainty of order 10 km.

    Cells where the catalog is empty (``+inf`` distance) are scored 1.0;
    NaN distances propagate to NaN scores.

    Args:
        distance_km: DataArray of distances in km from
            :func:`distance_to_scarps`.
        midpoint_km: Distance at which the score crosses 0.5. Must be
            strictly positive.
        steepness_km: Logistic scale (km). Larger values produce a
            gentler transition; smaller values approach a step
            function. Must be strictly positive.

    Returns:
        DataArray of ``[0, 1]`` scores aligned with ``distance_km``.

    Raises:
        ValueError: If either parameter is non-positive.
    """
    if midpoint_km <= 0:
        raise ValueError(f"midpoint_km must be positive, got {midpoint_km!r}")
    if steepness_km <= 0:
        raise ValueError(f"steepness_km must be positive, got {steepness_km!r}")

    arr = distance_km.to_numpy().astype(np.float64)
    # +inf (empty catalog or unreachable cell) → 1.0 directly; the
    # logistic asymptotes there anyway, but skipping the exponent keeps
    # numerical hygiene clean.
    finite_mask = np.isfinite(arr)
    score = np.full_like(arr, 1.0)
    if finite_mask.any():
        d = arr[finite_mask]
        score[finite_mask] = 1.0 / (1.0 + np.exp(-(d - midpoint_km) / steepness_km))
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
