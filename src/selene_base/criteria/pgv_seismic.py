"""PGV-style seismic shaking attenuation criterion (v2.2+).

Replaces v1.8's distance-to-nearest-scarp criterion (logistic on the
single closest scarp) with a per-cell **cumulative PGV (peak ground
velocity) attenuation kernel** that aggregates contributions from
*every* mapped scarp within a cutoff radius:

.. code-block:: text

    cum_pgv(cell) = Σ_{scarp s within cutoff} exp(-d_s / L)

where ``d_s`` is the cell-to-scarp perpendicular distance (i.e. the
shortest distance from the cell to any point on scarp ``s``'s
polyline) and ``L`` is the attenuation length scale. Defaults track
NASA-aligned seismic-hazard framing:

  - **L = 50 km** — matches the Watters (2024) PSJ Mw 5.3 strong-
    shaking distance and is consistent with Civilini et al. (2023)
    relocated-Apollo-moonquake epicentral uncertainty (~10 km).
  - **cutoff = 5L = 250 km** — captures > 99 % of any single scarp's
    contribution (``exp(-5) ≈ 0.0067``); cells outside the cutoff
    skip the per-scarp distance computation entirely.
  - **activity probability = 1.0** — uniform across all mapped
    scarps. Activity-weighted scoring (using Mishra & Kumar 2022's
    145-scarp dated subset) is deferred to v2.3+.

The cumulative PGV is mapped to a ``[0, 1]`` safety score via a
logistic sigmoid such that low ``cum_pgv`` (far from all scarps)
scores near 1.0 and high ``cum_pgv`` (near multiple scarps) scores
near 0:

.. code-block:: text

    score(cum_pgv) = 1 / (1 + exp((cum_pgv - midpoint) / scale))

The default midpoint anchor is **the single-strong-shaking event**
``cum_pgv = exp(-40 km / 50 km) ≈ 0.4493`` — i.e. a cell experiencing
ground motion equivalent to one Mw 5.3 source at the Watters (2024)
strong-shaking distance scores 0.5. The default scale (0.20) gives
a transition width comparable to the anchor's amplitude, so cells
with substantially more cumulative shaking saturate toward 0 and
cells with substantially less saturate toward 1.

**Why this replaces v1.8's distance-to-nearest-scarp logistic.** The
v1.8 metric collapses *which* scarp is closest into a single distance
threshold and scores accordingly, treating a cell 25 km from one
scarp identically to a cell 25 km from a *cluster* of five scarps.
The PGV kernel correctly accumulates contributions from multiple
nearby scarps — sites in dense south-polar scarp clusters score
lower (riskier) than isolated sites at the same nearest-scarp
distance, which is the right inversion for site-selection.

Per-cell diagnostic components persisted alongside the score:

- ``cum_pgv`` — raw cumulative attenuation sum.
- ``nearest_scarp_distance_km`` — distance to the closest
  contributing scarp (matches v1.8's primary metric for back-
  compatibility / direct comparison).
- ``n_contributing_scarps`` — count of scarps inside the cutoff
  for that cell. Surfaces the multi-scarp-cluster signal that v1.8
  collapsed.

Reference:
- Mishra, S., & Kumar, R. (2022). *Spatial and Temporal Distribution
  of Lobate Scarps in the Lunar South Polar Region.* GRL 49,
  e2022GL098505. (Bundled in-repo since v1.8 — 704 polar-region
  segments at ``criteria/data/scarps_mishra_kumar_2022/``.)
- Watters, T. R. (2024). *Recent active thrust faulting on the Moon.*
  Planetary Science Journal — documents Mw ~5.3 strong-shaking
  PGV at 40 km from active scarps; source of the 50 km L default
  and the 0.45 sigmoid-midpoint anchor.
- Watters, T. R. et al. (2015). *Global thrust faulting on the Moon
  and the influence of tidal stresses.* Geology 43, 851. (Watters
  2015 polar-segment subset bundled at
  ``criteria/data/scarps_watters_2015_polar/`` as the attribution
  anchor for Mishra & Kumar 2022's south-polar extension.)
- Civilini, F. et al. (2023). *Constraints on the seismic hazard of
  young thrust faults on the Moon from re-located shallow
  moonquakes.* (Physical motivation; ~10–50 km moonquake-to-scarp
  clustering distance.)

Filled in v2.2.
"""

from __future__ import annotations

from collections import defaultdict

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString

from selene_base.criteria.seismic import (
    BUNDLED_MISHRA_KUMAR_2022,
    BUNDLED_WATTERS_2015_POLAR,
    load_bundled_catalog,
)

DEFAULT_ATTENUATION_L_KM = 50.0
DEFAULT_CUTOFF_RADIUS_KM = 250.0  # 5L
DEFAULT_DENSIFY_SPACING_M = 1000.0  # matches v1.8's seismic.distance_to_scarps

# Sigmoid calibration anchored to Watters (2024) strong-shaking framing:
# cum_pgv = exp(-40 km / 50 km) ≈ 0.4493 corresponds to one Mw 5.3 source
# at the documented strong-shaking distance, which we score at 0.5.
DEFAULT_SIGMOID_MIDPOINT = 0.45
DEFAULT_SIGMOID_SCALE = 0.20

__all__ = [
    "DEFAULT_ATTENUATION_L_KM",
    "DEFAULT_CUTOFF_RADIUS_KM",
    "DEFAULT_DENSIFY_SPACING_M",
    "DEFAULT_SIGMOID_MIDPOINT",
    "DEFAULT_SIGMOID_SCALE",
    "BUNDLED_MISHRA_KUMAR_2022",
    "BUNDLED_WATTERS_2015_POLAR",
    "load_bundled_catalog",
    "compute_components",
    "compute",
]


def _densify_with_scarp_ids(
    scarps: gpd.GeoDataFrame,
    target_crs: str,
    *,
    spacing_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Reproject + sample dense points along each scarp, tracking source ID.

    Returns ``(points, scarp_ids)``:

    - ``points``: ``(N, 2)`` array of ``(x, y)`` in ``target_crs``
      units.
    - ``scarp_ids``: ``(N,)`` int64 array; the source-frame row index
      for each point. Used downstream to bucket per-cell distances by
      scarp so the per-scarp minimum becomes the contribution
      distance for that scarp.
    """
    transformer = Transformer.from_crs(scarps.crs, target_crs, always_xy=True)
    pts: list[tuple[float, float]] = []
    ids: list[int] = []
    for row_idx, geom in enumerate(scarps.geometry):
        if geom is None or geom.is_empty:
            continue
        if isinstance(geom, MultiLineString):
            lines: list[LineString] = list(geom.geoms)
        elif isinstance(geom, LineString):
            lines = [geom]
        else:  # point/polygon → centroid only
            x, y = transformer.transform(geom.centroid.x, geom.centroid.y)
            pts.append((x, y))
            ids.append(row_idx)
            continue
        for line in lines:
            length = line.length
            if length <= 0:
                continue
            n_steps = max(2, int(np.ceil(length / max(spacing_m, 1e-3))) + 1)
            for s in np.linspace(0.0, length, n_steps):
                px, py = transformer.transform(*line.interpolate(s).coords[0])
                pts.append((px, py))
                ids.append(row_idx)
    if not pts:
        return np.empty((0, 2), dtype=np.float64), np.empty((0,), dtype=np.int64)
    return np.asarray(pts, dtype=np.float64), np.asarray(ids, dtype=np.int64)


def _per_scarp_min_distance(
    cell_xy: np.ndarray,
    tree: cKDTree,
    scarp_ids: np.ndarray,
    cutoff_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-cell: per-scarp minimum distance + cumulative PGV terms.

    For each cell, queries every densified point within ``cutoff_m``
    via :meth:`cKDTree.query_ball_point`, buckets by scarp id, and
    returns the per-scarp minimum distance for the *contributing*
    scarps only (scarps that have ≥ 1 densified point inside the
    cutoff).

    Returns ``(cum_indices, cum_distances_m, n_contributing)``:

    - ``cum_indices``: ``(K,)`` int64; cell index for each
      (cell, scarp) contribution.
    - ``cum_distances_m``: ``(K,)`` float32; the per-(cell, scarp)
      minimum distance in metres.
    - ``n_contributing``: ``(M,)`` int32 per-cell count of
      contributing scarps.

    The caller maps ``cum_distances_m`` to ``exp(-d/L)`` and reduces
    by ``cum_indices`` to compute ``cum_pgv`` per cell.
    """
    n_cells = cell_xy.shape[0]
    n_contributing = np.zeros(n_cells, dtype=np.int32)
    out_indices: list[np.ndarray] = []
    out_distances: list[np.ndarray] = []
    # cKDTree.query_ball_point with return_sorted is sorted by distance,
    # but we still need distances themselves; query+filter is simpler.
    # For per-cell processing, do a per-cell loop with dict bucketing
    # (n_scarps is small, ~700, so dict-per-cell is cheap).
    for i in range(n_cells):
        idxs = tree.query_ball_point(cell_xy[i], r=cutoff_m)
        if not idxs:
            continue
        idxs_arr = np.asarray(idxs, dtype=np.int64)
        d_arr = np.linalg.norm(tree.data[idxs_arr] - cell_xy[i], axis=1)
        sids = scarp_ids[idxs_arr]
        # Per-scarp minimum: bucket by sids, take min(d_arr) per bucket.
        min_per: dict[int, float] = defaultdict(lambda: np.inf)
        for j in range(idxs_arr.size):
            s = int(sids[j])
            d = float(d_arr[j])
            if d < min_per[s]:
                min_per[s] = d
        n_contributing[i] = len(min_per)
        if min_per:
            mins = np.fromiter(min_per.values(), dtype=np.float32, count=len(min_per))
            out_distances.append(mins)
            out_indices.append(np.full(mins.size, i, dtype=np.int64))
    if not out_distances:
        return (
            np.empty((0,), dtype=np.int64),
            np.empty((0,), dtype=np.float32),
            n_contributing,
        )
    return (
        np.concatenate(out_indices),
        np.concatenate(out_distances),
        n_contributing,
    )


def compute_components(
    scarps: gpd.GeoDataFrame,
    target_grid: xr.DataArray,
    *,
    attenuation_l_km: float = DEFAULT_ATTENUATION_L_KM,
    cutoff_radius_km: float = DEFAULT_CUTOFF_RADIUS_KM,
    densify_spacing_m: float = DEFAULT_DENSIFY_SPACING_M,
    sigmoid_midpoint: float = DEFAULT_SIGMOID_MIDPOINT,
    sigmoid_scale: float = DEFAULT_SIGMOID_SCALE,
    chunk_rows: int = 256,
) -> dict[str, xr.DataArray]:
    """Per-cell PGV-style seismic-shaking attenuation score with diagnostics.

    For every pixel of ``target_grid``:

    1. Densify each scarp polyline into points at ``densify_spacing_m``
       along its length, tagged with the originating scarp's row
       index.
    2. Build a single :class:`scipy.spatial.cKDTree` over the
       densified points and query, per cell, every point within
       ``cutoff_radius_km``.
    3. Bucket the in-range densified points by scarp id; the per-
       scarp minimum distance is the contribution distance for that
       scarp.
    4. Sum ``exp(-d_s / L)`` over the contributing scarps to get the
       cell's ``cum_pgv``.
    5. Map ``cum_pgv`` to a [0, 1] safety score via the logistic
       ``score = 1 / (1 + exp((cum_pgv - midpoint) / scale))``.

    Returns a dict with four DataArrays aligned to ``target_grid``:

    - ``"score"`` — combined safety score in [0, 1]; high = safer.
    - ``"cum_pgv"`` — cumulative attenuation sum (raw, dimensionless).
    - ``"nearest_scarp_distance_km"`` — distance to the closest
      contributing scarp (km). Cells with no scarp inside the cutoff
      are ``+inf``; for direct comparison to v1.8's distance-to-
      nearest metric.
    - ``"n_contributing_scarps"`` — count of scarps inside the cutoff
      for each cell. Surfaces the multi-scarp-cluster signal that
      v1.8's distance-to-nearest collapsed.

    Args:
        scarps: GeoDataFrame of mapped lobate-scarp polylines with a
            CRS set. The bundled Mishra & Kumar (2022) catalog at
            :data:`BUNDLED_MISHRA_KUMAR_2022` is the project default
            (loaded via :func:`load_bundled_catalog`).
        target_grid: DataArray on the projected analysis grid; values
            unused, only the grid coords + ``rio.transform()`` /
            ``rio.crs`` matter.
        attenuation_l_km: Characteristic decay length L (km) in
            ``exp(-d/L)``. Default 50 km matches Watters (2024)
            strong-shaking distance.
        cutoff_radius_km: Maximum scarp-to-cell distance (km) at which
            a scarp still contributes. Default 250 km = 5 L (captures
            > 99 % of any single scarp's contribution).
        densify_spacing_m: Sampling spacing along each scarp polyline
            (m). Default 1000 m matches v1.8's seismic.distance_to_
            scarps; finer spacing slows the tree query, coarser
            spacing erodes the per-scarp-min approximation.
        sigmoid_midpoint: ``cum_pgv`` value mapped to score = 0.5.
            Default 0.45 anchors the score at one Mw 5.3 source at
            40 km (Watters 2024).
        sigmoid_scale: Logistic scale for the cum_pgv → score
            mapping. Default 0.20.
        chunk_rows: Pixel rows queried at once (memory knob; matches
            v1.8 seismic.distance_to_scarps).

    Returns:
        Mapping ``{"score", "cum_pgv", "nearest_scarp_distance_km",
        "n_contributing_scarps"}`` to DataArrays on ``target_grid``'s
        grid.

    Raises:
        ValueError: If the grid lacks ``y`` / ``x`` dims, the scarp
            CRS is missing, or any parameter is non-positive.
    """
    if target_grid.rio.crs is None:
        raise ValueError("target_grid has no CRS; set one before computing PGV")
    if scarps.crs is None:
        raise ValueError("scarps has no CRS; set one before computing PGV")
    if "y" not in target_grid.dims or "x" not in target_grid.dims:
        raise ValueError(f"target_grid must have ('y', 'x') dims, got {target_grid.dims!r}")
    if attenuation_l_km <= 0:
        raise ValueError(f"attenuation_l_km must be positive, got {attenuation_l_km!r}")
    if cutoff_radius_km <= 0:
        raise ValueError(f"cutoff_radius_km must be positive, got {cutoff_radius_km!r}")
    if densify_spacing_m <= 0:
        raise ValueError(f"densify_spacing_m must be positive, got {densify_spacing_m!r}")
    if sigmoid_scale <= 0:
        raise ValueError(f"sigmoid_scale must be positive, got {sigmoid_scale!r}")

    target_crs = str(target_grid.rio.crs.to_proj4())
    pts, scarp_ids = _densify_with_scarp_ids(scarps, target_crs, spacing_m=densify_spacing_m)

    height = target_grid.sizes["y"]
    width = target_grid.sizes["x"]

    cum_pgv = np.zeros((height, width), dtype=np.float32)
    nearest_km = np.full((height, width), np.inf, dtype=np.float32)
    n_contrib = np.zeros((height, width), dtype=np.int32)

    if pts.size > 0:
        tree = cKDTree(pts)
        cutoff_m = cutoff_radius_km * 1000.0
        l_m = attenuation_l_km * 1000.0
        transform = target_grid.rio.transform()
        cols = np.arange(width)
        x_row = transform.c + transform.a * (cols + 0.5)

        for r0 in range(0, height, chunk_rows):
            r1 = min(r0 + chunk_rows, height)
            rows = np.arange(r0, r1)
            ys = transform.f + transform.e * (rows + 0.5)
            xs = np.tile(x_row, len(rows))
            ys_full = np.repeat(ys, width)
            cell_xy = np.column_stack([xs, ys_full])

            cum_idx, cum_d, n_per = _per_scarp_min_distance(cell_xy, tree, scarp_ids, cutoff_m)
            n_contrib[r0:r1, :] = n_per.reshape(r1 - r0, width)
            if cum_d.size > 0:
                # Per (cell, scarp) contribution: exp(-d/L)
                contrib = np.exp(-cum_d / l_m).astype(np.float32)
                # Reduce by cell index for cum_pgv; min for nearest.
                chunk_cum = np.zeros(cell_xy.shape[0], dtype=np.float32)
                chunk_nearest_m = np.full(cell_xy.shape[0], np.inf, dtype=np.float32)
                np.add.at(chunk_cum, cum_idx, contrib)
                np.minimum.at(chunk_nearest_m, cum_idx, cum_d)
                cum_pgv[r0:r1, :] = chunk_cum.reshape(r1 - r0, width)
                nearest_km[r0:r1, :] = (chunk_nearest_m / 1000.0).reshape(r1 - r0, width)

    # Sigmoid mapping: low cum_pgv → high score (safe), high cum_pgv → low.
    score = 1.0 / (1.0 + np.exp((cum_pgv.astype(np.float64) - sigmoid_midpoint) / sigmoid_scale))
    score = score.astype(np.float64)
    score = np.clip(score, 0.0, 1.0)

    coords = target_grid.transpose("y", "x").coords

    def _wrap(values: np.ndarray, name: str) -> xr.DataArray:
        out = xr.DataArray(values, coords=coords, dims=("y", "x"), name=name)
        return out.rio.write_crs(target_grid.rio.crs, inplace=False)

    return {
        "score": _wrap(score, "pgv_seismic_score"),
        "cum_pgv": _wrap(cum_pgv.astype(np.float64), "pgv_seismic_cum_pgv"),
        "nearest_scarp_distance_km": _wrap(
            nearest_km.astype(np.float64), "pgv_seismic_nearest_scarp_distance_km"
        ),
        "n_contributing_scarps": _wrap(
            n_contrib.astype(np.int32), "pgv_seismic_n_contributing_scarps"
        ),
    }


def compute(
    scarps: gpd.GeoDataFrame,
    target_grid: xr.DataArray,
    *,
    attenuation_l_km: float = DEFAULT_ATTENUATION_L_KM,
    cutoff_radius_km: float = DEFAULT_CUTOFF_RADIUS_KM,
    densify_spacing_m: float = DEFAULT_DENSIFY_SPACING_M,
    sigmoid_midpoint: float = DEFAULT_SIGMOID_MIDPOINT,
    sigmoid_scale: float = DEFAULT_SIGMOID_SCALE,
    chunk_rows: int = 256,
) -> xr.DataArray:
    """Per-cell PGV-style seismic-shaking score (combined).

    Convenience wrapper over :func:`compute_components` returning only
    the combined ``[0, 1]`` score. Use :func:`compute_components` when
    the per-cell ``cum_pgv``, ``nearest_scarp_distance_km``, or
    ``n_contributing_scarps`` diagnostics are needed.
    """
    return compute_components(
        scarps,
        target_grid,
        attenuation_l_km=attenuation_l_km,
        cutoff_radius_km=cutoff_radius_km,
        densify_spacing_m=densify_spacing_m,
        sigmoid_midpoint=sigmoid_midpoint,
        sigmoid_scale=sigmoid_scale,
        chunk_rows=chunk_rows,
    )["score"]
