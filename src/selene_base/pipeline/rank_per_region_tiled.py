"""Per-region tiled ranking at high resolution (v1.5).

Pairs with :mod:`selene_base.pipeline.preprocess_tiled` to deliver the
v1.5 "Wueller-class resolution" result: HLS filters applied at 20 m
inside each USGS region polygon, then the ranker picks the top-N by
aggregate score.

Why this lives in a separate module from the v1.4 ranker:

- The 240 m global ranker reads one ``slope_deg``, ``illumination``,
  ``los_visibility`` COG and one ``score_southpole.tif``, all on the
  same global ±304 km grid.
- The tiled ranker has no global 20 m grid (would be ~30 GB per
  criterion); it reads the v1.5 per-tile horizon-profile NPZ produced by
  :func:`selene_base.pipeline.preprocess_tiled.run_tiled_per_region`
  plus the 20 m DEM source, derives slope and Earth-LOS visibility on
  the per-tile grid, resamples the existing 240 m illumination and
  aggregate-score COGs onto the tile grid, applies the four HLS hard
  filters at 20 m, and runs greedy NMS within each polygon.

What v1.5 changes vs v1.4.2:

- The HLS *buffer-from-steep* constraint becomes meaningful. At 240 m,
  ``hls_buffer_m=100`` rounds to ``100/240 = 0.42`` pixels — looser than
  the published threshold. At 20 m it is 5 pixels and the constraint
  actually applies as written.
- The horizon profile is computed from 20 m elevation, so narrow
  horizon-blocking features (massif edges, crater rims) are resolved
  rather than averaged into a 240 m cell.

What v1.5 does *not* change:

- The aggregate score map. The ranker upsamples
  ``data/outputs/score_southpole.tif`` (the v1.4 240 m product) onto
  each tile via bilinear reproject_match. Within a small polygon the
  240 m score map varies slowly, so this upsample is a reasonable
  proxy for "score per buildable 20 m cell"; computing per-criterion
  scores at 20 m globally is out of scope for v1.5.
- The four HLS thresholds — these are NASA's published values and are
  not tuned from v1.5's result.
"""

from __future__ import annotations

import gc
import json
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401
import typer
import xarray as xr
from pyproj import Transformer
from rasterio.features import geometry_mask
from scipy.ndimage import distance_transform_edt
from shapely.geometry import Point

from selene_base.criteria import los_to_earth
from selene_base.criteria.slope import derive_slope_degrees
from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS, LUNAR_SOUTH_POLAR_CRS
from selene_base.pipeline.preprocess_tiled import (
    DEFAULT_BUFFER_M,
    TileSpec,
    _load_lola_source,
    compute_tile_specs,
    horizon_npz_path,
    reproject_to_tile,
    resolve_lola_source,
)
from selene_base.scoring.ranking import (
    DEFAULT_CRITERIA,
    HLS_BUFFER_M,
    HLS_DTE_VISIBILITY_MIN,
    HLS_ILLUMINATION_MIN,
    HLS_SLOPE_MAX_DEG,
)
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PER_REGION_TILED_SUBDIR = "per_region_tiled"
SCORED_SUBDIR = "scored"


@dataclass
class TileRankResult:
    """Per-tile ranking output."""

    region_name: str
    region_code: str
    n_sites: int
    eligible_area_km2: float
    polygon_cell_area_km2: float
    best_score: float | None
    mean_score: float | None
    elapsed_s: float


def _open_global_cog(path: Path) -> xr.DataArray:
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


def _make_tile_template(tile: TileSpec, resolution_m: float, target_crs: str) -> xr.DataArray:
    """Empty DataArray on the tile grid, used as a reproject_match template."""
    height, width = tile.shape(resolution_m)
    xs = np.linspace(tile.xmin + resolution_m / 2, tile.xmax - resolution_m / 2, width)
    ys = np.linspace(tile.ymax - resolution_m / 2, tile.ymin + resolution_m / 2, height)
    template = xr.DataArray(
        np.zeros((height, width), dtype=np.float32),
        dims=("y", "x"),
        coords={"y": ys, "x": xs},
    ).rio.write_crs(target_crs, inplace=False)
    return template


def _resample_global_to_tile_grid(global_da: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    """Reproject ``global_da`` onto the tile template grid (bilinear)."""
    return global_da.rio.reproject_match(template, resampling=1)  # 1 = bilinear


def _load_horizon_profile_npz(
    npz_path: Path,
    tile_template: xr.DataArray,
    target_crs: str,
) -> xr.DataArray:
    """Load a per-tile horizon-profile NPZ and wrap it in an xarray DataArray.

    The NPZ stores raw arrays without CRS metadata; we attach the tile
    template's coords and CRS so the result can be passed straight into
    :func:`compute_earth_visibility_fraction`.
    """
    arr = np.load(npz_path, allow_pickle=False)
    horizon = arr["horizon_profile_deg"]
    azimuth_deg = arr["azimuth_deg"]
    if horizon.shape[1:] != tile_template.shape:
        raise ValueError(
            f"horizon NPZ shape {horizon.shape!r} does not match tile template "
            f"{tile_template.shape!r} — was the NPZ produced at this resolution?"
        )
    da = xr.DataArray(
        horizon,
        dims=("azimuth", "y", "x"),
        coords={
            "azimuth": azimuth_deg,
            "y": tile_template.coords["y"],
            "x": tile_template.coords["x"],
        },
        name="horizon_profile_deg",
    )
    return da.rio.write_crs(target_crs, inplace=False)


def _compute_lat_lon_gamma(
    template: xr.DataArray, target_crs: str
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """Lat / lon / grid-convergence on the tile pixel grid.

    Mirrors the global preprocess (preprocess.py LOS section): compute
    geographic lat/lon by transforming projected (x, y) and grid
    convergence ``γ = atan2(x_p, y_p)`` directly.
    """
    xs = template["x"].to_numpy()
    ys = template["y"].to_numpy()
    xx, yy = np.meshgrid(xs, ys)
    transformer = Transformer.from_crs(
        target_crs, "+proj=longlat +R=1737400 +no_defs +type=crs", always_xy=True
    )
    lons, lats = transformer.transform(xx, yy)
    pixel_lat = xr.DataArray(
        lats.astype(np.float64), dims=template.dims, coords=template.coords
    ).rio.write_crs(target_crs, inplace=False)
    pixel_lon = xr.DataArray(
        lons.astype(np.float64), dims=template.dims, coords=template.coords
    ).rio.write_crs(target_crs, inplace=False)
    gamma = xr.DataArray(
        np.arctan2(xx, yy).astype(np.float64), dims=template.dims, coords=template.coords
    ).rio.write_crs(target_crs, inplace=False)
    return pixel_lat, pixel_lon, gamma


def _free_gpu_memory() -> None:
    try:
        import cupy as cp  # type: ignore[import-not-found]

        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass
    gc.collect()


def process_tile(
    spec: TileSpec,
    *,
    elevation_source: xr.DataArray,
    horizon_npz: Path,
    illumination_global: xr.DataArray,
    score_global: xr.DataArray,
    sub_scores_global: dict[str, xr.DataArray],
    polygon: object,  # shapely Polygon in target_crs
    resolution_m: float,
    target_crs: str,
    n_per_region: int,
    min_distance_m: float,
    hls_slope_max_deg: float,
    hls_buffer_m: float,
    hls_illumination_min: float,
    hls_dte_visibility_min: float,
    echo: Callable[[str], None],
) -> tuple[list[dict[str, object]], TileRankResult]:
    """Window + derive + rank one polygon. Returns site rows + summary."""
    import time

    t0 = time.perf_counter()

    template = _make_tile_template(spec, resolution_m, target_crs)
    height, width = template.shape
    cell_area_km2 = (resolution_m * resolution_m) / 1e6

    echo(
        f"[tile-rank] {spec.region_name} ({spec.region_code}) "
        f"{height}x{width} = {height * width:,} px"
    )

    # ---- per-tile elevation + slope (derived at the analysis resolution) ----
    elev_tile = reproject_to_tile(
        elevation_source,
        spec,
        target_crs=target_crs,
        resolution_m=resolution_m,
        resampling="bilinear",
    ).rename("elevation_m")
    slope_tile = derive_slope_degrees(elev_tile, pixel_size_m=resolution_m)

    # ---- per-tile LOS visibility from the v1.5 horizon NPZ ----
    horizon_tile = _load_horizon_profile_npz(horizon_npz, template, target_crs)
    pixel_lat, pixel_lon, gamma = _compute_lat_lon_gamma(template, target_crs)
    # Try the GPU path first; fall back silently on hosts without CuPy
    # (e.g. CI). The visibility sweep is the dominant CPU cost at 20 m.
    try:
        los_tile = los_to_earth.compute_earth_visibility_fraction(
            horizon_tile, pixel_lat, pixel_lon, gamma, use_gpu=True
        )
    except RuntimeError:
        los_tile = los_to_earth.compute_earth_visibility_fraction(
            horizon_tile, pixel_lat, pixel_lon, gamma
        )
    # Drop the heaviest temporaries (~17 GB float32 horizon profile +
    # ~3 GB lat/lon/gamma) immediately. They are not needed past this
    # point; on GB10 unified memory they otherwise crowd both numpy and
    # cupy out of DRAM during the remaining HLS / NMS work.
    del horizon_tile, pixel_lat, pixel_lon, gamma
    _free_gpu_memory()

    # ---- 240 m global rasters resampled onto the tile ----
    illum_tile = _resample_global_to_tile_grid(illumination_global, template)
    score_tile = _resample_global_to_tile_grid(score_global, template)

    # ---- HLS filters at 20 m ----
    score_arr = score_tile.to_numpy().astype(np.float64)
    slope_arr = slope_tile.to_numpy().astype(np.float64)
    illum_arr = illum_tile.to_numpy().astype(np.float64)
    los_arr = los_tile.to_numpy().astype(np.float64)
    if not (slope_arr.shape == illum_arr.shape == los_arr.shape == score_arr.shape):
        raise ValueError(
            f"shape mismatch on tile {spec.region_code}: "
            f"score={score_arr.shape!r} slope={slope_arr.shape!r} "
            f"illum={illum_arr.shape!r} los={los_arr.shape!r}"
        )

    # Polygon mask + HLS distance-to-steep buffer.
    transform = template.rio.transform()
    polygon_mask = geometry_mask(
        [polygon],
        out_shape=(height, width),
        transform=transform,
        invert=True,
        all_touched=True,
    )
    polygon_cell_area_km2 = float(polygon_mask.sum()) * cell_area_km2

    safe_slope_mask = np.isfinite(slope_arr) & (slope_arr <= hls_slope_max_deg)
    if safe_slope_mask.any():
        distance_to_steep_pix = distance_transform_edt(safe_slope_mask)
        distance_to_steep_m = distance_to_steep_pix.astype(np.float64) * resolution_m
    else:
        distance_to_steep_m = np.zeros_like(slope_arr)

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
    eligible_area_km2 = float(compliant_mask.sum()) * cell_area_km2

    rows: list[dict[str, object]] = []
    if not compliant_mask.any():
        echo(
            f"[tile-rank] {spec.region_code}: no HLS-compliant cells in polygon "
            f"({polygon_cell_area_km2:.1f} km² polygon, {eligible_area_km2:.2f} km² eligible)"
        )
        elapsed = time.perf_counter() - t0
        # Match the late-return cleanup so an empty-polygon tile does not
        # leak cupy + numpy state into the next iteration. (horizon_tile,
        # pixel_lat/lon, gamma were already released after the LOS call.)
        del los_tile, illum_tile, score_tile, slope_tile, elev_tile
        del score_arr, slope_arr, illum_arr, los_arr
        del polygon_mask, safe_slope_mask, distance_to_steep_m, compliant_mask
        _free_gpu_memory()
        return rows, TileRankResult(
            region_name=spec.region_name,
            region_code=spec.region_code,
            n_sites=0,
            eligible_area_km2=eligible_area_km2,
            polygon_cell_area_km2=polygon_cell_area_km2,
            best_score=None,
            mean_score=None,
            elapsed_s=elapsed,
        )

    candidate_rows, candidate_cols = np.where(compliant_mask)
    candidate_scores = score_arr[candidate_rows, candidate_cols]
    order = np.argsort(-candidate_scores, kind="stable")
    candidate_rows = candidate_rows[order]
    candidate_cols = candidate_cols[order]
    candidate_scores = candidate_scores[order]

    min_distance_pix_sq = (min_distance_m / resolution_m) ** 2
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

    rows_arr = np.asarray(accepted_rows)
    cols_arr = np.asarray(accepted_cols)
    xs = transform.c + transform.a * (cols_arr + 0.5) + transform.b * (rows_arr + 0.5)
    ys = transform.f + transform.d * (cols_arr + 0.5) + transform.e * (rows_arr + 0.5)
    transformer = Transformer.from_crs(target_crs, LUNAR_GEOGRAPHIC_CRS, always_xy=True)
    lons, lats = transformer.transform(xs, ys)

    sub_score_tiles: dict[str, np.ndarray] = {}
    for crit, sub_global in sub_scores_global.items():
        sub_tile = _resample_global_to_tile_grid(sub_global, template)
        sub_score_tiles[crit] = sub_tile.to_numpy()

    for k, (r, c, s, x, y, lon, lat) in enumerate(
        zip(accepted_rows, accepted_cols, accepted_scores, xs, ys, lons, lats, strict=True)
    ):
        row: dict[str, object] = {
            "region_name": spec.region_name,
            "region_code": spec.region_code.upper(),  # match v1.4 casing in output
            "rank_in_region": k + 1,
            "score": float(s),
            "lat": float(lat),
            "lon": float(lon),
            "x_m": float(x),
            "y_m": float(y),
            "hls_compliant": True,
            "resolution_m": float(resolution_m),
        }
        for crit in DEFAULT_CRITERIA:
            col = f"score_{crit}"
            if crit in sub_score_tiles:
                v = float(sub_score_tiles[crit][r, c])
                row[col] = v if np.isfinite(v) else float("nan")
            else:
                row[col] = float("nan")
        row["geometry"] = Point(lon, lat)
        rows.append(row)

    elapsed = time.perf_counter() - t0
    echo(
        f"[tile-rank] {spec.region_code}: {len(rows)} site(s), "
        f"eligible {eligible_area_km2:.2f}/{polygon_cell_area_km2:.1f} km² "
        f"({100.0 * eligible_area_km2 / max(polygon_cell_area_km2, 1e-9):.2f} %), "
        f"{elapsed:.1f} s"
    )

    # Drop the per-tile arrays before next region. (horizon_tile,
    # pixel_lat/lon, gamma were released after the LOS call.)
    del los_tile, illum_tile, score_tile, slope_tile, elev_tile
    del score_arr, slope_arr, illum_arr, los_arr
    del polygon_mask, safe_slope_mask, distance_to_steep_m, compliant_mask
    del sub_score_tiles
    _free_gpu_memory()

    return rows, TileRankResult(
        region_name=spec.region_name,
        region_code=spec.region_code,
        n_sites=len(rows),
        eligible_area_km2=eligible_area_km2,
        polygon_cell_area_km2=polygon_cell_area_km2,
        best_score=max(accepted_scores) if accepted_scores else None,
        mean_score=float(np.mean(accepted_scores)) if accepted_scores else None,
        elapsed_s=elapsed,
    )


def run(
    *,
    resolution_m: float = 20.0,
    region_codes: Iterable[str] | None = None,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    score_map_path: Path | None = None,
    target_crs: str = str(LUNAR_SOUTH_POLAR_CRS),
    buffer_m: float = DEFAULT_BUFFER_M,
    n_per_region: int = 10,
    min_distance_km: float = 2.0,
    hls_slope_max_deg: float = HLS_SLOPE_MAX_DEG,
    hls_buffer_m: float = HLS_BUFFER_M,
    hls_illumination_min: float = HLS_ILLUMINATION_MIN,
    hls_dte_visibility_min: float = HLS_DTE_VISIBILITY_MIN,
    per_region_subdir: str = DEFAULT_PER_REGION_TILED_SUBDIR,
    source_path: Path | None = None,
    echo: Callable[[str], None] = typer.echo,
) -> gpd.GeoDataFrame:
    """Per-region tiled ranking at ``resolution_m`` (v1.5 driver).

    For each USGS polygon, derives slope and Earth-LOS visibility on a
    polygon-bbox + ``buffer_m`` tile, applies the four HLS hard filters
    at the high resolution, and runs greedy NMS within the polygon.

    Returns:
        GeoDataFrame of accepted sites across all regions, columns
        identical to :func:`selene_base.pipeline.rank_per_region.run` plus
        a ``resolution_m`` column.
    """
    processed_dir = Path(processed_dir)
    outputs_dir = Path(outputs_dir)
    raw_dir = Path(raw_dir)
    per_region_dir = outputs_dir / per_region_subdir
    per_region_dir.mkdir(parents=True, exist_ok=True)

    if score_map_path is None:
        score_map_path = outputs_dir / "score_southpole.tif"
    illum_cog = processed_dir / "illumination_southpole_240m.tif"
    for label, path in [
        ("aggregate score (240 m)", score_map_path),
        ("illumination (240 m)", illum_cog),
    ]:
        if not Path(path).exists():
            raise FileNotFoundError(
                f"{label} input not found at {path}; run `selene preprocess` "
                "and `selene score` first."
            )

    score_global = _open_global_cog(Path(score_map_path))
    illumination_global = _open_global_cog(illum_cog)

    scored_dir = processed_dir / SCORED_SUBDIR
    sub_scores_global: dict[str, xr.DataArray] = {}
    if scored_dir.exists():
        for crit in DEFAULT_CRITERIA:
            sub_path = scored_dir / f"{crit}_score_southpole_240m.tif"
            if sub_path.exists():
                sub_scores_global[crit] = _open_global_cog(sub_path)

    src_path = (
        source_path
        if source_path is not None
        else resolve_lola_source(raw_dir, prefer_resolution_m=int(round(resolution_m)))
    )
    echo(f"[tile-rank] LOLA source: {src_path}")
    elevation_source = _load_lola_source(src_path)

    polygons = regions_polygons_to_geodataframe(target_crs=target_crs)
    if region_codes is not None:
        wanted = {c.upper() for c in region_codes}
        polygons = polygons[polygons["RegionCode"].str.upper().isin(wanted)]

    specs = compute_tile_specs(
        target_crs=target_crs,
        buffer_m=buffer_m,
        resolution_m=resolution_m,
        region_codes=region_codes,
    )
    spec_by_code = {s.region_code: s for s in specs}

    echo(
        f"[tile-rank] {len(polygons)} region(s) at {resolution_m:g} m, "
        f"top-{n_per_region} per region, NMS={min_distance_km:g} km"
    )
    echo(
        f"  HLS filters: slope <= {hls_slope_max_deg} deg, buffer >= {hls_buffer_m:.0f} m, "
        f"illumination >= {hls_illumination_min}, DTE visibility >= {hls_dte_visibility_min}"
    )

    all_rows: list[dict[str, object]] = []
    summaries: list[TileRankResult] = []
    for _, region in polygons.iterrows():
        code = str(region["RegionCode"]).lower()
        spec = spec_by_code.get(code)
        if spec is None:
            echo(f"[tile-rank] {region['Region']}: missing tile spec (skipped)")
            continue
        npz = horizon_npz_path(processed_dir, resolution_m=resolution_m, region_code=code)
        if not npz.exists():
            echo(
                f"[tile-rank] {region['Region']}: horizon NPZ missing at {npz}; "
                "run `selene preprocess --tiled-per-region --resolution "
                f"{int(round(resolution_m))}` first."
            )
            summaries.append(
                TileRankResult(
                    region_name=str(region["Region"]),
                    region_code=code,
                    n_sites=0,
                    eligible_area_km2=0.0,
                    polygon_cell_area_km2=0.0,
                    best_score=None,
                    mean_score=None,
                    elapsed_s=0.0,
                )
            )
            continue
        try:
            rows, summary = process_tile(
                spec,
                elevation_source=elevation_source,
                horizon_npz=npz,
                illumination_global=illumination_global,
                score_global=score_global,
                sub_scores_global=sub_scores_global,
                polygon=region.geometry,
                resolution_m=resolution_m,
                target_crs=target_crs,
                n_per_region=n_per_region,
                min_distance_m=min_distance_km * 1000.0,
                hls_slope_max_deg=hls_slope_max_deg,
                hls_buffer_m=hls_buffer_m,
                hls_illumination_min=hls_illumination_min,
                hls_dte_visibility_min=hls_dte_visibility_min,
                echo=echo,
            )
        except Exception as exc:
            # Don't let a single corrupt NPZ or transient OOM throw away
            # the work already done on the other 8 tiles. Record the
            # failure in the summary and continue.
            echo(f"[tile-rank] {region['Region']} ({code}): FAILED — {type(exc).__name__}: {exc}")
            summaries.append(
                TileRankResult(
                    region_name=str(region["Region"]),
                    region_code=code,
                    n_sites=0,
                    eligible_area_km2=0.0,
                    polygon_cell_area_km2=0.0,
                    best_score=None,
                    mean_score=None,
                    elapsed_s=0.0,
                )
            )
            continue
        finally:
            # Defence in depth: even if process_tile raises, drop the
            # cupy memory pool before moving to the next polygon. On GB10
            # unified memory, GPU allocations come from the same DRAM as
            # numpy, so leaked GPU blocks crash the next tile via OOM.
            _free_gpu_memory()
        all_rows.extend(rows)
        summaries.append(summary)

    # Re-id sites globally (1-indexed across all regions) so the output
    # schema matches the v1.4 per-region table.
    for i, row in enumerate(all_rows, start=1):
        row["site_id"] = i

    if all_rows:
        df = pd.DataFrame(all_rows)
        sites = gpd.GeoDataFrame(df, geometry="geometry", crs=LUNAR_GEOGRAPHIC_CRS)
    else:
        sites = gpd.GeoDataFrame(
            pd.DataFrame(columns=["site_id", "region_name", "region_code", "geometry"]),
            geometry=[],
            crs=LUNAR_GEOGRAPHIC_CRS,
        )

    geojson_path = per_region_dir / "sites.geojson"
    csv_path = per_region_dir / "sites.csv"
    summary_path = per_region_dir / "per_region_summary.json"

    sites.to_file(geojson_path, driver="GeoJSON")
    if len(sites) > 0:
        sites.drop(columns="geometry").to_csv(csv_path, index=False)
    else:
        Path(csv_path).write_text("", encoding="utf-8")

    summary_doc = {
        "n_sites_total": int(len(sites)),
        "n_regions_total": int(len(polygons)),
        "n_per_region": int(n_per_region),
        "min_distance_km": float(min_distance_km),
        "resolution_m": float(resolution_m),
        "buffer_m": float(buffer_m),
        "hls": {
            "slope_max_deg": float(hls_slope_max_deg),
            "buffer_m": float(hls_buffer_m),
            "illumination_min": float(hls_illumination_min),
            "dte_visibility_min": float(hls_dte_visibility_min),
        },
        "per_region": [
            {
                "name": s.region_name,
                "code": s.region_code.upper(),
                "n_sites": s.n_sites,
                "best_score": s.best_score,
                "mean_score": s.mean_score,
                "eligible_area_km2": s.eligible_area_km2,
                "polygon_cell_area_km2": s.polygon_cell_area_km2,
                "eligible_area_fraction": (
                    s.eligible_area_km2 / s.polygon_cell_area_km2
                    if s.polygon_cell_area_km2 > 0
                    else 0.0
                ),
                "elapsed_s": s.elapsed_s,
            }
            for s in summaries
        ],
    }
    summary_path.write_text(json.dumps(summary_doc, indent=2), encoding="utf-8")
    echo(f"[done] {len(sites)} site(s) -> {geojson_path}")
    echo(f"[done] per-region summary -> {summary_path}")
    echo("")
    echo(_format_summary_table(summaries))
    return sites


def _format_summary_table(summaries: list[TileRankResult]) -> str:
    header = (
        f"{'region':<22} {'code':<4} {'n':>3} {'best':>6} {'mean':>6} {'elig_pct':>10} {'sec':>7}"
    )
    lines = [header]
    for s in summaries:
        best = f"{s.best_score:>6.3f}" if s.best_score is not None else "    --"
        mean = f"{s.mean_score:>6.3f}" if s.mean_score is not None else "    --"
        pct = (
            100.0 * s.eligible_area_km2 / s.polygon_cell_area_km2
            if s.polygon_cell_area_km2 > 0
            else 0.0
        )
        lines.append(
            f"{s.region_name:<22} {s.region_code.upper():<4} {s.n_sites:>3} "
            f"{best} {mean} {pct:>9.2f}% {s.elapsed_s:>7.1f}"
        )
    return "\n".join(lines)
