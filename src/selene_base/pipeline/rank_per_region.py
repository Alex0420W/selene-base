"""Per-region NASA HLS-compliant ranking (week 11).

Within each USGS region polygon, applies the four NASA HLS hard
filters (slope ≤ 8°, 100 m buffer, illumination ≥ 33 %, DTE
visibility ≥ 50 %) and returns the top-N highest-scoring sites that
satisfy all of them. Sites are guaranteed inside their named polygon
by construction. See :func:`selene_base.scoring.ranking.top_n_sites_per_region`
for the algorithm.

This pipeline module produces:

- ``data/outputs/per_region/sites.geojson`` — full per-region site
  table with all attributes.
- ``data/outputs/per_region/sites.csv`` — flat human-friendly view.
- ``data/outputs/per_region/per_region_summary.json`` — per-region
  stats (eligible-area fraction, n_sites, best score, etc.) consumed
  by ``selene validate-per-region``.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401
import typer
import xarray as xr
from rasterio.features import geometry_mask
from scipy.ndimage import distance_transform_edt

from selene_base.scoring.ranking import (
    HLS_BUFFER_M,
    HLS_DTE_VISIBILITY_MIN,
    HLS_ILLUMINATION_MIN,
    HLS_SLOPE_MAX_DEG,
    load_sub_scores,
    top_n_sites_per_region,
)
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_PER_REGION_SUBDIR = "per_region"
SCORED_SUBDIR = "scored"


def _open_cog(path: Path) -> xr.DataArray:
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


def _eligible_area_km2_per_region(
    score_map: xr.DataArray,
    regions_polygons: gpd.GeoDataFrame,
    *,
    slope_deg: xr.DataArray,
    illumination: xr.DataArray,
    los_visibility: xr.DataArray,
    hls_slope_max_deg: float,
    hls_buffer_m: float,
    hls_illumination_min: float,
    hls_dte_visibility_min: float,
) -> tuple[dict[str, float], dict[str, float]]:
    """Compute per-region totals: HLS-eligible cells (km²) and total polygon-cell coverage (km²).

    Mirrors the eligibility logic in :func:`top_n_sites_per_region` so the
    summary's eligible-area fraction reflects exactly what the ranker
    sees. ``polygon_area_km2`` is computed *from the rasterised polygon
    on the analysis grid*, not from the polygon's published Area_km2,
    so the ratio gives a true cell-level fraction.
    """
    transform = score_map.rio.transform()
    pixel_size_m = float(abs(transform.a))
    cell_area_km2 = (pixel_size_m * pixel_size_m) / 1e6

    score_arr = score_map.to_numpy().astype(np.float64)
    slope_arr = slope_deg.to_numpy().astype(np.float64)
    illum_arr = illumination.to_numpy().astype(np.float64)
    los_arr = los_visibility.to_numpy().astype(np.float64)
    safe_slope_mask = np.isfinite(slope_arr) & (slope_arr <= hls_slope_max_deg)
    if safe_slope_mask.any():
        distance_to_steep_pix = distance_transform_edt(safe_slope_mask)
        distance_to_steep_m = distance_to_steep_pix.astype(np.float64) * pixel_size_m
    else:
        distance_to_steep_m = np.zeros_like(slope_arr)

    height, width = score_arr.shape
    polygons = regions_polygons.to_crs(score_map.rio.crs)
    eligible_area: dict[str, float] = {}
    polygon_cell_area: dict[str, float] = {}
    for _, region in polygons.iterrows():
        name = str(region.get("Region", region.get("name", "?")))
        polygon_mask = geometry_mask(
            [region.geometry],
            out_shape=(height, width),
            transform=transform,
            invert=True,
            all_touched=True,
        )
        polygon_cell_area[name] = float(polygon_mask.sum()) * cell_area_km2
        eligible_mask = (
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
        eligible_area[name] = float(eligible_mask.sum()) * cell_area_km2
    return eligible_area, polygon_cell_area


def run(
    *,
    score_map_path: Path | None = None,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    per_region_subdir: str = DEFAULT_PER_REGION_SUBDIR,
    n_per_region: int = 10,
    min_distance_km: float = 2.0,
    hls_slope_max_deg: float = HLS_SLOPE_MAX_DEG,
    hls_buffer_m: float = HLS_BUFFER_M,
    hls_illumination_min: float = HLS_ILLUMINATION_MIN,
    hls_dte_visibility_min: float = HLS_DTE_VISIBILITY_MIN,
    echo: Callable[[str], None] = typer.echo,
) -> gpd.GeoDataFrame:
    """Run per-region HLS-compliant ranking and write artefacts.

    Returns:
        The per-region site GeoDataFrame (also written to disk).

    Raises:
        FileNotFoundError: If the aggregate score COG or any of the
            three required HLS-filter inputs (slope, illumination, LOS
            visibility) are missing.
    """
    outputs_dir = Path(outputs_dir)
    per_region_dir = outputs_dir / per_region_subdir
    per_region_dir.mkdir(parents=True, exist_ok=True)
    if score_map_path is None:
        score_map_path = outputs_dir / "score_southpole.tif"

    processed_dir = Path(processed_dir)
    slope_cog = processed_dir / "lola_slope_deg_southpole_240m.tif"
    illum_cog = processed_dir / "illumination_southpole_240m.tif"
    los_cog = processed_dir / "los_visibility_fraction_southpole_240m.tif"

    for label, path in [
        ("aggregate score", score_map_path),
        ("slope", slope_cog),
        ("illumination", illum_cog),
        ("LOS visibility", los_cog),
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"{label} input not found at {path}; run `selene preprocess` "
                "and `selene score` first."
            )

    score_map = _open_cog(score_map_path)
    slope_deg = _open_cog(slope_cog)
    illumination = _open_cog(illum_cog)
    los_visibility = _open_cog(los_cog)

    scored_dir = processed_dir / SCORED_SUBDIR
    sub_scores = load_sub_scores(scored_dir) if scored_dir.exists() else {}

    polygons = regions_polygons_to_geodataframe()

    echo(
        f"[rank-per-region] {len(polygons)} USGS regions, top-{n_per_region} per region, "
        f"min_distance={min_distance_km} km"
    )
    echo(
        f"  HLS filters: slope <= {hls_slope_max_deg} deg, buffer >= {hls_buffer_m:.0f} m, "
        f"illumination >= {hls_illumination_min}, DTE visibility >= {hls_dte_visibility_min}"
    )
    sites = top_n_sites_per_region(
        score_map,
        polygons,
        slope_deg=slope_deg,
        illumination=illumination,
        los_visibility=los_visibility,
        sub_scores=sub_scores,
        n_per_region=n_per_region,
        min_distance_m=min_distance_km * 1000.0,
        hls_slope_max_deg=hls_slope_max_deg,
        hls_buffer_m=hls_buffer_m,
        hls_illumination_min=hls_illumination_min,
        hls_dte_visibility_min=hls_dte_visibility_min,
    )

    eligible_area, polygon_cell_area = _eligible_area_km2_per_region(
        score_map,
        polygons,
        slope_deg=slope_deg,
        illumination=illumination,
        los_visibility=los_visibility,
        hls_slope_max_deg=hls_slope_max_deg,
        hls_buffer_m=hls_buffer_m,
        hls_illumination_min=hls_illumination_min,
        hls_dte_visibility_min=hls_dte_visibility_min,
    )

    geojson_path = per_region_dir / "sites.geojson"
    csv_path = per_region_dir / "sites.csv"
    summary_path = per_region_dir / "per_region_summary.json"

    if len(sites) == 0:
        echo("[rank-per-region] no HLS-compliant sites found in any region")
    sites.to_file(geojson_path, driver="GeoJSON")
    flat = sites.drop(columns="geometry").copy()
    flat.to_csv(csv_path, index=False)
    echo(f"[done] {len(sites)} site(s) -> {geojson_path}")
    echo(f"[done] {len(sites)} site(s) -> {csv_path}")

    summary_rows = []
    for _, region in polygons.iterrows():
        name = str(region.get("Region", "?"))
        code = str(region.get("RegionCode", ""))
        in_region = sites[sites["region_name"] == name] if len(sites) > 0 else sites
        summary_rows.append(
            {
                "name": name,
                "code": code,
                "n_sites": int(len(in_region)),
                "best_score": float(in_region["score"].max()) if len(in_region) > 0 else None,
                "mean_score": float(in_region["score"].mean()) if len(in_region) > 0 else None,
                "eligible_area_km2": float(eligible_area.get(name, 0.0)),
                "polygon_cell_area_km2": float(polygon_cell_area.get(name, 0.0)),
                "eligible_area_fraction": (
                    float(eligible_area.get(name, 0.0) / polygon_cell_area.get(name, 1.0))
                    if polygon_cell_area.get(name, 0.0) > 0
                    else 0.0
                ),
            }
        )
    summary_path.write_text(
        json.dumps(
            {
                "n_sites_total": int(len(sites)),
                "n_regions_total": int(len(polygons)),
                "n_per_region": int(n_per_region),
                "min_distance_km": float(min_distance_km),
                "hls": {
                    "slope_max_deg": float(hls_slope_max_deg),
                    "buffer_m": float(hls_buffer_m),
                    "illumination_min": float(hls_illumination_min),
                    "dte_visibility_min": float(hls_dte_visibility_min),
                },
                "per_region": summary_rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    echo(f"[done] per-region summary -> {summary_path}")
    echo("")
    echo(_format_summary_table(summary_rows))
    return sites


def _format_summary_table(rows: list[dict]) -> str:
    lines = [
        f"{'region':<22} {'code':<4} {'n':>3} {'best':>6} {'mean':>6} {'eligible_pct':>13}",
    ]
    for row in rows:
        n = row["n_sites"]
        best = f"{row['best_score']:>6.3f}" if row["best_score"] is not None else "    --"
        mean = f"{row['mean_score']:>6.3f}" if row["mean_score"] is not None else "    --"
        eligible = row.get("eligible_area_fraction", 0.0) * 100.0
        lines.append(f"{row['name']:<22} {row['code']:<4} {n:>3} {best} {mean} {eligible:>12.2f}%")
    return "\n".join(lines)
