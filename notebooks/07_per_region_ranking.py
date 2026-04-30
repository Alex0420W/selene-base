# %% [markdown]
# # Week 11 — per-region ranking with NASA HLS hard filters
#
# Reframes the analysis: rather than ranking globally and validating
# against NASA polygons, we rank *within* each USGS Artemis III polygon
# while applying NASA's published HLS hard-constraint filters as a
# precondition. Sites are guaranteed inside their named polygon by
# construction. Generates the headline overlay plot for the v1.3.0
# README at `data/outputs/sanity/per_region_ranking.png`.
#
# Run after:
# ```
# selene preprocess
# selene score
# selene rank-per-region --n-per-region 3
# ```
#
# HLS thresholds (NASA HLS spec; Gracy & Lee 2024 LPSC #1695; Wueller
# et al. 2026 JGR Planets):
#
# - slope <= 8 deg
# - distance to nearest steeper-than-8-deg cell >= 100 m
# - direct illumination >= 33 %
# - direct-to-Earth visibility >= 50 %

# %%
from __future__ import annotations

import os

# Same lunar/Earth ellipsoid override as week 10's USGS-polygon notebook.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

import json
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401  (registers .rio accessor)
from matplotlib.patches import Patch
from rasterio.features import geometry_mask
from scipy.ndimage import distance_transform_edt

from selene_base.scoring.ranking import (
    HLS_BUFFER_M,
    HLS_DTE_VISIBILITY_MIN,
    HLS_ILLUMINATION_MIN,
    HLS_SLOPE_MAX_DEG,
)
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)

OUT_DIR = Path("data/outputs/sanity")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SCORE_COG = Path("data/outputs/score_southpole.tif")
SLOPE_COG = Path("data/processed/lola_slope_deg_southpole_240m.tif")
ILLUM_COG = Path("data/processed/illumination_southpole_240m.tif")
LOS_COG = Path("data/processed/los_visibility_fraction_southpole_240m.tif")
SITES_GEOJSON = Path("data/outputs/per_region/sites.geojson")
SUMMARY_JSON = Path("data/outputs/per_region/per_region_summary.json")

# %%
# Load rasters; force the polar lunar CRS so pyproj doesn't reject the
# Moon-vs-Earth ellipsoid mismatch when reprojecting overlays.
score_da = rioxarray.open_rasterio(SCORE_COG, masked=True).squeeze("band", drop=True)
score_da = score_da.rio.write_crs(POLAR_PROJ, inplace=False)
slope_da = rioxarray.open_rasterio(SLOPE_COG, masked=True).squeeze("band", drop=True)
slope_da = slope_da.rio.write_crs(POLAR_PROJ, inplace=False)
illum_da = rioxarray.open_rasterio(ILLUM_COG, masked=True).squeeze("band", drop=True)
illum_da = illum_da.rio.write_crs(POLAR_PROJ, inplace=False)
los_da = rioxarray.open_rasterio(LOS_COG, masked=True).squeeze("band", drop=True)
los_da = los_da.rio.write_crs(POLAR_PROJ, inplace=False)

# Derive the global eligibility mask (cells that pass *every* HLS hard
# filter). This is what the per-region ranker selects from inside each
# polygon. The slope-buffer is per-cell distance to the nearest cell
# with slope > 8°.
score_arr = score_da.to_numpy().astype(np.float64)
slope_arr = slope_da.to_numpy().astype(np.float64)
illum_arr = illum_da.to_numpy().astype(np.float64)
los_arr = los_da.to_numpy().astype(np.float64)
pixel_size_m = float(abs(score_da.rio.transform().a))
safe_slope = np.isfinite(slope_arr) & (slope_arr <= HLS_SLOPE_MAX_DEG)
distance_to_steep_m = distance_transform_edt(safe_slope) * pixel_size_m
hls_eligible = (
    np.isfinite(score_arr)
    & np.isfinite(slope_arr)
    & np.isfinite(illum_arr)
    & np.isfinite(los_arr)
    & (slope_arr <= HLS_SLOPE_MAX_DEG)
    & (distance_to_steep_m >= HLS_BUFFER_M)
    & (illum_arr >= HLS_ILLUMINATION_MIN)
    & (los_arr >= HLS_DTE_VISIBILITY_MIN)
)
print(
    f"globally HLS-eligible cells: {hls_eligible.sum():,} / {hls_eligible.size:,} "
    f"({100.0 * hls_eligible.sum() / hls_eligible.size:.2f}%)"
)

# %%
sites = gpd.read_file(SITES_GEOJSON).to_crs(POLAR_PROJ)
polygons = regions_polygons_to_geodataframe(target_crs=POLAR_PROJ)
summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
print(
    f"{summary['n_sites_total']} HLS-compliant sites across "
    f"{sum(1 for r in summary['per_region'] if r['n_sites'] > 0)} / "
    f"{summary['n_regions_total']} USGS regions"
)
for entry in summary["per_region"]:
    print(
        f"  {entry['name']:<22} {entry['code']:<3} n={entry['n_sites']}  "
        f"best={entry['best_score']!s:>8}  "
        f"eligible={entry['eligible_area_fraction'] * 100:.2f}%"
    )

# %%
# Headline overlay: score raster + USGS polygons (red) + HLS-eligible
# cells (green) + per-region top sites (yellow).
fig, ax = plt.subplots(figsize=(11, 11))
score_da.plot.imshow(
    ax=ax,
    cmap="viridis",
    vmin=0.0,
    vmax=score_da.quantile(0.999).item(),
    add_colorbar=True,
    cbar_kwargs={"label": "aggregate suitability score", "shrink": 0.7},
)

# Overlay HLS eligibility as a translucent green mask. We draw it
# only where the eligibility mask is true to keep the underlying
# raster visible elsewhere.
eligible_xy = np.where(hls_eligible)
if eligible_xy[0].size > 0:
    transform = score_da.rio.transform()
    elig_xs = transform.c + transform.a * (eligible_xy[1] + 0.5)
    elig_ys = transform.f + transform.e * (eligible_xy[0] + 0.5)
    ax.scatter(
        elig_xs,
        elig_ys,
        s=0.5,
        c="#4ade80",
        alpha=0.35,
        edgecolors="none",
        rasterized=True,
        label="HLS-eligible cell",
    )

polygons.boundary.plot(ax=ax, color="#d6232b", linewidth=2.0, linestyle="-")
for _, row in polygons.iterrows():
    centroid = row.geometry.centroid
    ax.annotate(
        row["RegionCode"],
        xy=(centroid.x, centroid.y),
        xytext=(0, 6),
        textcoords="offset points",
        ha="center",
        fontsize=9,
        color="#d6232b",
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.7},
    )

ax.scatter(
    sites.geometry.x,
    sites.geometry.y,
    s=64,
    facecolor="#ffd400",
    edgecolor="black",
    linewidth=1.0,
    zorder=10,
    label=f"per-region top sites ({len(sites)})",
)

legend_handles = [
    Patch(facecolor="none", edgecolor="#d6232b", linewidth=2.0, label="USGS polygon"),
    Patch(facecolor="#4ade80", edgecolor="#4ade80", alpha=0.5, label="HLS-eligible cell"),
    plt.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor="#ffd400",
        markeredgecolor="black",
        markersize=8,
        label=f"per-region top site ({len(sites)})",
    ),
]
ax.legend(handles=legend_handles, loc="lower right", fontsize=9)
ax.set_title(
    f"Per-region HLS-compliant ranking: "
    f"{summary['n_sites_total']} sites across "
    f"{sum(1 for r in summary['per_region'] if r['n_sites'] > 0)}/"
    f"{summary['n_regions_total']} USGS regions"
)
ax.set_xlabel("polar stereographic x (m)")
ax.set_ylabel("polar stereographic y (m)")
ax.set_aspect("equal")
fig.tight_layout()
out_path = OUT_DIR / "per_region_ranking.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")

# %%
# Per-region eligibility maps: zoom into each polygon and render the
# HLS-eligible cells inside it, with the chosen sites overlaid.
n_regions = len(polygons)
n_cols = 3
n_rows = (n_regions + n_cols - 1) // n_cols
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
axes = np.atleast_2d(axes)
for ax_idx, (_, region) in enumerate(polygons.iterrows()):
    ax = axes[ax_idx // n_cols, ax_idx % n_cols]
    name = str(region["Region"])
    code = str(region["RegionCode"])
    height, width = score_arr.shape
    polygon_mask = geometry_mask(
        [region.geometry],
        out_shape=(height, width),
        transform=score_da.rio.transform(),
        invert=True,
        all_touched=True,
    )
    rows, cols = np.where(polygon_mask)
    if rows.size == 0:
        ax.set_visible(False)
        continue
    r_min, r_max = rows.min(), rows.max() + 1
    c_min, c_max = cols.min(), cols.max() + 1
    pad = 5
    r_min = max(0, r_min - pad)
    r_max = min(height, r_max + pad)
    c_min = max(0, c_min - pad)
    c_max = min(width, c_max + pad)
    score_zoom = score_arr[r_min:r_max, c_min:c_max]
    eligible_zoom = hls_eligible[r_min:r_max, c_min:c_max]

    transform = score_da.rio.transform()
    extent = (
        transform.c + transform.a * c_min,
        transform.c + transform.a * c_max,
        transform.f + transform.e * r_max,
        transform.f + transform.e * r_min,
    )
    ax.imshow(
        score_zoom,
        cmap="viridis",
        extent=extent,
        origin="upper",
        vmin=0.0,
        vmax=np.nanquantile(score_arr, 0.999),
    )
    eligible_overlay = np.where(eligible_zoom, 1.0, np.nan)
    ax.imshow(
        eligible_overlay,
        cmap="Greens",
        extent=extent,
        origin="upper",
        alpha=0.55,
        vmin=0.0,
        vmax=1.0,
    )
    polygons.iloc[[ax_idx]].boundary.plot(ax=ax, color="#d6232b", linewidth=1.6)

    region_sites = sites[sites["region_name"] == name]
    if len(region_sites) > 0:
        ax.scatter(
            region_sites.geometry.x,
            region_sites.geometry.y,
            s=80,
            facecolor="#ffd400",
            edgecolor="black",
            linewidth=1.0,
            zorder=10,
        )
    summary_entry = next(s for s in summary["per_region"] if s["name"] == name)
    ax.set_title(
        f"{name} ({code})\nn={len(region_sites)}, "
        f"eligible {summary_entry['eligible_area_fraction'] * 100:.2f}%"
    )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
# hide any unused axes
for ax_idx in range(n_regions, n_rows * n_cols):
    axes[ax_idx // n_cols, ax_idx % n_cols].set_visible(False)
fig.tight_layout()
out_path = OUT_DIR / "per_region_eligibility.png"
fig.savefig(out_path, dpi=130, bbox_inches="tight")
print(f"saved {out_path}")
