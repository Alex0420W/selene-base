# %% [markdown]
# # Week 10 — USGS polygon validation
#
# Replaces the 15 km disk approximations of NASA's Artemis III
# candidate regions with the USGS-published simplified region envelopes
# (DOI 10.5066/P1MEQ6UK). Generates the headline overlay plot for the
# v1.2.0 README at `data/outputs/sanity/usgs_polygon_validation.png`.
#
# The 9 USGS polygons are simplified 4-vertex envelopes of NASA's LROC
# QuickMap region definitions; together they total ~8000 km² but
# Mons Mouton Plateau alone accounts for 4452 km² (the eight other
# polygons average ~400 km² each). Disk approximations from v1.1.0
# average 707 km² regardless of region — systematically larger for
# the small regions and ~6× smaller than the actual Mons Mouton
# Plateau footprint.
#
# Run after `selene preprocess && selene score && selene rank`.

# %%
from __future__ import annotations

import os

# Suppress pyproj's Earth/Moon celestial-body check. The score COG and
# the USGS GeoJSON both use the lunar sphere R=1737400 m, but they
# describe it through different WKT/PROJ shapes — pyproj refuses to
# operate across them without this override. Setting it here keeps the
# notebook self-contained.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray  # noqa: F401  (registers the .rio accessor)
from matplotlib.patches import Patch

from selene_base.validation.comparison import proximity_analysis
from selene_base.validation.nasa_regions import (
    regions_polygons_to_geodataframe,
    regions_to_geodataframe,
)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)

OUT_DIR = Path("data/outputs/sanity")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SCORE_COG = Path("data/outputs/score_southpole.tif")
TOP_SITES_GEOJSON = Path("data/outputs/top_sites.geojson")

# %%
# Force the score COG's CRS to the canonical lunar polar string. The
# raster is already in this projection — but the COG's WKT writes the
# datum as "unknown" rather than tagging it as Moon, which causes
# pyproj to refuse a transform between the (unknown→Earth-by-default)
# WKT and our explicit lunar CRS. Overwriting the CRS in-memory makes
# both sides agree they're talking about the Moon.
score_da = rioxarray.open_rasterio(SCORE_COG, masked=True).squeeze("band", drop=True)
score_da = score_da.rio.write_crs(POLAR_PROJ, inplace=False)
top_sites = gpd.read_file(TOP_SITES_GEOJSON).to_crs(POLAR_PROJ)
print(f"score raster: {score_da.shape} on {POLAR_PROJ}")
print(f"top sites: {len(top_sites)} loaded")

# %%
# Load both region representations (legacy disks + USGS polygons),
# matched to the raster's CRS for overlay.
disks = regions_to_geodataframe(target_crs=POLAR_PROJ)
usgs = regions_polygons_to_geodataframe(target_crs=POLAR_PROJ)
print(f"disks: {len(disks)} ({list(disks['name'])})")
print(f"USGS polygons: {len(usgs)}")
print(f"USGS total area: {usgs['Area_km2'].sum():.1f} km²")

# %%
# Compute the proximity result so we know which (if any) sites land
# inside a USGS polygon — those get a brighter highlight on the plot.
result = proximity_analysis(
    top_sites.to_crs("+proj=longlat +R=1737400 +no_defs +type=crs"),
    regions_to_geodataframe(),
    nasa_regions_polygons=regions_polygons_to_geodataframe(),
)
print(
    f"sites inside any USGS polygon: {result['sites_inside_any_usgs_polygon']} / {len(top_sites)}"
)
print(
    "USGS regions containing a top site: "
    f"{result['regions_with_top_site_inside_usgs_polygon']} / {result['n_usgs_regions']}"
)
print(
    "median distance from top site to nearest USGS polygon: "
    f"{result['median_distance_to_nearest_usgs_polygon_km']:.1f} km"
)

inside_ids = {row["site_id"] for row in result["per_site_usgs"] if row["inside_any_usgs_polygon"]}

# %%
# Headline overlay plot — score raster + both disk approximations and
# USGS polygons + top-20 sites, with any USGS-polygon-inside sites
# highlighted.
fig, ax = plt.subplots(figsize=(11, 11))
score_da.plot.imshow(
    ax=ax,
    cmap="viridis",
    vmin=0.0,
    vmax=score_da.quantile(0.999).item(),
    add_colorbar=True,
    cbar_kwargs={"label": "aggregate suitability score", "shrink": 0.7},
)

# Legacy disk outlines (light grey, dashed) — context for v1.0.0 / v1.1.0
disks.boundary.plot(ax=ax, color="#cccccc", linewidth=1.0, linestyle="--")

# USGS polygons — solid red outline; Mons Mouton Plateau gets a slightly
# heavier line because it's the giant one.
usgs_named = usgs.copy()
usgs_named["is_plateau"] = usgs_named["Region"] == "Mons Mouton Plateau"
usgs_named[~usgs_named["is_plateau"]].boundary.plot(
    ax=ax, color="#d6232b", linewidth=2.0, linestyle="-"
)
usgs_named[usgs_named["is_plateau"]].boundary.plot(
    ax=ax, color="#d6232b", linewidth=2.5, linestyle="-"
)

# Region name labels at each USGS polygon centroid.
for _, row in usgs.iterrows():
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

# Top sites — yellow circles; sites inside USGS polygons get a brighter
# fill, others stay open.
top_outside = top_sites[~top_sites["site_id"].isin(inside_ids)]
top_inside = top_sites[top_sites["site_id"].isin(inside_ids)]
ax.scatter(
    top_outside.geometry.x,
    top_outside.geometry.y,
    s=24,
    facecolor="none",
    edgecolor="#ffd400",
    linewidth=1.6,
    label=f"top sites outside USGS polys ({len(top_outside)})",
)
if len(top_inside) > 0:
    ax.scatter(
        top_inside.geometry.x,
        top_inside.geometry.y,
        s=64,
        facecolor="#ffd400",
        edgecolor="black",
        linewidth=1.2,
        label=f"top sites INSIDE USGS polys ({len(top_inside)})",
        zorder=10,
    )

# Custom legend.
legend_handles = [
    Patch(facecolor="none", edgecolor="#d6232b", linewidth=2.0, label="USGS polygon (week 10)"),
    Patch(facecolor="none", edgecolor="#cccccc", linewidth=1.0, label="15 km disk (legacy)"),
]
legend_labels = [h.get_label() for h in legend_handles]
ax.legend(
    handles=legend_handles
    + [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="none",
            markeredgecolor="#ffd400",
            markersize=8,
            label="top-20 site",
        )
    ],
    labels=legend_labels + ["top-20 site"],
    loc="lower right",
    fontsize=9,
)
ax.set_title(
    f"USGS-polygon validation: {result['sites_inside_any_usgs_polygon']}/{len(top_sites)} "
    f"top sites inside any USGS polygon  •  "
    f"median distance {result['median_distance_to_nearest_usgs_polygon_km']:.1f} km"
)
ax.set_xlabel("polar stereographic x (m)")
ax.set_ylabel("polar stereographic y (m)")
ax.set_aspect("equal")
fig.tight_layout()

out_path = OUT_DIR / "usgs_polygon_validation.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")
