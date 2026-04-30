# %% [markdown]
# # Week 12 — quantitative comparison vs Wueller et al. 2026 (v1.4.0 framework)
#
# Wueller, F., et al. (2026), JGR Planets, doi:10.1029/2025JE009434,
# published 130 candidate Artemis III landing sites identified by the
# same outer methodology selene-base implements (NASA HLS hard filters
# followed by within-region selection). v1.4.0 ships the comparison
# **framework**; the real Wueller catalog is currently gated behind
# AGU/Wiley and no open data release has been located, so the bundled
# `wueller_2026_sites.csv` is a synthetic 5-row placeholder.
#
# This notebook produces the v1.4 visualisation set against whichever
# CSV is at `src/selene_base/validation/data/wueller_2026_sites.csv` —
# placeholder or real. The placeholder run renders cleanly and is
# explicitly labelled "SYNTHETIC PLACEHOLDER" on every plot; once the
# real CSV ships, rerunning this notebook regenerates the figures with
# real numbers without code changes.
#
# Run after:
#
# ```
# selene preprocess && selene score
# selene rank-per-region --n-per-region 3
# ```

# %%
from __future__ import annotations

import os

# Same lunar/Earth ellipsoid override as the week 10 / 11 notebooks.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401  (registers .rio accessor)
from matplotlib.lines import Line2D

from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe
from selene_base.validation.wueller_comparison import (
    DEFAULT_MATCH_THRESHOLD_KM,
    compare_sites,
    is_synthetic_placeholder,
    load_wueller_sites,
)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)

OUT_DIR = Path("data/outputs/sanity")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SELENE_SITES_GEOJSON = Path("data/outputs/per_region/sites.geojson")
SCORE_COG = Path("data/outputs/score_southpole.tif")

# %%
selene_sites = gpd.read_file(SELENE_SITES_GEOJSON).to_crs(POLAR_PROJ)
wueller_sites = load_wueller_sites(target_crs=POLAR_PROJ)
polygons = regions_polygons_to_geodataframe(target_crs=POLAR_PROJ)
result = compare_sites(
    selene_sites.to_crs("+proj=longlat +R=1737400 +no_defs +type=crs"),
    wueller_sites.to_crs("+proj=longlat +R=1737400 +no_defs +type=crs"),
    match_threshold_km=DEFAULT_MATCH_THRESHOLD_KM,
)
USING_PLACEHOLDER = result["using_synthetic_placeholder"]
print(
    f"selene={result['n_selene_sites']}, wueller={result['n_wueller_sites']}, "
    f"matched(selene)={result['n_selene_matched']}, "
    f"matched(wueller)={result['n_wueller_matched']}, "
    f"placeholder={USING_PLACEHOLDER}"
)

# %%
# Score raster optional — fall back to a plain map background if missing.
try:
    score_da = (
        rioxarray.open_rasterio(SCORE_COG, masked=True)
        .squeeze("band", drop=True)
        .rio.write_crs(POLAR_PROJ, inplace=False)
    )
except Exception:  # noqa: BLE001 — best-effort backdrop
    score_da = None

# %%
# Headline overlay: USGS polygons + selene sites (yellow) + Wueller sites
# (red) + lines connecting matched pairs.
fig, ax = plt.subplots(figsize=(11, 11))
if score_da is not None:
    score_da.plot.imshow(
        ax=ax,
        cmap="viridis",
        vmin=0.0,
        vmax=score_da.quantile(0.999).item(),
        add_colorbar=True,
        cbar_kwargs={"label": "aggregate suitability score", "shrink": 0.7},
    )
polygons.boundary.plot(ax=ax, color="#d6232b", linewidth=1.5)
ax.scatter(
    selene_sites.geometry.x,
    selene_sites.geometry.y,
    s=64,
    facecolor="#ffd400",
    edgecolor="black",
    linewidth=1.0,
    zorder=10,
)
ax.scatter(
    wueller_sites.geometry.x,
    wueller_sites.geometry.y,
    s=64,
    facecolor="#3b82f6",
    edgecolor="black",
    linewidth=1.0,
    marker="^",
    zorder=10,
)

# Connect matched pairs with thin grey segments.
selene_xy = selene_sites.geometry
wueller_xy = wueller_sites.geometry
sel_id_to_xy = {
    str(r["site_id"]): (geom.x, geom.y)
    for (_, r), geom in zip(selene_sites.iterrows(), selene_xy, strict=True)
}
w_id_to_xy = {
    str(r["wueller_site_id"]): (geom.x, geom.y)
    for (_, r), geom in zip(wueller_sites.iterrows(), wueller_xy, strict=True)
}
for entry in result["per_selene_site"]:
    if not entry["matched"]:
        continue
    sx, sy = sel_id_to_xy[entry["site_id"]]
    wx, wy = w_id_to_xy[entry["nearest_wueller_id"]]
    ax.plot([sx, wx], [sy, wy], color="#444444", linewidth=0.8, zorder=8, alpha=0.8)

ax.set_aspect("equal")
ax.set_xlabel("polar stereographic x (m)")
ax.set_ylabel("polar stereographic y (m)")
banner = " (SYNTHETIC PLACEHOLDER)" if USING_PLACEHOLDER else ""
ax.set_title(
    f"selene-base v1.3 vs Wueller 2026{banner} -- "
    f"{result['n_selene_matched']}/{result['n_selene_sites']} selene matched, "
    f"{result['n_wueller_matched']}/{result['n_wueller_sites']} Wueller matched "
    f"(threshold = {result['match_threshold_km']:.1f} km)"
)
legend_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor="#ffd400",
        markeredgecolor="black",
        markersize=10,
        label=f"selene v1.3 ({result['n_selene_sites']} sites)",
    ),
    Line2D(
        [0],
        [0],
        marker="^",
        color="w",
        markerfacecolor="#3b82f6",
        markeredgecolor="black",
        markersize=10,
        label=f"Wueller 2026 ({result['n_wueller_sites']} sites)"
        + (" — placeholder" if USING_PLACEHOLDER else ""),
    ),
    Line2D([0], [0], color="#444444", linewidth=0.8, label="matched pair"),
    Line2D(
        [0],
        [0],
        color="#d6232b",
        linewidth=1.5,
        label="USGS Artemis III polygon",
    ),
]
ax.legend(handles=legend_handles, loc="lower right", fontsize=9)
fig.tight_layout()
out_path = OUT_DIR / "selene_vs_wueller.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")

# %%
# Distance histogram: from each selene site to its nearest Wueller match.
distances = np.array([e["distance_km"] for e in result["per_selene_site"]])
fig, ax = plt.subplots(figsize=(7, 4))
if distances.size > 0 and np.isfinite(distances).any():
    ax.hist(
        distances[np.isfinite(distances)],
        bins=20,
        color="#3b82f6",
        edgecolor="black",
    )
    ax.axvline(
        result["match_threshold_km"],
        color="#d6232b",
        linestyle="--",
        linewidth=1.0,
        label=f"match threshold ({result['match_threshold_km']:.1f} km)",
    )
    ax.legend(loc="upper right")
ax.set_xlabel("distance from selene site to nearest Wueller site (km)")
ax.set_ylabel("number of selene sites")
ax.set_title(f"selene → Wueller nearest-distance distribution{banner}")
fig.tight_layout()
out_path = OUT_DIR / "wueller_distance_hist.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")

# %%
# Per-region match-count bar chart.
regions = [e["region"] for e in result["per_region"]]
n_selene = [e["n_selene"] for e in result["per_region"]]
n_wueller = [e["n_wueller"] for e in result["per_region"]]
n_matched = [e["n_matched"] for e in result["per_region"]]
fig, ax = plt.subplots(figsize=(11, 5))
x = np.arange(len(regions))
w = 0.27
ax.bar(x - w, n_selene, width=w, color="#ffd400", edgecolor="black", label="selene n")
ax.bar(x, n_wueller, width=w, color="#3b82f6", edgecolor="black", label="Wueller n")
ax.bar(x + w, n_matched, width=w, color="#22c55e", edgecolor="black", label="matched")
ax.set_xticks(x)
ax.set_xticklabels(regions, rotation=30, ha="right")
ax.set_ylabel("site count")
ax.set_title(f"Per-region site counts and matches{banner}")
ax.legend()
fig.tight_layout()
out_path = OUT_DIR / "wueller_per_region_bars.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")
