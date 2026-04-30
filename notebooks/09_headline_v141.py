# %% [markdown]
# # v1.4.1 headline image (combined per-region catalog + Wueller agreement)
#
# Single-purpose notebook: produces `docs/img/headline_v141.png`, the
# README's top-of-page visual. Combines the v1.3 per-region selene
# catalog (23 sites across 8/9 USGS regions), the in-scope Wueller
# 2026 sites (73 of 130), the USGS polygons, and match lines for the
# 18 selene sites within 5 km of a Wueller site.

# %%
from __future__ import annotations

import os

# Same lunar/Earth ellipsoid override as the week 10–12 notebooks.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe
from selene_base.validation.wueller_comparison import (
    DEFAULT_MATCH_THRESHOLD_KM,
    compare_sites,
    load_wueller_sites,
)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)
GEOG_PROJ = "+proj=longlat +R=1737400 +no_defs +type=crs"

OUT_DIR = Path("docs/img")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SELENE_SITES_GEOJSON = Path("data/outputs/per_region/sites.geojson")

# Plot bounds: ±200 km, slightly tighter than the ±304 km full grid
# since every candidate site clusters near the pole.
EXTENT_KM = 200

# Manual label offsets (metres) for the central cluster, where four
# polygons overlap visually. Anchored at the polygon centroid; arrow
# drawn from the centroid to the offset label position. Offsets aim
# the labels into the empty quadrants of the +200 km extent.
LABEL_OFFSETS_M: dict[str, tuple[float, float]] = {
    # Top cluster (+75 to +90 km x, +100 to +150 km y) — fan labels out.
    "Malapert Massif": (-50_000, +60_000),
    "Mons Mouton": (-30_000, +60_000),
    "Mons Mouton Plateau": (+85_000, +35_000),
    "Nobile Rim 1": (+90_000, -10_000),
    "Nobile Rim 2": (+30_000, -25_000),
    # The four outliers — fan into empty space.
    "Haworth": (-90_000, +25_000),
    "de Gerlache Rim 2": (-110_000, +5_000),
    "Peak Near Cabeus B": (+10_000, +30_000),
    "Slater Plain": (+60_000, -25_000),
}

# %%
selene_sites_polar = gpd.read_file(SELENE_SITES_GEOJSON).to_crs(POLAR_PROJ)
selene_sites_geog = selene_sites_polar.to_crs(GEOG_PROJ)
wueller_all = load_wueller_sites(target_crs=POLAR_PROJ)
wueller_in_scope_polar = wueller_all[wueller_all["in_usgs_scope"]].reset_index(drop=True)
wueller_in_scope_geog = wueller_in_scope_polar.to_crs(GEOG_PROJ)
polygons = regions_polygons_to_geodataframe(target_crs=POLAR_PROJ)

result = compare_sites(
    selene_sites_geog,
    wueller_in_scope_geog,
    match_threshold_km=DEFAULT_MATCH_THRESHOLD_KM,
    filter_to_usgs_scope=False,  # already pre-filtered above
)
print(
    f"selene={result['n_selene_sites']}, wueller_in_scope={result['n_wueller_sites']}, "
    f"matched={result['n_selene_matched']}, "
    f"median_km={result['median_match_distance_km']:.2f}, "
    f"placeholder={result['using_synthetic_placeholder']}"
)

matched_selene_ids = {e["site_id"] for e in result["per_selene_site"] if e["matched"]}
selene_to_nearest_wueller_id = {
    e["site_id"]: e["nearest_wueller_id"] for e in result["per_selene_site"] if e["matched"]
}

sel_id_to_xy = {
    str(r["site_id"]): (geom.x, geom.y)
    for (_, r), geom in zip(selene_sites_polar.iterrows(), selene_sites_polar.geometry, strict=True)
}
w_id_to_xy = {
    str(r["wueller_site_id"]): (geom.x, geom.y)
    for (_, r), geom in zip(
        wueller_in_scope_polar.iterrows(), wueller_in_scope_polar.geometry, strict=True
    )
}

# %%
plt.style.use("default")
fig, ax = plt.subplots(figsize=(13, 11))
fig.patch.set_facecolor("#1a1a1a")
ax.set_facecolor("#1a1a1a")

# Subtle radial graticule: range rings every 50 km out to 200 km, plus
# axes through the pole.
for r_km in (50, 100, 150, 200):
    circle = plt.Circle(
        (0, 0),
        r_km * 1000,
        fill=False,
        edgecolor="#3a3a3a",
        linewidth=0.6,
        linestyle="--",
        zorder=1,
    )
    ax.add_patch(circle)
ax.axhline(0, color="#3a3a3a", linewidth=0.5, zorder=1)
ax.axvline(0, color="#3a3a3a", linewidth=0.5, zorder=1)

# Layer 2: USGS polygons (red outline, no fill).
polygons.boundary.plot(ax=ax, color="#e63946", linewidth=1.6, zorder=4)

# Region labels with manual offset + connector lines.
for _, row in polygons.iterrows():
    name = str(row["Region"])
    cx = row.geometry.centroid.x
    cy = row.geometry.centroid.y
    dx, dy = LABEL_OFFSETS_M.get(name, (15_000, 15_000))
    label_x = cx + dx
    label_y = cy + dy
    ax.annotate(
        name,
        xy=(label_x, label_y),
        ha="center",
        va="center",
        fontsize=8.5,
        color="#fca5a5",
        fontweight="bold",
        zorder=11,
    )
    ax.add_patch(
        FancyArrowPatch(
            (cx, cy),
            (label_x, label_y),
            arrowstyle="-",
            color="#7a2727",
            linewidth=0.6,
            zorder=3,
        )
    )

# Layer 3: Wueller in-scope sites (yellow circles, small).
ax.scatter(
    wueller_in_scope_polar.geometry.x,
    wueller_in_scope_polar.geometry.y,
    s=28,
    facecolor="#facc15",
    edgecolor="#78350f",
    linewidth=0.5,
    zorder=6,
)

# Layer 5 (drawn before selene sites so dots sit on top of lines):
# match lines connecting each matched selene site to its nearest in-scope
# Wueller site.
for entry in result["per_selene_site"]:
    if not entry["matched"]:
        continue
    sx, sy = sel_id_to_xy[str(entry["site_id"])]
    wx, wy = w_id_to_xy[str(entry["nearest_wueller_id"])]
    ax.plot([sx, wx], [sy, wy], color="#9ca3af", linewidth=0.9, zorder=7, alpha=0.85)

# Layer 4: selene sites — solid cyan for matched, outlined for unmatched.
matched_mask = selene_sites_polar["site_id"].astype(str).isin(matched_selene_ids)
matched_sel = selene_sites_polar[matched_mask]
unmatched_sel = selene_sites_polar[~matched_mask]
ax.scatter(
    matched_sel.geometry.x,
    matched_sel.geometry.y,
    s=70,
    facecolor="#06d6a0",
    edgecolor="white",
    linewidth=1.2,
    zorder=10,
)
ax.scatter(
    unmatched_sel.geometry.x,
    unmatched_sel.geometry.y,
    s=70,
    facecolor="none",
    edgecolor="#06d6a0",
    linewidth=1.6,
    zorder=10,
)
for _, row in selene_sites_polar.iterrows():
    ax.annotate(
        str(int(row["site_id"])),
        (row.geometry.x, row.geometry.y),
        xytext=(5, 4),
        textcoords="offset points",
        color="white",
        fontsize=7.5,
        fontweight="bold",
        zorder=11,
    )

# Layer 1: south pole marker, drawn last so it sits visually on top.
ax.scatter(
    [0],
    [0],
    s=180,
    marker="*",
    facecolor="white",
    edgecolor="black",
    linewidth=0.8,
    zorder=12,
)
ax.annotate(
    "South Pole",
    xy=(0, 0),
    xytext=(8, -10),
    textcoords="offset points",
    fontsize=8,
    color="white",
    zorder=12,
)

# Bounds: ±200 km.
ax.set_xlim(-EXTENT_KM * 1000, EXTENT_KM * 1000)
ax.set_ylim(-EXTENT_KM * 1000, EXTENT_KM * 1000)
ax.set_aspect("equal")

# Tick labels in km for readability.
tick_m = [-200_000, -100_000, 0, 100_000, 200_000]
ax.set_xticks(tick_m)
ax.set_yticks(tick_m)
ax.set_xticklabels([f"{t // 1000:+d} km" for t in tick_m], color="#d1d5db")
ax.set_yticklabels([f"{t // 1000:+d} km" for t in tick_m], color="#d1d5db")
ax.tick_params(colors="#d1d5db")
for spine in ax.spines.values():
    spine.set_color("#3a3a3a")

ax.set_xlabel("polar stereographic x", color="#d1d5db")
ax.set_ylabel("polar stereographic y", color="#d1d5db")

n_matched = result["n_selene_matched"]
n_total = result["n_selene_sites"]
n_wueller = result["n_wueller_sites"]
median_km = result["median_match_distance_km"]
pct = round(100 * n_matched / n_total)

fig.suptitle(
    f"selene-base v1.4.1: {n_total} sites across 8/9 USGS regions; "
    f"{pct}% match within 5 km of Wueller 2026 (median {median_km:.2f} km)",
    fontsize=15,
    color="white",
    y=0.975,
)
ax.set_title(
    "Comparison against Wueller et al. 2026 (JGR Planets, doi:10.1029/2025JE009434), "
    "CC-BY 4.0 data deposit",
    fontsize=10,
    color="#d1d5db",
    pad=14,
)

# Legend, top-right.
legend_handles = [
    Line2D(
        [0],
        [0],
        marker="None",
        color="#e63946",
        linewidth=1.6,
        label="NASA Artemis III regions (USGS polygons)",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor="#facc15",
        markeredgecolor="#78350f",
        markersize=7,
        linestyle="None",
        label=f"Wueller 2026 in-scope sites (n={n_wueller})",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor="#06d6a0",
        markeredgecolor="white",
        markersize=10,
        linestyle="None",
        label=f"selene-base matched ({n_matched})",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor="none",
        markeredgecolor="#06d6a0",
        markeredgewidth=1.6,
        markersize=10,
        linestyle="None",
        label=f"selene-base unmatched ({n_total - n_matched})",
    ),
    Line2D(
        [0],
        [0],
        color="#9ca3af",
        linewidth=0.9,
        label="matched-pair distance ≤ 5 km",
    ),
]
legend = ax.legend(
    handles=legend_handles,
    loc="lower left",
    fontsize=9,
    framealpha=0.92,
    facecolor="#262626",
    edgecolor="#3a3a3a",
    labelcolor="white",
)
for text in legend.get_texts():
    text.set_color("white")

fig.tight_layout(rect=(0, 0, 1, 0.94))
out_path = OUT_DIR / "headline_v141.png"
fig.savefig(out_path, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print(f"saved {out_path}")
