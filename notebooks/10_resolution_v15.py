# %% [markdown]
# # Resolution analysis (v1.5 — 20 m Wueller-class)
#
# Produces two plots:
#
# - `docs/img/selene_vs_wueller_20m.png` — the v1.4.2 headline overlay
#   regenerated from the per-region tiled (20 m) catalog.
# - `docs/img/resolution_sensitivity_v15.png` — per-region side-by-side
#   of the v1.4.2 (240 m) and v1.5 (20 m) Wueller-comparison numbers.
#
# Run after `selene preprocess --tiled-per-region --resolution 20`,
# `selene rank-per-region --tiled-per-region --resolution 20`, and
# `selene compare-wueller --sites data/outputs/per_region_tiled/sites.geojson
# --outputs-dir data/outputs/v15`.

# %%
from __future__ import annotations

import json
import os

os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
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

V14_SITES = Path("data/outputs/per_region/sites.geojson")
V15_SITES = Path("data/outputs/per_region_tiled/sites.geojson")
V14_COMPARE = Path("data/outputs/wueller_comparison.json")
V15_COMPARE = Path("data/outputs/v15/wueller_comparison.json")

LABEL_OFFSETS_M: dict[str, tuple[float, float]] = {
    "Malapert Massif": (+27_000, 0),
    "Mons Mouton": (0, -22_000),
    "Mons Mouton Plateau": (0, +25_000),
    "Haworth": (-30_000, 0),
    "Nobile Rim 1": (+22_000, -18_000),
    "Nobile Rim 2": (+8_000, +18_000),
    "Peak Near Cabeus B": (-15_000, -12_000),
    "Slater Plain": (+55_000, -22_000),
    "de Gerlache Rim 2": (-110_000, +5_000),
}
PLOT_XMIN, PLOT_XMAX = -200_000, +180_000
PLOT_YMIN, PLOT_YMAX = -85_000, +200_000


# %% [markdown]
# ## Plot 1 — selene-base v1.5 (20 m) vs Wueller 2026

# %%
selene_polar = gpd.read_file(V15_SITES).to_crs(POLAR_PROJ)
selene_geog = selene_polar.to_crs(GEOG_PROJ)
wueller_all = load_wueller_sites(target_crs=POLAR_PROJ)
wueller_in_scope_polar = wueller_all[wueller_all["in_usgs_scope"]].reset_index(drop=True)
wueller_in_scope_geog = wueller_in_scope_polar.to_crs(GEOG_PROJ)
polygons = regions_polygons_to_geodataframe(target_crs=POLAR_PROJ)

result = compare_sites(
    selene_geog,
    wueller_in_scope_geog,
    match_threshold_km=DEFAULT_MATCH_THRESHOLD_KM,
    filter_to_usgs_scope=False,
)
print(
    f"v1.5: selene={result['n_selene_sites']}, wueller_in_scope={result['n_wueller_sites']}, "
    f"matched={result['n_selene_matched']}, "
    f"median_km={result['median_match_distance_km']:.2f}"
)

matched_ids = {e["site_id"] for e in result["per_selene_site"] if e["matched"]}
sel_xy = {
    str(r["site_id"]): (g.x, g.y)
    for (_, r), g in zip(selene_polar.iterrows(), selene_polar.geometry, strict=True)
}
w_xy = {
    str(r["wueller_site_id"]): (g.x, g.y)
    for (_, r), g in zip(
        wueller_in_scope_polar.iterrows(), wueller_in_scope_polar.geometry, strict=True
    )
}

# %%
plt.style.use("default")
fig, ax = plt.subplots(figsize=(13, 10))
fig.patch.set_facecolor("#1a1a1a")
ax.set_facecolor("#1a1a1a")

for r_km in (50, 100, 150):
    ax.add_patch(
        plt.Circle(
            (0, 0),
            r_km * 1000,
            fill=False,
            edgecolor="#3a3a3a",
            linewidth=0.6,
            linestyle="--",
            zorder=1,
        )
    )
ax.axhline(0, color="#3a3a3a", linewidth=0.5, zorder=1)
ax.axvline(0, color="#3a3a3a", linewidth=0.5, zorder=1)

polygons.boundary.plot(ax=ax, color="#e63946", linewidth=1.6, zorder=4)
for _, row in polygons.iterrows():
    name = str(row["Region"])
    cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
    dx, dy = LABEL_OFFSETS_M.get(name, (15_000, 15_000))
    ax.annotate(
        name,
        xy=(cx + dx, cy + dy),
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
            (cx + dx, cy + dy),
            arrowstyle="-",
            color="#7a2727",
            linewidth=0.6,
            zorder=3,
        )
    )

ax.scatter(
    wueller_in_scope_polar.geometry.x,
    wueller_in_scope_polar.geometry.y,
    s=28,
    facecolor="#facc15",
    edgecolor="#78350f",
    linewidth=0.5,
    zorder=6,
)

for entry in result["per_selene_site"]:
    if not entry["matched"]:
        continue
    sx, sy = sel_xy[str(entry["site_id"])]
    wx, wy = w_xy[str(entry["nearest_wueller_id"])]
    ax.plot([sx, wx], [sy, wy], color="#cbd5e1", linewidth=1.4, zorder=7, alpha=0.95)

matched_mask = selene_polar["site_id"].astype(str).isin(matched_ids)
matched_sel = selene_polar[matched_mask]
unmatched_sel = selene_polar[~matched_mask]
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
for _, row in selene_polar.iterrows():
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

ax.scatter(
    [0], [0], s=180, marker="*",
    facecolor="white", edgecolor="black", linewidth=0.8, zorder=12,
)
ax.annotate(
    "South Pole", xy=(0, 0), xytext=(8, -10),
    textcoords="offset points", fontsize=8, color="white", zorder=12,
)

ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
ax.set_ylim(PLOT_YMIN, PLOT_YMAX)
ax.set_aspect("equal")

x_ticks = [-200_000, -100_000, 0, 100_000]
y_ticks = [-50_000, 0, 50_000, 100_000, 150_000, 200_000]
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xticklabels([f"{t // 1000:+d} km" for t in x_ticks], color="#d1d5db")
ax.set_yticklabels([f"{t // 1000:+d} km" for t in y_ticks], color="#d1d5db")
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
    f"selene-base v1.5 (20 m, per-region tiled): {n_total} sites; "
    f"{pct}% match within 5 km of Wueller 2026 (median {median_km:.2f} km)",
    fontsize=20,
    fontweight="bold",
    color="#ffffff",
    y=0.985,
)
ax.set_title(
    "Per-region tiled HLS filtering at 20 m. "
    "Reference: Wueller et al. 2026 (JGR Planets, doi:10.1029/2025JE009434).",
    fontsize=11,
    color="#e5e7eb",
    pad=14,
)

legend = ax.legend(
    handles=[
        Line2D([0], [0], color="#e63946", linewidth=1.6,
               label="NASA Artemis III regions (USGS polygons)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#facc15",
               markeredgecolor="#78350f", markersize=7, linestyle="None",
               label=f"Wueller 2026 in-scope sites (n={n_wueller})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#06d6a0",
               markeredgecolor="white", markersize=10, linestyle="None",
               label=f"selene-base matched ({n_matched})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
               markeredgecolor="#06d6a0", markeredgewidth=1.6, markersize=10,
               linestyle="None",
               label=f"selene-base unmatched ({n_total - n_matched})"),
        Line2D([0], [0], color="#9ca3af", linewidth=0.9,
               label="matched-pair distance ≤ 5 km"),
    ],
    loc="lower left",
    fontsize=9,
    framealpha=0.92,
    facecolor="#262626",
    edgecolor="#3a3a3a",
    labelcolor="white",
)
for text in legend.get_texts():
    text.set_color("white")

fig.tight_layout(rect=(0, 0, 1, 0.92))
out_path = OUT_DIR / "selene_vs_wueller_20m.png"
fig.savefig(out_path, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print(f"wrote {out_path}")


# %% [markdown]
# ## Plot 2 — resolution sensitivity v1.4.2 (240 m) vs v1.5 (20 m)

# %%
v14 = json.loads(V14_COMPARE.read_text())
v15 = json.loads(V15_COMPARE.read_text())

regions = sorted({p["region"] for p in v14["per_region"] + v15["per_region"]})
v14_by_region = {p["region"]: p for p in v14["per_region"]}
v15_by_region = {p["region"]: p for p in v15["per_region"]}


def _extract(per: dict | None, key: str, default: float = float("nan")) -> float:
    if per is None or per.get(key) is None:
        return default
    return float(per[key])


v14_match_pct = []
v15_match_pct = []
v14_med_km = []
v15_med_km = []
for r in regions:
    p14 = v14_by_region.get(r)
    p15 = v15_by_region.get(r)
    n14 = _extract(p14, "n_selene", 0.0)
    m14 = _extract(p14, "n_matched", 0.0)
    n15 = _extract(p15, "n_selene", 0.0)
    m15 = _extract(p15, "n_matched", 0.0)
    v14_match_pct.append(100.0 * m14 / n14 if n14 > 0 else 0.0)
    v15_match_pct.append(100.0 * m15 / n15 if n15 > 0 else 0.0)
    v14_med_km.append(_extract(p14, "median_distance_km"))
    v15_med_km.append(_extract(p15, "median_distance_km"))

fig, (ax_pct, ax_dist) = plt.subplots(1, 2, figsize=(15, 6))
fig.patch.set_facecolor("#1a1a1a")
for ax in (ax_pct, ax_dist):
    ax.set_facecolor("#1a1a1a")
    for spine in ax.spines.values():
        spine.set_color("#3a3a3a")
    ax.tick_params(colors="#d1d5db")

x = np.arange(len(regions))
width = 0.36

ax_pct.bar(x - width / 2, v14_match_pct, width,
           color="#94a3b8", edgecolor="#1a1a1a", label="v1.4.2 (240 m)")
ax_pct.bar(x + width / 2, v15_match_pct, width,
           color="#06d6a0", edgecolor="#1a1a1a", label="v1.5 (20 m)")
ax_pct.set_xticks(x)
ax_pct.set_xticklabels(
    [r.replace(" ", "\n") for r in regions], fontsize=8.5, color="#d1d5db",
)
ax_pct.set_ylabel("selene matched (%)", color="#d1d5db")
ax_pct.set_title(
    "Per-region agreement: % of selene sites matched within 5 km of Wueller",
    color="#ffffff", fontsize=11, pad=10,
)
ax_pct.set_ylim(0, 105)
leg = ax_pct.legend(facecolor="#262626", edgecolor="#3a3a3a", labelcolor="white")
for t in leg.get_texts():
    t.set_color("white")
ax_pct.grid(axis="y", color="#3a3a3a", linewidth=0.5, alpha=0.5)

# Right panel: median matched-pair distance
v14_dist_plot = [d if d == d else 0.0 for d in v14_med_km]
v15_dist_plot = [d if d == d else 0.0 for d in v15_med_km]
ax_dist.bar(x - width / 2, v14_dist_plot, width,
            color="#94a3b8", edgecolor="#1a1a1a", label="v1.4.2 (240 m)")
ax_dist.bar(x + width / 2, v15_dist_plot, width,
            color="#06d6a0", edgecolor="#1a1a1a", label="v1.5 (20 m)")
ax_dist.set_xticks(x)
ax_dist.set_xticklabels(
    [r.replace(" ", "\n") for r in regions], fontsize=8.5, color="#d1d5db",
)
ax_dist.set_ylabel("median matched-pair distance (km)", color="#d1d5db")
ax_dist.set_title(
    "Per-region median matched-pair distance",
    color="#ffffff", fontsize=11, pad=10,
)
leg = ax_dist.legend(facecolor="#262626", edgecolor="#3a3a3a", labelcolor="white")
for t in leg.get_texts():
    t.set_color("white")
ax_dist.grid(axis="y", color="#3a3a3a", linewidth=0.5, alpha=0.5)

# Headline numbers across versions
hl = (
    f"v1.4.2: {v14['n_selene_matched']}/{v14['n_selene_sites']} matched "
    f"(median {v14['median_match_distance_km']:.2f} km)   →   "
    f"v1.5: {v15['n_selene_matched']}/{v15['n_selene_sites']} matched "
    f"(median {v15['median_match_distance_km']:.2f} km)"
)
fig.suptitle(
    "Resolution sensitivity: 240 m → 20 m methodology converges",
    color="#ffffff", fontsize=18, fontweight="bold", y=0.995,
)
fig.text(0.5, 0.93, hl, ha="center", color="#e5e7eb", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.91))
out_path = OUT_DIR / "resolution_sensitivity_v15.png"
fig.savefig(out_path, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print(f"wrote {out_path}")
