# %% [markdown]
# # v1.5 catalog report — plots
#
# Generates the three plots embedded in `docs/v1.5_catalog_report.md`:
#
# - `docs/img/v15_catalog_score_distributions.png` — per-criterion
#   score distribution, selene's 69 sites vs Wueller's 73 in-scope sites
#   (boxplots, both evaluated under selene's 240 m criterion rasters).
# - `docs/img/v15_catalog_aggregate_histogram.png` — aggregate score
#   histogram, selene 69 vs Wueller 73 overlaid.
# - `docs/img/v15_catalog_region_heatmap.png` — per-USGS-region
#   agreement matrix: matched %, mean selene aggregate, mean Wueller-side
#   aggregate (under selene's criteria), median pair distance.
#
# Run after `selene rank-per-region --tiled-per-region --resolution 20`,
# `selene compare-wueller --sites data/outputs/per_region_tiled/sites.geojson
# --outputs-dir data/outputs/v15`, and `selene score-wueller-sites`.

# %%
from __future__ import annotations

import json
import os

os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = Path("docs/img")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SELENE_MAIN = Path("data/outputs/v1.5_catalog_main_table.csv")
WUELLER_EVAL = Path("data/outputs/v1.5_catalog_wueller_evaluation.csv")
V15_COMPARE = Path("data/outputs/v15/wueller_comparison.json")

CRITERIA = (
    "slope",
    "illumination",
    "coupling",
    "thermal",
    "ice",
    "los_to_earth",
)
PRETTY = {
    "slope": "Slope",
    "illumination": "Illumination",
    "coupling": "Coupling",
    "thermal": "Thermal",
    "ice": "Ice",
    "los_to_earth": "LOS-to-Earth",
}

selene = pd.read_csv(SELENE_MAIN)
wueller = pd.read_csv(WUELLER_EVAL)
v15 = json.loads(V15_COMPARE.read_text())


# %% [markdown]
# ## Plot 1 — per-criterion score distribution

# %%
fig, axes = plt.subplots(1, len(CRITERIA), figsize=(15, 5), sharey=True)
fig.patch.set_facecolor("#1a1a1a")
for ax, crit in zip(axes, CRITERIA, strict=True):
    ax.set_facecolor("#1a1a1a")
    sel_vals = selene[f"score_{crit}"].dropna().to_numpy()
    wue_vals = wueller[f"score_{crit}"].dropna().to_numpy()
    bp = ax.boxplot(
        [sel_vals, wue_vals],
        positions=[0, 1],
        widths=0.55,
        patch_artist=True,
        showfliers=True,
        medianprops=dict(color="white", linewidth=1.5),
        flierprops=dict(marker="o", markersize=3, markerfacecolor="#888", markeredgecolor="none"),
    )
    for patch, color in zip(bp["boxes"], ["#06d6a0", "#facc15"], strict=True):
        patch.set_facecolor(color)
        patch.set_edgecolor("#1a1a1a")
        patch.set_alpha(0.85)
    for whisker in bp["whiskers"]:
        whisker.set_color("#9ca3af")
    for cap in bp["caps"]:
        cap.set_color("#9ca3af")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["selene\n(n=69)", "Wueller\n(n=73)"], color="#d1d5db", fontsize=9)
    ax.set_title(PRETTY[crit], color="#ffffff", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(colors="#d1d5db")
    for spine in ax.spines.values():
        spine.set_color("#3a3a3a")
    ax.grid(axis="y", color="#3a3a3a", linewidth=0.4, alpha=0.4)
axes[0].set_ylabel("score", color="#d1d5db")
fig.suptitle(
    "Per-criterion score distribution (240 m): selene's 69 sites vs Wueller's 73 in-scope sites",
    color="#ffffff",
    fontsize=15,
    fontweight="bold",
    y=0.995,
)
fig.text(
    0.5,
    0.93,
    "Both populations evaluated against selene's 240 m criterion rasters. Higher = better.",
    ha="center",
    color="#e5e7eb",
    fontsize=10,
)
fig.tight_layout(rect=(0, 0, 1, 0.91))
out_path = OUT_DIR / "v15_catalog_score_distributions.png"
fig.savefig(out_path, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print(f"wrote {out_path}")


# %% [markdown]
# ## Plot 2 — aggregate score histogram

# %%
fig, ax = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor("#1a1a1a")
ax.set_facecolor("#1a1a1a")
bins = np.linspace(0.0, 1.0, 26)
ax.hist(
    selene["score"].to_numpy(),
    bins=bins,
    alpha=0.75,
    color="#06d6a0",
    edgecolor="#1a1a1a",
    label=f"selene 69 sites (median {selene['score'].median():.3f})",
)
ax.hist(
    wueller["aggregate_score"].dropna().to_numpy(),
    bins=bins,
    alpha=0.65,
    color="#facc15",
    edgecolor="#1a1a1a",
    label=f"Wueller 73 in-scope (median {wueller['aggregate_score'].median():.3f})",
)
ax.axvline(
    selene["score"].median(),
    color="#06d6a0",
    linewidth=1.4,
    linestyle="--",
    alpha=0.9,
)
ax.axvline(
    wueller["aggregate_score"].median(),
    color="#facc15",
    linewidth=1.4,
    linestyle="--",
    alpha=0.9,
)
ax.set_xlabel("aggregate score (240 m, weighted sum across active criteria)", color="#d1d5db")
ax.set_ylabel("number of sites", color="#d1d5db")
ax.tick_params(colors="#d1d5db")
for spine in ax.spines.values():
    spine.set_color("#3a3a3a")
ax.grid(axis="y", color="#3a3a3a", linewidth=0.4, alpha=0.4)
ax.set_xlim(0.0, 1.0)
leg = ax.legend(facecolor="#262626", edgecolor="#3a3a3a", labelcolor="white", loc="upper left")
for t in leg.get_texts():
    t.set_color("white")
fig.suptitle(
    "Aggregate score: selene's catalog vs Wueller's catalog at Wueller's coordinates",
    color="#ffffff",
    fontsize=15,
    fontweight="bold",
    y=0.985,
)
ax.set_title(
    "Both evaluated against selene's 240 m aggregate score raster.",
    color="#e5e7eb",
    fontsize=10,
    pad=10,
)
fig.tight_layout(rect=(0, 0, 1, 0.92))
out_path = OUT_DIR / "v15_catalog_aggregate_histogram.png"
fig.savefig(out_path, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print(f"wrote {out_path}")


# %% [markdown]
# ## Plot 3 — region-level agreement heatmap

# %%
per_region = {p["region"]: p for p in v15["per_region"]}
all_regions = sorted(set(per_region) | set(selene["region_name"].unique()))

rows = []
for region in all_regions:
    sel_in_region = selene[selene["region_name"] == region]
    wue_in_region = wueller[wueller["region"] == region]
    pr = per_region.get(region)
    n_sel = int(len(sel_in_region))
    n_wue = int(len(wue_in_region))
    matched = int(pr["n_matched"]) if pr else 0
    match_pct = 100.0 * matched / n_sel if n_sel else 0.0
    median_dist = float(pr["median_distance_km"]) if pr and pr.get("median_distance_km") else np.nan
    sel_mean_agg = float(sel_in_region["score"].mean()) if n_sel else np.nan
    wue_mean_agg = float(wue_in_region["aggregate_score"].mean()) if n_wue else np.nan
    rows.append(
        {
            "region": region,
            "match_pct": match_pct,
            "selene_mean_agg": sel_mean_agg,
            "wueller_mean_agg": wue_mean_agg,
            "median_dist_km": median_dist,
        }
    )
heat = pd.DataFrame(rows).set_index("region")

fig, ax = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor("#1a1a1a")
ax.set_facecolor("#1a1a1a")

cell_grid = heat[["match_pct", "selene_mean_agg", "wueller_mean_agg", "median_dist_km"]].to_numpy()
# Each column has its own scale, so show colour per column normalised within column.
norm_cells = np.zeros_like(cell_grid)
for j in range(cell_grid.shape[1]):
    col = cell_grid[:, j]
    finite = np.isfinite(col)
    if finite.sum() == 0:
        norm_cells[:, j] = 0.0
        continue
    cmin = np.nanmin(col)
    cmax = np.nanmax(col)
    if cmax == cmin:
        norm_cells[:, j] = 0.5
    else:
        norm_cells[:, j] = (col - cmin) / (cmax - cmin)

# match_pct, selene_mean_agg, wueller_mean_agg: higher = better (green high)
# median_dist_km: lower = better (invert)
norm_cells[:, 3] = 1.0 - norm_cells[:, 3]
norm_cells = np.where(np.isnan(cell_grid), np.nan, norm_cells)

ax.imshow(norm_cells, cmap="viridis", aspect="auto", vmin=0.0, vmax=1.0)
ax.set_yticks(range(len(heat.index)))
ax.set_yticklabels(heat.index, color="#d1d5db", fontsize=10)
ax.set_xticks(range(4))
ax.set_xticklabels(
    ["matched %", "selene mean agg", "Wueller mean agg\n(at their coords)", "median dist (km)"],
    color="#d1d5db",
    fontsize=10,
)
ax.tick_params(colors="#d1d5db")
for spine in ax.spines.values():
    spine.set_color("#3a3a3a")

# Annotate
for i in range(cell_grid.shape[0]):
    for j in range(cell_grid.shape[1]):
        v = cell_grid[i, j]
        if not np.isfinite(v):
            txt = "—"
        elif j == 0:
            txt = f"{v:.0f} %"
        elif j == 3:
            txt = f"{v:.2f}"
        else:
            txt = f"{v:.3f}"
        nv = norm_cells[i, j]
        color = "#1a1a1a" if (np.isfinite(nv) and nv > 0.55) else "#ffffff"
        ax.text(j, i, txt, ha="center", va="center", color=color, fontsize=10, fontweight="bold")

fig.suptitle(
    "Per-region agreement: spatial match + per-criterion methodology agreement",
    color="#ffffff",
    fontsize=15,
    fontweight="bold",
    y=0.99,
)
ax.set_title(
    "Brighter = better. Selene/Wueller mean aggregate computed against "
    "selene's 240 m score raster.",
    color="#e5e7eb",
    fontsize=10,
    pad=10,
)
fig.tight_layout(rect=(0, 0, 1, 0.93))
out_path = OUT_DIR / "v15_catalog_region_heatmap.png"
fig.savefig(out_path, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print(f"wrote {out_path}")
