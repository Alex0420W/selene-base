# %% [markdown]
# # Week 3 — full pipeline preview
#
# Loads the cached per-criterion score COGs that `selene score` writes,
# the aggregate score COG, and the ranked top sites. Saves three figures
# to `data/outputs/sanity/`:
#
# 1. `criteria_grid.png` — every available per-criterion score map.
# 2. `aggregate.png` — the final weighted-sum score map.
# 3. `top_sites.png` — top sites overlaid on the aggregate.
#
# Run after:
# ```
# selene download lola && selene download illumination && selene download robbins
# selene preprocess
# selene score
# selene rank --top-n 20
# ```

# %%
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

from selene_base.scoring.ranking import DEFAULT_CRITERIA

PROCESSED = Path("data/processed")
OUTPUTS = Path("data/outputs")
SANITY = OUTPUTS / "sanity"
SANITY.mkdir(parents=True, exist_ok=True)


def _open(path: Path) -> xr.DataArray:
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


def _coarsen(da: xr.DataArray, factor: int) -> np.ndarray:
    return da.coarsen(y=factor, x=factor, boundary="trim").mean().to_numpy()


# %% [markdown]
# ## 1. Per-criterion score maps
#
# We render a panel for each criterion whose score COG is on disk, and
# leave the rest blank with a "data not yet downloaded" annotation. That
# keeps the figure layout stable as week 4 fills in the remaining
# criteria.

# %%
COARSE = 6
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for ax, crit in zip(axes.flat, DEFAULT_CRITERIA, strict=True):
    cog = PROCESSED / "scored" / f"{crit}_score_southpole_240m.tif"
    if not cog.exists():
        ax.text(
            0.5,
            0.5,
            f"{crit}\n(score map not present)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="grey",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(crit)
        continue
    da = _open(cog)
    im = ax.imshow(_coarsen(da, COARSE), cmap="viridis", origin="upper", vmin=0, vmax=1)
    ax.set_title(f"{crit} score")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, label="score")

fig.suptitle("Per-criterion score maps (week 3)")
fig.tight_layout()
out = SANITY / "criteria_grid.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"[saved] {out}")


# %% [markdown]
# ## 2. Aggregate score map
# Weighted sum across the criteria that ran (renormalised across the present set).

# %%
agg_path = OUTPUTS / "score_southpole.tif"
if not agg_path.exists():
    raise SystemExit(f"missing {agg_path}; run `selene score`")
agg = _open(agg_path)
arr = _coarsen(agg, COARSE)

fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(arr, cmap="inferno", origin="upper", vmin=0, vmax=1)
ax.set_title("Aggregate suitability score (south pole, 240 m)")
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im, ax=ax, label="score")
fig.tight_layout()
out = SANITY / "aggregate.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"[saved] {out}")


# %% [markdown]
# ## 3. Top sites overlaid on the aggregate
# Markers in the projected grid; sites pulled straight from
# `data/outputs/top_sites.geojson`.

# %%
sites_path = OUTPUTS / "top_sites.geojson"
if not sites_path.exists():
    print(f"[skip] {sites_path}; run `selene rank`")
else:
    sites = gpd.read_file(sites_path)
    print(f"loaded {len(sites)} site(s)")

    transform = agg.rio.transform()
    pixel_size = float(abs(transform.a))
    height, width = agg.sizes["y"], agg.sizes["x"]

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(arr, cmap="inferno", origin="upper", vmin=0, vmax=1)
    # Convert site x_m/y_m back to coarsened-pixel coordinates.
    cx = (sites["x_m"].to_numpy() - transform.c) / transform.a / COARSE
    cy = (sites["y_m"].to_numpy() - transform.f) / transform.e / COARSE
    ax.scatter(cx, cy, s=80, edgecolor="white", facecolor="cyan", linewidths=1.5)
    for _, row in sites.iterrows():
        rx = (row["x_m"] - transform.c) / transform.a / COARSE
        ry = (row["y_m"] - transform.f) / transform.e / COARSE
        ax.annotate(
            str(int(row["rank"])),
            (rx, ry),
            xytext=(4, 4),
            textcoords="offset points",
            color="white",
            fontsize=8,
            fontweight="bold",
        )
    ax.set_title(f"Top {len(sites)} candidate sites")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    out = SANITY / "top_sites.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[saved] {out}")
