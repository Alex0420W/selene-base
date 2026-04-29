# %% [markdown]
# # Week 7 — coupling visualisation
#
# Renders the inputs and output of the spatial-coupling criterion on
# real data, plus the new aggregate score and top-20 sites overlaid
# with NASA Artemis III candidate regions.
#
# Saves three figures into `data/outputs/sanity/`:
#
# 1. `coupling_inputs.png` — distance-to-PSR + distance-to-sunlit-ridge.
# 2. `coupling_score.png` — the product map, log-scaled so the rim
#    band shows up.
# 3. `coupling_overlay.png` — coupling score with NASA centroids and
#    our top-20 ranked sites overlaid.
#
# Run after:
# ```
# selene preprocess && selene score && selene rank
# ```

# %%
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401
from matplotlib.colors import LogNorm
from pyproj import Transformer

from selene_base.criteria import coupling as coupling_criterion
from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.nasa_regions import regions_to_geodataframe

PROCESSED = Path("data/processed")
OUTPUTS = Path("data/outputs")
SANITY = OUTPUTS / "sanity"
SANITY.mkdir(parents=True, exist_ok=True)


def _open(path: Path):
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


def _coarsen(da, factor: int = 4):
    return da.coarsen(y=factor, x=factor, boundary="trim").mean().to_numpy()


# %% [markdown]
# ## 1. Coupling inputs
# Per-pixel distance to the nearest PSR (left) and to the nearest
# "sunlit ridge" — high-illumination cell with slope in the 5°–25° band
# (right). The rim band lights up where both distances are small.

# %%
illum = _open(PROCESSED / "illumination_southpole_240m.tif")
slope_deg = _open(PROCESSED / "lola_slope_deg_southpole_240m.tif")

distance_psr = coupling_criterion.derive_distance_to_psr(illum, pixel_size_m=240.0)
distance_ridge = coupling_criterion.derive_distance_to_sunlit_ridge(
    illum, slope_deg, pixel_size_m=240.0
)
print(f"distance_to_psr   median {np.nanmedian(distance_psr.to_numpy()) / 1000:.1f} km")
print(f"distance_to_ridge median {np.nanmedian(distance_ridge.to_numpy()) / 1000:.1f} km")

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, da, title in zip(
    axes,
    (distance_psr, distance_ridge),
    ("distance to nearest PSR (km)", "distance to nearest sunlit ridge (km)"),
    strict=True,
):
    arr = _coarsen(da) / 1000.0
    im = ax.imshow(arr, cmap="magma_r", origin="upper", vmin=0, vmax=50)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, label="km")
fig.tight_layout()
out = SANITY / "coupling_inputs.png"
fig.savefig(out, dpi=170)
plt.close(fig)
print(f"[saved] {out}")


# %% [markdown]
# ## 2. Coupling score
# Product of the two falloffs at coupling_distance_km = 5. The polar
# rim band shows up; the rest of the cap is essentially zero.

# %%
coupling_score = _open(PROCESSED / "scored" / "coupling_score_southpole_240m.tif")
arr = _coarsen(coupling_score)
finite = arr[np.isfinite(arr)]
print(
    f"coupling: max {finite.max():.3f}, median {np.median(finite):.4f}, "
    f"pct > 0.1 {100 * (finite > 0.1).sum() / finite.size:.2f}%"
)

fig, ax = plt.subplots(figsize=(8, 8))
arr_clip = np.where(arr > 0.001, arr, 0.001)  # avoid log(0)
im = ax.imshow(arr_clip, cmap="viridis", origin="upper", norm=LogNorm(vmin=0.001, vmax=1.0))
ax.set_title("Coupling score (log-scaled): high where near-PSR AND near-rim")
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im, ax=ax, label="coupling score (log)")
fig.tight_layout()
out = SANITY / "coupling_score.png"
fig.savefig(out, dpi=170)
plt.close(fig)
print(f"[saved] {out}")


# %% [markdown]
# ## 3. Coupling score with NASA candidates and our top-20 overlaid

# %%
nasa = regions_to_geodataframe()
sites = gpd.read_file(OUTPUTS / "top_sites.geojson")

to_polar = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, coupling_score.rio.crs, always_xy=True)
nasa_x, nasa_y = to_polar.transform(nasa["lon"].to_numpy(), nasa["lat"].to_numpy())
site_x, site_y = sites["x_m"].to_numpy(), sites["y_m"].to_numpy()

transform = coupling_score.rio.transform()
extent = [
    transform.c,
    transform.c + transform.a * coupling_score.sizes["x"],
    transform.f + transform.e * coupling_score.sizes["y"],
    transform.f,
]

fig, ax = plt.subplots(figsize=(9, 9))
arr_clip = np.where(arr > 0.001, arr, 0.001)
im = ax.imshow(
    arr_clip,
    cmap="viridis",
    origin="upper",
    extent=extent,
    norm=LogNorm(vmin=0.001, vmax=1.0),
)
fig.colorbar(im, ax=ax, label="coupling score (log)", fraction=0.045, pad=0.04)

ax.scatter(
    nasa_x,
    nasa_y,
    s=700,
    facecolors="none",
    edgecolors="#e63946",
    linewidths=2.0,
    label="NASA Artemis III candidates",
)
for x, y, name in zip(nasa_x, nasa_y, nasa["name"], strict=True):
    ax.annotate(
        name,
        (x, y),
        xytext=(8, -8),
        textcoords="offset points",
        fontsize=8,
        color="#e63946",
        weight="bold",
    )

ax.scatter(
    site_x,
    site_y,
    s=80,
    color="#06d6a0",
    edgecolor="white",
    linewidths=1.2,
    label="selene-base top sites (6-crit)",
    zorder=5,
)
for _, row in sites.iterrows():
    ax.annotate(
        str(int(row["rank"])),
        (row["x_m"], row["y_m"]),
        xytext=(4, 4),
        textcoords="offset points",
        color="white",
        fontsize=8,
        fontweight="bold",
        zorder=6,
    )

ax.set_xlabel("x (m, lunar south polar stereographic)")
ax.set_ylabel("y (m, lunar south polar stereographic)")
ax.set_title("Coupling score with NASA candidates and our top-20 overlaid")
ax.legend(loc="lower left", fontsize=9)
fig.tight_layout()
out = SANITY / "coupling_overlay.png"
fig.savefig(out, dpi=170)
plt.close(fig)
print(f"[saved] {out}")
