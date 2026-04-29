# %% [markdown]
# # Week 2 — slope first pass
#
# Loads the cached LOLA COG (written by `selene preprocess`), derives slope
# in degrees on the common 240 m grid, applies the slope criterion, and
# saves three side-by-side plots to `data/outputs/sanity/slope.png`:
#
# 1. elevation (m)
# 2. slope (°)
# 3. slope criterion score, [0, 1]
#
# Run after:
# ```
# selene download lola
# selene preprocess
# ```

# %%
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401
import xarray as xr
import yaml

from selene_base.criteria import slope as slope_criterion

PROCESSED = Path("data/processed")
SANITY_DIR = Path("data/outputs/sanity")
SANITY_DIR.mkdir(parents=True, exist_ok=True)

REGION_CFG = Path("config/region_southpole.yaml")
with REGION_CFG.open() as fh:
    cfg = yaml.safe_load(fh)
PIXEL_M = float(cfg["resolution_m"])


# %% [markdown]
# ## Load the cached LOLA elevation COG

# %%
lola_cog = PROCESSED / "lola_southpole_240m.tif"
if not lola_cog.exists():
    raise SystemExit(f"missing {lola_cog}; run `selene download lola` then `selene preprocess`")

elevation = (
    rioxarray.open_rasterio(lola_cog, masked=True).squeeze("band", drop=True).rename("elevation_m")
)
print(f"shape: {dict(elevation.sizes)}")
print(f"crs:   {elevation.rio.crs}")
e_arr = elevation.to_numpy()
finite = e_arr[np.isfinite(e_arr)]
print(
    f"elevation min/median/max: {finite.min():.0f} / {np.median(finite):.0f} / {finite.max():.0f} m"
)


# %% [markdown]
# ## Derive slope and compute the score

# %%
slope_deg = slope_criterion.derive_slope_degrees(elevation, pixel_size_m=PIXEL_M)
score = slope_criterion.compute(slope_deg)

s_arr = slope_deg.to_numpy()
sc_arr = score.to_numpy()
finite_s = s_arr[np.isfinite(s_arr)]
finite_sc = sc_arr[np.isfinite(sc_arr)]
print(
    f"slope (°) min/median/max: "
    f"{finite_s.min():.2f} / {np.median(finite_s):.2f} / {finite_s.max():.2f}"
)
print(
    f"score min/median/max: "
    f"{finite_sc.min():.3f} / {np.median(finite_sc):.3f} / {finite_sc.max():.3f}"
)
print(f"cells with score > 0.7: {100.0 * (finite_sc > 0.7).sum() / finite_sc.size:.1f}%")


# %% [markdown]
# ## Side-by-side plots
# Coarsen each panel before plotting so the figure renders quickly.


# %%
def _coarsen(da: xr.DataArray, factor: int) -> np.ndarray:
    return da.coarsen(y=factor, x=factor, boundary="trim").mean().to_numpy()


COARSE = 4

fig, axes = plt.subplots(1, 3, figsize=(16, 6))

ax = axes[0]
im = ax.imshow(_coarsen(elevation, COARSE), cmap="terrain", origin="upper")
ax.set_title("LOLA elevation (m)")
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im, ax=ax, label="m")

ax = axes[1]
im = ax.imshow(_coarsen(slope_deg, COARSE), cmap="magma", origin="upper", vmin=0, vmax=30)
ax.set_title("slope (°)")
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im, ax=ax, label="°")

ax = axes[2]
im = ax.imshow(_coarsen(score, COARSE), cmap="viridis", origin="upper", vmin=0, vmax=1)
ax.set_title("slope criterion score (max=15°)")
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(im, ax=ax, label="score")

fig.tight_layout()
out = SANITY_DIR / "slope.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"[saved] {out}")
