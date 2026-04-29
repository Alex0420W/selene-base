# %% [markdown]
# # Week 1 data inventory
#
# Loads each dataset, prints summary stats, and saves a sanity-check plot
# to `data/outputs/sanity/<dataset>.png`. Run after `selene download all`.
#
# This file is in jupytext "percent" format — open it as a notebook in
# VS Code's interactive window, or run it as a script:
#
# ```bash
# python notebooks/01_data_inventory.py
# ```
#
# Each cell skips gracefully if the corresponding raw file is not
# present, so partial downloads still produce partial plots.

# %%
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

from selene_base.data.load import (
    LUNAR_GEOGRAPHIC_CRS,
    LUNAR_SOUTH_POLAR_CRS,
    load_crater_catalog,
    load_diviner,
    load_illumination,
    load_lend,
    load_lola_ldem,
)

RAW = Path("data/raw")
SANITY_DIR = Path("data/outputs/sanity")
SANITY_DIR.mkdir(parents=True, exist_ok=True)
print(f"writing sanity plots to {SANITY_DIR.resolve()}")


# %% [markdown]
# ## 1. Robbins crater catalog
# Scatter of south-polar crater centres in stereographic coordinates,
# coloured by diameter. Confirms the lat ≤ -75° filter took effect.

# %%
robbins_path = RAW / "robbins" / "robbins_southpole.csv.gz"
if not robbins_path.exists():
    print(f"[skip] {robbins_path} — run `selene download robbins`")
else:
    gdf = load_crater_catalog(robbins_path)
    print(f"rows: {len(gdf):,}")
    print(gdf[["lat", "lon", "diam_km"]].describe())

    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, LUNAR_SOUTH_POLAR_CRS, always_xy=True)
    x, y = transformer.transform(gdf["lon"].to_numpy(), gdf["lat"].to_numpy())

    fig, ax = plt.subplots(figsize=(7, 7))
    sc = ax.scatter(
        x / 1000.0,
        y / 1000.0,
        c=np.log10(gdf["diam_km"].to_numpy().clip(min=1e-3)),
        s=2,
        cmap="viridis",
        alpha=0.6,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x (km, south-polar stereographic)")
    ax.set_ylabel("y (km, south-polar stereographic)")
    ax.set_title(f"Robbins south-polar craters (n={len(gdf):,})")
    fig.colorbar(sc, ax=ax, label="log10 diameter (km)")
    fig.tight_layout()
    out = SANITY_DIR / "robbins.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[saved] {out}")


# %% [markdown]
# ## 2. LOLA LDEM elevation
# Colormap of elevation. Expect a roughly bowl-shaped polar profile —
# the south pole is in the South Pole–Aitken basin, so the map skews
# negative.

# %%
lola_path = RAW / "lola" / "ldem_80s_80m.img"
if not lola_path.exists():
    print(f"[skip] {lola_path} — run `selene download lola`")
else:
    da = load_lola_ldem(lola_path)
    print(f"shape: {dict(da.sizes)}")
    print(f"crs:   {da.rio.crs}")

    sub = da.coarsen(y=8, x=8, boundary="trim").mean().values
    finite = sub[np.isfinite(sub)]
    if finite.size > 0:
        print(
            f"elevation min/median/max: "
            f"{finite.min():.0f} / {np.median(finite):.0f} / {finite.max():.0f} m"
        )

    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(sub, cmap="terrain", origin="upper")
    ax.set_title("LOLA LDEM 80°S 80 m (×8 coarsened)")
    ax.set_xlabel("x pixels")
    ax.set_ylabel("y pixels")
    fig.colorbar(im, ax=ax, label="elevation (m)")
    fig.tight_layout()
    out = SANITY_DIR / "lola.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[saved] {out}")


# %% [markdown]
# ## 3. Diviner annual Tbol max / min
# Side-by-side maximum and minimum bolometric temperature. Tmin should
# show the cold traps inside permanently shadowed craters (≪ 100 K).

# %%
tmax_path = RAW / "diviner" / "diviner_tbol_max_sp.tif"
tmin_path = RAW / "diviner" / "diviner_tbol_min_sp.tif"
if not (tmax_path.exists() and tmin_path.exists()):
    print(f"[skip] {tmax_path} / {tmin_path} — run `selene download diviner`")
else:
    ds = load_diviner(tmax_path, tmin_path)
    for var in ("tbol_max", "tbol_min"):
        arr = ds[var].values
        finite = arr[np.isfinite(arr)]
        print(
            f"{var}: shape={arr.shape} "
            f"min={finite.min():.1f} K  median={np.median(finite):.1f} K  max={finite.max():.1f} K"
        )

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, var, cmap in zip(axes, ("tbol_max", "tbol_min"), ("inferno", "cividis"), strict=True):
        im = ax.imshow(ds[var].values, cmap=cmap, origin="upper")
        ax.set_title(f"Diviner {var}")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, label="K")
    fig.tight_layout()
    out = SANITY_DIR / "diviner.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[saved] {out}")


# %% [markdown]
# ## 4. Mazarico average illumination
# Bright Shackleton/Mons Mouton rims and dark crater floors are the
# expected pattern.

# %%
illum_path = RAW / "illumination" / "avgvisib_65s_240m_201608.img"
if not illum_path.exists():
    print(f"[skip] {illum_path} — run `selene download illumination`")
else:
    da = load_illumination(illum_path)
    print(f"shape: {dict(da.sizes)}")
    sub = da.coarsen(y=4, x=4, boundary="trim").mean().values
    finite = sub[np.isfinite(sub)]
    if finite.size > 0:
        print(
            f"illumination fraction min/median/max: "
            f"{finite.min():.3f} / {np.median(finite):.3f} / {finite.max():.3f}"
        )

    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(sub * 100.0, cmap="magma", origin="upper", vmin=0, vmax=100)
    ax.set_title("Mazarico average illumination (% of time)")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, label="illumination (%)")
    fig.tight_layout()
    out = SANITY_DIR / "illumination.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[saved] {out}")


# %% [markdown]
# ## 5. LEND epithermal neutron flux
# Low-flux regions correlate with hydrogen / putative water-ice. Expect
# faint suppression in/near PSRs; the absolute scale is instrument
# count rate, not a unit conversion.

# %%
lend_path = RAW / "lend" / "lend_csetn_sp.img"
if not lend_path.exists():
    print(f"[skip] {lend_path} — run `selene download lend`")
else:
    da = load_lend(lend_path)
    print(f"shape: {dict(da.sizes)}")
    arr = da.values
    finite = arr[np.isfinite(arr)]
    if finite.size > 0:
        lo, hi = np.percentile(finite, [2, 98])
        print(f"flux 2/98 percentiles: {lo:.3f} / {hi:.3f}")

    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(arr, cmap="viridis_r", origin="upper")
    ax.set_title("LEND epithermal neutron flux (low = H-rich)")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, label="flux / counts (native units)")
    fig.tight_layout()
    out = SANITY_DIR / "lend.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[saved] {out}")
