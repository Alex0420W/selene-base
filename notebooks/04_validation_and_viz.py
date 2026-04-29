# %% [markdown]
# # Week 4 — validation and visualisation
#
# Runs the full `selene validate` + `selene viz` pipeline on whatever
# the previous weeks left in `data/processed/` and `data/outputs/`, and
# emits the headline figures used by the README:
#
# - `data/outputs/sanity/webmap_screenshot.png` — static screenshot of
#   the aggregate score with NASA candidates and ranked sites overlaid.
# - `data/outputs/sanity/validation_table.png` — per-region distance
#   table.
# - `data/outputs/validation.json` — full proximity result.
# - `data/outputs/webmap.html` — interactive web map.
# - `data/outputs/sites/index.html` + 20 per-site HTML reports.
#
# Run after:
# ```
# selene download lola && selene download illumination && selene download robbins
# selene preprocess && selene score && selene rank --top-n 20 --min-distance-km 25
# ```

# %%
from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray  # noqa: F401
from pyproj import Transformer

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.pipeline import validate as _validate
from selene_base.pipeline import viz as _viz
from selene_base.validation.comparison import render_summary
from selene_base.validation.nasa_regions import regions_to_geodataframe

OUTPUTS = Path("data/outputs")
PROCESSED = Path("data/processed")
SANITY = OUTPUTS / "sanity"
SANITY.mkdir(parents=True, exist_ok=True)


# %% [markdown]
# ## 1. Validation against NASA's nine candidates

# %%
result = _validate.run(outputs_dir=OUTPUTS)
print()
print(render_summary(result))


# %% [markdown]
# ## 2. Static map screenshot for the README
#
# Folium's HTML doesn't render in markdown previews, so we build an
# equivalent matplotlib screenshot — same data, same colour scheme —
# and save it as a PNG.

# %%
score = rioxarray.open_rasterio(OUTPUTS / "score_southpole.tif", masked=True).squeeze(
    "band", drop=True
)
sites = gpd.read_file(OUTPUTS / "top_sites.geojson")
nasa = regions_to_geodataframe()

# Project NASA centroids into the score map's CRS for plotting.
to_polar = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, score.rio.crs, always_xy=True)
nasa_x, nasa_y = to_polar.transform(nasa["lon"].to_numpy(), nasa["lat"].to_numpy())

# Coarsen for a fast plot.
factor = 6
arr = score.coarsen(y=factor, x=factor, boundary="trim").mean().to_numpy()
transform = score.rio.transform()
extent = [
    transform.c,
    transform.c + transform.a * score.sizes["x"],
    transform.f + transform.e * score.sizes["y"],
    transform.f,
]

fig, ax = plt.subplots(figsize=(9, 9))
im = ax.imshow(arr, cmap="plasma", origin="upper", extent=extent, vmin=0, vmax=1)
fig.colorbar(im, ax=ax, label="aggregate suitability score", fraction=0.045, pad=0.04)

# NASA Artemis III candidates (red rings).
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

# Top sites (cyan dots with rank labels).
ax.scatter(
    sites["x_m"],
    sites["y_m"],
    s=80,
    color="#06d6a0",
    edgecolor="white",
    linewidths=1.2,
    label="selene-base top sites",
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
ax.set_title("selene-base top 20 candidates vs NASA Artemis III regions")
ax.legend(loc="lower left", fontsize=9)
fig.tight_layout()
screenshot_path = SANITY / "webmap_screenshot.png"
fig.savefig(screenshot_path, dpi=170)
plt.close(fig)
print(f"[saved] {screenshot_path}")


# %% [markdown]
# ## 3. Per-region distance table screenshot

# %%
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.axis("off")
table_data = [
    [
        r["name"],
        r["nearest_site_id"],
        f"{r['distance_km']:.1f}",
        "yes" if r["contains_top_site"] else "no",
    ]
    for r in result["per_region"]
]
table = ax.table(
    cellText=table_data,
    colLabels=["NASA candidate", "nearest site", "distance (km)", "contains top site?"],
    loc="center",
    cellLoc="center",
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.0, 1.5)
ax.set_title("Distance from each NASA Artemis III candidate to nearest selene-base top site")
fig.tight_layout()
table_path = SANITY / "validation_table.png"
fig.savefig(table_path, dpi=170, bbox_inches="tight")
plt.close(fig)
print(f"[saved] {table_path}")


# %% [markdown]
# ## 4. Build interactive web map + per-site reports

# %%
artefacts = _viz.run(outputs_dir=OUTPUTS, processed_dir=PROCESSED)
print()
for name, path in artefacts.items():
    print(f"  {name:<14} {path}")


# %% [markdown]
# ## 5. Numbers for the README abstract

# %%
inside = result["sites_within_any_region"]
near = result["sites_within_25km_of_region"]
matched = result["regions_with_a_top_site"]
n_sites = result["n_top_sites"]
n_regions = result["n_nasa_regions"]
print(
    f"top-{n_sites} sites: {inside}/{n_sites} inside any NASA region; "
    f"{near}/{n_sites} within 25 km of any centroid; "
    f"{matched}/{n_regions} regions matched"
)
nearest = sorted(result["per_region"], key=lambda r: r["distance_km"])[:3]
for r in nearest:
    print(f"  closest: {r['name']:<22} -> {r['nearest_site_id']} ({r['distance_km']:.1f} km)")

with (OUTPUTS / "validation.json").open() as fh:
    saved = json.load(fh)
assert saved == result, "validation.json round-trip mismatch"
print("validation.json round-trip OK")
