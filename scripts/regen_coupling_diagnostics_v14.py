"""Regenerate the two coupling-score diagnostic plots at v1.4.0 state.

Outputs:

- ``docs/img/coupling_overlay.png``: coupling score (log-scaled)
  background, USGS-published Artemis III polygons (red outlines),
  selene-base v1.4 global-ranking top-20 (cyan dots, numbered).
- ``docs/img/coupling_score.png``: coupling score with USGS polygon
  outlines and an explicit south-pole marker for spatial context, and
  the sparsity caption rendered into the image.

The week-7 ``notebooks/05_coupling_visualization.py`` notebook is left
untouched as a documentation snapshot of the v1.0-era state; this
script is the v1.4 regenerator.

Run from the repo root:

    python scripts/regen_coupling_diagnostics_v14.py

Both PNGs land directly in ``docs/img/`` (where the README references
them), not the throwaway ``data/outputs/sanity/`` directory.
"""

from __future__ import annotations

import os

# Same lunar/Earth ellipsoid override as the week 10/11/12 notebooks —
# the coupling COG carries an "unknown" ellipsoid wrapped around
# R = 1737400, and the USGS polygons use an explicit lunar PROJ string;
# pyproj refuses to operate across them otherwise.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from pyproj import Transformer

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

REPO_ROOT = Path(__file__).resolve().parent.parent
COUPLING_SCORE_TIF = (
    REPO_ROOT / "data" / "processed" / "scored" / "coupling_score_southpole_240m.tif"
)
TOP_SITES_GEOJSON = REPO_ROOT / "data" / "outputs" / "top_sites.geojson"
DOCS_IMG = REPO_ROOT / "docs" / "img"
DOCS_IMG.mkdir(parents=True, exist_ok=True)

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


def _open(path: Path):
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


def _coarsen(da, factor: int = 4) -> np.ndarray:
    return da.coarsen(y=factor, x=factor, boundary="trim").mean().to_numpy()


# Manual label offsets, in projected metres, for the four regions
# clustered near the pole. Positive dx pushes east on screen; positive
# dy pushes north on screen. Connector lines are drawn for any region
# whose offset exceeds CONNECTOR_THRESHOLD_M.
LABEL_OFFSETS_M: dict[str, tuple[float, float]] = {
    # Five polygons cluster within ~80 km of each other near the pole;
    # place each label along a different azimuth from its polygon
    # centroid so the boxes don't collide.
    # MM polygon centre ≈ (90 km, 115 km): label south-west, well below
    # the cluster.
    "Mons Mouton": (-55_000.0, -85_000.0),
    "Mons Mouton Plateau": (-25_000.0, 110_000.0),  # north (large polygon)
    "Malapert Massif": (-95_000.0, -25_000.0),  # far west
    # N1 polygon centre ≈ (110 km, 110 km): south-east of polygon, but
    # 40 km below MM's label so the two don't touch.
    "Nobile Rim 1": (75_000.0, -75_000.0),
    "Nobile Rim 2": (95_000.0, 30_000.0),  # far east
    "Haworth": (-110_000.0, 25_000.0),
    "Peak Near Cabeus B": (-25_000.0, 40_000.0),
    "de Gerlache Rim 2": (-75_000.0, -55_000.0),
    "Slater Plain": (45_000.0, -45_000.0),
}
CONNECTOR_THRESHOLD_M = 25_000.0


def _draw_polygon_outlines(ax, polygons: gpd.GeoDataFrame, *, color: str, linewidth: float) -> None:
    """Plot polygon boundaries (no fill) on a matplotlib axis."""
    polygons.boundary.plot(ax=ax, color=color, linewidth=linewidth)


def _label_polygons(ax, polygons: gpd.GeoDataFrame, *, color: str) -> None:
    """Place a region-name label per polygon, with manual offsets and
    connector lines for crowded regions near the pole."""
    for _, row in polygons.iterrows():
        name = str(row["Region"])
        centroid = row.geometry.centroid
        cx, cy = centroid.x, centroid.y
        dx, dy = LABEL_OFFSETS_M.get(name, (0.0, 0.0))
        label_x, label_y = cx + dx, cy + dy
        offset_m = (dx * dx + dy * dy) ** 0.5
        if offset_m > CONNECTOR_THRESHOLD_M:
            ax.plot(
                [cx, label_x],
                [cy, label_y],
                color=color,
                linewidth=0.5,
                alpha=0.6,
                zorder=4,
            )
        ax.text(
            label_x,
            label_y,
            name,
            color=color,
            fontsize=8,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "none", "alpha": 0.85},
            zorder=5,
        )


def main() -> None:
    coupling = _open(COUPLING_SCORE_TIF).rio.write_crs(POLAR_PROJ, inplace=False)
    arr = _coarsen(coupling)
    arr_clip = np.where(arr > 1e-3, arr, 1e-3)

    # Sparsity numbers are computed on the *full-resolution* raster
    # (not the coarsened one used for plotting) so they match the
    # README and the criteria/coupling.py module docstring; mean-
    # coarsening would attenuate the very-low-but-nonzero counts.
    full = coupling.to_numpy()
    finite_full = full[np.isfinite(full)]
    pct_above_0 = 100.0 * (finite_full > 0.0).sum() / finite_full.size
    pct_above_p1 = 100.0 * (finite_full > 0.1).sum() / finite_full.size
    print(
        f"coupling sparsity (full-resolution): {pct_above_0:.2f}% > 0.0, {pct_above_p1:.2f}% > 0.1"
    )

    transform = coupling.rio.transform()
    extent = (
        transform.c,
        transform.c + transform.a * coupling.sizes["x"],
        transform.f + transform.e * coupling.sizes["y"],
        transform.f,
    )

    polygons = regions_polygons_to_geodataframe(target_crs=POLAR_PROJ)

    # ---------------- coupling_overlay.png (top-20 + USGS polygons) ----------------
    sites = gpd.read_file(TOP_SITES_GEOJSON)
    to_polar = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    if "x_m" in sites.columns and "y_m" in sites.columns:
        site_x = sites["x_m"].to_numpy()
        site_y = sites["y_m"].to_numpy()
    else:
        site_x, site_y = to_polar.transform(sites["lon"].to_numpy(), sites["lat"].to_numpy())

    fig, ax = plt.subplots(figsize=(9, 9))
    im = ax.imshow(
        arr_clip,
        cmap="viridis",
        origin="upper",
        extent=extent,
        norm=LogNorm(vmin=1e-3, vmax=1.0),
    )
    fig.colorbar(im, ax=ax, label="coupling score (log)", fraction=0.045, pad=0.04)

    _draw_polygon_outlines(ax, polygons, color="#e63946", linewidth=1.8)
    _label_polygons(ax, polygons, color="#e63946")

    ax.scatter(
        site_x,
        site_y,
        s=70,
        color="#22d3ee",
        edgecolor="black",
        linewidths=0.7,
        zorder=10,
    )
    for _, row in sites.iterrows():
        if "x_m" in row:
            rx, ry = float(row["x_m"]), float(row["y_m"])
        else:
            rx, ry = to_polar.transform(row["lon"], row["lat"])
            rx, ry = float(rx), float(ry)
        ax.annotate(
            str(int(row["rank"])),
            (rx, ry),
            xytext=(5, 5),
            textcoords="offset points",
            color="white",
            fontsize=7,
            fontweight="bold",
            zorder=11,
            path_effects=None,
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color="#e63946",
            linewidth=1.8,
            label="NASA Artemis III regions (USGS polygons)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#22d3ee",
            markeredgecolor="black",
            markersize=8,
            label="selene-base top-20 (7-criterion, global ranking)",
        ),
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=9, framealpha=0.92)
    ax.set_xlabel("x (m, lunar south polar stereographic)")
    ax.set_ylabel("y (m, lunar south polar stereographic)")
    ax.set_title("Coupling score with USGS polygons and selene-base v1.4 top-20 overlaid")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    fig.tight_layout()
    out = DOCS_IMG / "coupling_overlay.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")

    # ---------------- coupling_score.png (sparsity diagnostic) ----------------
    fig, ax = plt.subplots(figsize=(9, 9))
    im = ax.imshow(
        arr_clip,
        cmap="viridis",
        origin="upper",
        extent=extent,
        norm=LogNorm(vmin=1e-3, vmax=1.0),
    )
    fig.colorbar(im, ax=ax, label="coupling score (log)", fraction=0.045, pad=0.04)

    _draw_polygon_outlines(ax, polygons, color="#e63946", linewidth=1.0)

    # South-pole marker.
    ax.scatter(
        [0.0],
        [0.0],
        marker="*",
        s=200,
        color="white",
        edgecolor="black",
        linewidths=1.0,
        zorder=10,
    )
    ax.annotate(
        "South Pole",
        (0.0, 0.0),
        xytext=(10, 10),
        textcoords="offset points",
        color="white",
        fontsize=10,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.25", "fc": "black", "ec": "white", "alpha": 0.7},
        zorder=11,
    )

    ax.set_xlabel("x (m, lunar south polar stereographic)")
    ax.set_ylabel("y (m, lunar south polar stereographic)")
    ax.set_title("Coupling score: sparse rim band where near-PSR meets near-sunlit-ridge")
    # Sparsity caption rendered into the image, lower-right corner.
    ax.text(
        0.98,
        0.02,
        f"{pct_above_0:.2f}% of cells > 0.0   |   {pct_above_p1:.2f}% > 0.1",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="white",
        bbox={"boxstyle": "round,pad=0.4", "fc": "black", "ec": "white", "alpha": 0.75},
    )
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    fig.tight_layout()
    out = DOCS_IMG / "coupling_score.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")


if __name__ == "__main__":
    main()
