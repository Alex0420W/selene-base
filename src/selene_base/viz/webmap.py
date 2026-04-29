"""Folium-based standalone web map of ranked sites + NASA candidates.

Renders the aggregate suitability score as a translucent PNG overlay,
plots the top-N ranked sites with popups, draws NASA's nine Artemis III
candidate disks, and exposes everything through a single
:class:`folium.LayerControl` so reviewers can toggle layers in the
browser. Output is a single self-contained HTML file — open it offline,
no server, no CDN.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path

import folium
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401
import xarray as xr
from folium import FeatureGroup, LayerControl
from folium.raster_layers import ImageOverlay
from pyproj import Transformer

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)
DEFAULT_CENTER_LATLON = (-87.0, 0.0)
DEFAULT_ZOOM = 5


def _open_cog(path: Path) -> xr.DataArray:
    return rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)


def _png_data_uri(arr: np.ndarray, *, cmap: str, vmin: float, vmax: float) -> str:
    """Render ``arr`` as a transparent PNG data URI."""
    fig, ax = plt.subplots(figsize=(8, 8), dpi=120)
    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.imshow(arr, cmap=cmap, origin="upper", vmin=vmin, vmax=vmax, interpolation="nearest")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", transparent=True, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _polar_bounds_to_latlon(da: xr.DataArray) -> list[list[float]]:
    """Compute folium ``[[south, west], [north, east]]`` lat/lon bounds.

    The score COG is in polar metres; we sample its corner pixel centres
    in projected space and reproject to lon/lat for folium.ImageOverlay.
    """
    transform = da.rio.transform()
    h, w = da.sizes["y"], da.sizes["x"]
    xs = [transform.c, transform.c + transform.a * w]
    ys = [transform.f, transform.f + transform.e * h]
    transformer = Transformer.from_crs(POLAR_PROJ, LUNAR_GEOGRAPHIC_CRS, always_xy=True)
    lons, lats = transformer.transform(
        [xs[0], xs[1], xs[1], xs[0]],
        [ys[0], ys[0], ys[1], ys[1]],
    )
    south = float(np.min(lats))
    north = float(np.max(lats))
    west = float(np.min(lons))
    east = float(np.max(lons))
    return [[south, west], [north, east]]


def _coarsen_for_overlay(da: xr.DataArray, factor: int = 4) -> np.ndarray:
    return da.coarsen(y=factor, x=factor, boundary="trim").mean().to_numpy()


def _site_popup_html(row: gpd.GeoSeries) -> str:
    crit_cols = [c for c in row.index if c.startswith("score_")]
    parts: list[str] = []
    for c in crit_cols:
        value = row[c]
        if value is None or (isinstance(value, float) and value != value):
            cell = "n/a"
        else:
            cell = f"{float(value):.3f}"
        name = c.removeprefix("score_")
        parts.append(f"<tr><td>{name}</td><td style='text-align:right'>{cell}</td></tr>")
    rows = "".join(parts)
    return (
        "<div style='font-family:sans-serif;font-size:12px'>"
        f"<b>{row['site_id']}</b> &middot; rank {int(row['rank'])}<br>"
        f"lat {float(row['lat']):.3f} deg, lon {float(row['lon']):.3f} deg<br>"
        f"<b>score: {float(row['score']):.3f}</b>"
        f"<table style='margin-top:4px'>{rows}</table></div>"
    )


def _add_score_layer(
    m: folium.Map,
    cog_path: Path,
    *,
    name: str,
    cmap: str,
    bounds: list[list[float]] | None = None,
    show: bool = True,
) -> list[list[float]]:
    """Render the COG as an ImageOverlay layer and return its bounds."""
    da = _open_cog(cog_path)
    overlay_bounds = bounds if bounds is not None else _polar_bounds_to_latlon(da)
    arr = _coarsen_for_overlay(da, factor=4)
    finite = arr[np.isfinite(arr)]
    vmin = float(np.nanpercentile(finite, 1)) if finite.size else 0.0
    vmax = float(np.nanpercentile(finite, 99)) if finite.size else 1.0
    image_uri = _png_data_uri(arr, cmap=cmap, vmin=vmin, vmax=vmax)
    overlay = ImageOverlay(
        image=image_uri,
        bounds=overlay_bounds,
        opacity=0.75,
        name=name,
        show=show,
    )
    overlay.add_to(m)
    return overlay_bounds


def build_map(
    score_cog: Path,
    top_sites: gpd.GeoDataFrame,
    nasa_regions: gpd.GeoDataFrame,
    output_path: Path,
    *,
    processed_dir: Path | None = None,
) -> Path:
    """Generate an interactive Leaflet web map and save as standalone HTML.

    Layers:

    - **Aggregate suitability score** (raster overlay, plasma colormap,
      visible by default).
    - **Top-N candidate sites** — numbered cyan circle markers with a
      per-criterion popup table.
    - **NASA Artemis III candidate regions** — red disks with name
      tooltips.
    - **Per-criterion score layers** (slope, illumination, hazard) when
      their COGs are present in ``processed_dir / 'scored'``; off by
      default to keep the map fast on first open.

    Args:
        score_cog: Aggregate score COG path.
        top_sites: GeoDataFrame in lunar geographic CRS with at minimum
            ``site_id``, ``rank``, ``lat``, ``lon``, ``score`` columns.
        nasa_regions: GeoDataFrame from
            :func:`selene_base.validation.nasa_regions.regions_to_geodataframe`.
        output_path: Destination HTML path.
        processed_dir: Where per-criterion score COGs live; pass
            ``None`` to skip the per-criterion layers.

    Returns:
        ``output_path``.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Lunar Moon — folium can't fetch tiles from real-Earth servers and
    # render them on lunar geometry, so we use no base tiles and a plain
    # dark background. The score raster overlay is the basemap.
    m = folium.Map(
        location=DEFAULT_CENTER_LATLON,
        zoom_start=DEFAULT_ZOOM,
        tiles=None,
        control_scale=True,
        prefer_canvas=True,
    )
    folium.TileLayer(
        tiles="cartodbpositronnolabels",
        name="Earth basemap (positional reference only)",
        attr="© CartoDB | Moon data © NASA / LRO / selene-base",
        show=False,
    ).add_to(m)

    bounds = _add_score_layer(
        m,
        score_cog,
        name="Aggregate suitability score",
        cmap="plasma",
        show=True,
    )

    # Optional per-criterion overlays.
    scored_dir: Path | None = None
    if processed_dir is not None:
        scored_dir = Path(processed_dir) / "scored"
    if scored_dir is not None and scored_dir.exists():
        for crit in ("slope", "illumination", "hazard"):
            crit_path = scored_dir / f"{crit}_score_southpole_240m.tif"
            if crit_path.exists():
                _add_score_layer(
                    m,
                    crit_path,
                    name=f"{crit} score",
                    cmap="viridis",
                    bounds=bounds,
                    show=False,
                )

    # NASA candidate regions (in lon/lat already).
    nasa_geo = nasa_regions.to_crs(LUNAR_GEOGRAPHIC_CRS)
    nasa_group = FeatureGroup(name="NASA Artemis III candidates", show=True)
    for _, region in nasa_geo.iterrows():
        coords = [(lat, lon) for lon, lat in region.geometry.exterior.coords]
        folium.Polygon(
            locations=coords,
            color="#e63946",
            weight=2,
            fill=True,
            fill_color="#e63946",
            fill_opacity=0.20,
            tooltip=f"{region['name']} (NASA Artemis III)",
        ).add_to(nasa_group)
        folium.CircleMarker(
            location=(region["lat"], region["lon"]),
            radius=3,
            color="#e63946",
            fill=True,
            fill_color="#e63946",
            tooltip=region["name"],
        ).add_to(nasa_group)
    nasa_group.add_to(m)

    # Top sites.
    sites_group = FeatureGroup(name="selene-base top sites", show=True)
    for _, site in top_sites.iterrows():
        folium.CircleMarker(
            location=(float(site["lat"]), float(site["lon"])),
            radius=7,
            color="#06d6a0",
            weight=2,
            fill=True,
            fill_color="#06d6a0",
            fill_opacity=0.85,
            tooltip=f"#{int(site['rank'])} · {site['site_id']} · {site['score']:.3f}",
            popup=folium.Popup(_site_popup_html(site), max_width=300),
        ).add_to(sites_group)
        folium.map.Marker(
            (float(site["lat"]), float(site["lon"])),
            icon=folium.DivIcon(
                icon_size=(20, 20),
                icon_anchor=(10, -8),
                html=(
                    f"<div style='font-size:11px;font-weight:bold;color:white;"
                    f"text-shadow:0 0 3px black'>{int(site['rank'])}</div>"
                ),
            ),
        ).add_to(sites_group)
    sites_group.add_to(m)

    LayerControl(collapsed=False).add_to(m)

    m.fit_bounds(bounds)
    m.save(str(output_path))
    return output_path
