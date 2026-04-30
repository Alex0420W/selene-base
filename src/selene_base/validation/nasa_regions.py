"""NASA's nine announced Artemis III candidate landing regions.

Two parallel representations live here:

- :func:`regions_to_geodataframe` builds 15 km-radius **disk
  approximations** around centroids derived from NASA's October 2024
  public announcement. This is the legacy week-4 representation and
  what the disk-based polygon-inside metric (week 8) compares against.
- :func:`regions_polygons_to_geodataframe` loads the USGS
  officially-published **simplified region envelopes** (4-vertex
  quadrilaterals) shipped in
  ``validation/data/nasa_regions_polygons_usgs.geojson``. This is the
  authoritative machine-readable approximation introduced in week 10.
  The shipped GeoJSON is canonical for the project and should not be
  re-downloaded.

Both representations remain available because the project's
five-stage validation history (v0.1 through v1.1) used the disk
approximations, and the corresponding metrics still ship for context
and continuity.

Source for the USGS polygons: USGS Data Release 10.5066/P1MEQ6UK
(McClernan, M.T., 2024, "Down Selected Artemis III Candidate Landing
Site Navigational Grids"). The polygons are the simplified envelopes
of NASA's LROC QuickMap region definitions, not the full operational
landing footprints.
"""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

import geopandas as gpd
import numpy as np
from pyproj import Transformer
from shapely.geometry import Polygon

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

DEFAULT_RADIUS_KM = 15.0

USGS_POLYGONS_GEOJSON = Path(__file__).parent / "data" / "nasa_regions_polygons_usgs.geojson"
USGS_POLYGONS_SOURCE_CRS = "+proj=longlat +R=1737400 +no_defs"
USGS_POLYGONS_DOI = "10.5066/P1MEQ6UK"
USGS_REGION_NAMES: tuple[str, ...] = (
    "Nobile Rim 2",
    "Mons Mouton",
    "Malapert Massif",
    "de Gerlache Rim 2",
    "Mons Mouton Plateau",
    "Slater Plain",
    "Peak Near Cabeus B",
    "Nobile Rim 1",
    "Haworth",
)


class CandidateRegion(TypedDict):
    """One NASA Artemis III candidate region."""

    name: str
    lat: float
    lon: float
    radius_km: float


ARTEMIS_III_CANDIDATE_REGIONS: list[CandidateRegion] = [
    {"name": "Cabeus B", "lat": -82.3, "lon": -53.3, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Haworth", "lat": -86.5, "lon": -5.0, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Malapert Massif", "lat": -85.9, "lon": 2.9, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Mons Mouton", "lat": -84.6, "lon": -28.6, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Mons Mouton Plateau", "lat": -85.4, "lon": -25.0, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Nobile Rim 1", "lat": -85.5, "lon": 35.0, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Nobile Rim 2", "lat": -84.9, "lon": 31.5, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "de Gerlache Rim 2", "lat": -88.4, "lon": -62.4, "radius_km": DEFAULT_RADIUS_KM},
    {"name": "Slater Plain", "lat": -88.1, "lon": -54.3, "radius_km": DEFAULT_RADIUS_KM},
]


def regions_to_geodataframe(target_crs: str | None = None) -> gpd.GeoDataFrame:
    """Materialise NASA candidate centroids as disk polygons.

    Each region's disk is built in lunar south-polar stereographic
    metres (so the radius is exact in km) and then optionally
    reprojected to ``target_crs``. When ``target_crs`` is ``None`` the
    GeoDataFrame stays in :data:`selene_base.data.load.LUNAR_GEOGRAPHIC_CRS`.

    Args:
        target_crs: Target CRS as a PROJ string or EPSG code. ``None``
            (default) returns lunar geographic lon/lat polygons.

    Returns:
        GeoDataFrame with columns ``name``, ``lat`` (centroid),
        ``lon`` (centroid), ``radius_km``, and a ``geometry`` polygon.
    """
    polar_crs = (
        "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 "
        "+x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
    )
    geo_to_polar = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, polar_crs, always_xy=True)
    polar_to_geo = Transformer.from_crs(polar_crs, LUNAR_GEOGRAPHIC_CRS, always_xy=True)

    rows: list[dict[str, object]] = []
    n_segments = 64  # smooth disk
    angles = np.linspace(0.0, 2.0 * np.pi, n_segments, endpoint=False)

    for region in ARTEMIS_III_CANDIDATE_REGIONS:
        cx, cy = geo_to_polar.transform(region["lon"], region["lat"])
        radius_m = region["radius_km"] * 1000.0
        polar_xs = cx + radius_m * np.cos(angles)
        polar_ys = cy + radius_m * np.sin(angles)
        # Convert disk vertices back to lon/lat so the polygon lives in a
        # CRS that pyproj/geopandas can reproject without surprises.
        lons, lats = polar_to_geo.transform(polar_xs, polar_ys)
        polygon = Polygon(zip(lons, lats, strict=True))
        rows.append(
            {
                "name": region["name"],
                "lat": region["lat"],
                "lon": region["lon"],
                "radius_km": region["radius_km"],
                "geometry": polygon,
            }
        )

    gdf = gpd.GeoDataFrame(rows, crs=LUNAR_GEOGRAPHIC_CRS)
    if target_crs is not None and target_crs != str(LUNAR_GEOGRAPHIC_CRS):
        gdf = gdf.to_crs(target_crs)
    return gdf


def regions_polygons_to_geodataframe(target_crs: str | None = None) -> gpd.GeoDataFrame:
    """Load USGS-published Artemis III candidate region envelopes.

    The polygons are the **simplified** USGS envelopes of NASA's nine
    candidate regions — 4-vertex quadrilaterals in lunar planetocentric
    lon/lat space, sourced from USGS Data Release 10.5066/P1MEQ6UK
    (McClernan 2024). These are the authoritative machine-readable
    approximation of NASA's region geometries; they are *not* the full
    operational landing footprints.

    The on-disk GeoJSON does not carry CRS metadata, so the source CRS
    is set explicitly to :data:`USGS_POLYGONS_SOURCE_CRS`
    (``+proj=longlat +R=1737400 +no_defs``).

    Args:
        target_crs: Target CRS as a PROJ string or EPSG code. ``None``
            (default) returns lunar geographic lon/lat polygons.

    Returns:
        GeoDataFrame with columns ``Region`` (str), ``RegionCode``
        (two-letter abbreviation), ``Area_km2`` (float), and
        ``geometry`` (Polygon). Reprojected to ``target_crs`` if given.

    Raises:
        FileNotFoundError: If the bundled GeoJSON is missing.
    """
    if not USGS_POLYGONS_GEOJSON.exists():
        raise FileNotFoundError(
            f"USGS polygon GeoJSON not found at {USGS_POLYGONS_GEOJSON}; "
            "this file ships with the package and should not be missing."
        )

    gdf = gpd.read_file(USGS_POLYGONS_GEOJSON)
    gdf = gdf.set_crs(USGS_POLYGONS_SOURCE_CRS, allow_override=True)
    keep = [c for c in ("Region", "RegionCode", "Area_km2", "geometry") if c in gdf.columns]
    gdf = gdf[keep]
    if target_crs is not None and target_crs != USGS_POLYGONS_SOURCE_CRS:
        gdf = gdf.to_crs(target_crs)
    return gdf
