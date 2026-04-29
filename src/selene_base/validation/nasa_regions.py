"""NASA's nine announced Artemis III candidate landing regions.

Coordinates are approximate centroids derived from NASA's public
site-selection materials (October 2024 announcement, current as of
2026). NASA's actual definitions are polygons that vary in shape and
size — we use disk approximations of ``radius_km = 15`` to compare
proximity, not for any claim of exact polygon match. The 15 km radius
is the publicly cited "operational region" scale around each centroid.

Don't use this list as authoritative geometry; it's a comparison-only
approximation for the validation step in :mod:`selene_base.validation.comparison`.
"""

from __future__ import annotations

from typing import TypedDict

import geopandas as gpd
import numpy as np
from pyproj import Transformer
from shapely.geometry import Polygon

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

DEFAULT_RADIUS_KM = 15.0


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
