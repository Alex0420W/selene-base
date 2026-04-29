"""Proximity metrics between selene-base ranked sites and NASA candidates.

Distances are computed in lunar south-polar stereographic metres (so a
"km" is a real km on the lunar sphere of R = 1737400 m), independent of
whichever CRS the input GeoDataFrames carry.
"""

from __future__ import annotations

from typing import TypedDict

import geopandas as gpd
import numpy as np
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import Point

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)
DEFAULT_NEAR_KM = 25.0


class PerRegionResult(TypedDict):
    name: str
    nearest_site_id: str
    nearest_site_rank: int
    distance_km: float
    contains_top_site: bool


class PerSiteResult(TypedDict):
    site_id: str
    rank: int
    nearest_region: str
    distance_km: float
    inside_region: bool


class ProximityResult(TypedDict):
    n_top_sites: int
    n_nasa_regions: int
    sites_within_any_region: int
    sites_within_25km_of_region: int
    regions_with_a_top_site: int
    near_km: float
    per_region: list[PerRegionResult]
    per_site: list[PerSiteResult]


def _project_xy(geom_lon: np.ndarray, geom_lat: np.ndarray) -> np.ndarray:
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    xs, ys = transformer.transform(geom_lon, geom_lat)
    return np.column_stack([xs, ys])


def proximity_analysis(
    top_sites: gpd.GeoDataFrame,
    nasa_regions: gpd.GeoDataFrame,
    *,
    near_km: float = DEFAULT_NEAR_KM,
) -> ProximityResult:
    """Compute alignment metrics between top sites and NASA candidates.

    For each top-N site we find the nearest NASA region (by centroid in
    polar metres) and check whether the site falls inside the region's
    polygon. For each NASA region we find the nearest top-N site and
    report whether at least one top site lies inside the region.

    Args:
        top_sites: GeoDataFrame with at minimum ``site_id``, ``rank``,
            ``lat``, ``lon`` columns. Geometry is unused; we use lat/lon
            so the result is independent of how the file was written.
        nasa_regions: GeoDataFrame from
            :func:`selene_base.validation.nasa_regions.regions_to_geodataframe`.
        near_km: Distance threshold (km) for the
            ``sites_within_25km_of_region`` headline number.

    Returns:
        :class:`ProximityResult` dict.
    """
    n_sites = len(top_sites)
    n_regions = len(nasa_regions)

    if n_sites == 0 or n_regions == 0:
        return {
            "n_top_sites": n_sites,
            "n_nasa_regions": n_regions,
            "sites_within_any_region": 0,
            "sites_within_25km_of_region": 0,
            "regions_with_a_top_site": 0,
            "near_km": near_km,
            "per_region": [],
            "per_site": [],
        }

    # Reproject to lon/lat for consistent containment checks. The NASA
    # GeoDataFrame is already in lunar geographic when fresh from
    # regions_to_geodataframe, but accept any CRS by reprojecting.
    nasa_geo = nasa_regions.to_crs(LUNAR_GEOGRAPHIC_CRS)
    sites_lat = top_sites["lat"].to_numpy()
    sites_lon = top_sites["lon"].to_numpy()
    sites_xy = _project_xy(sites_lon, sites_lat)
    region_centroids_xy = _project_xy(
        nasa_geo["lon"].to_numpy(),
        nasa_geo["lat"].to_numpy(),
    )

    site_tree = cKDTree(sites_xy)
    region_tree = cKDTree(region_centroids_xy)

    # Per-site: nearest region, distance, inside-any-region check.
    site_to_region_dist, site_to_region_idx = region_tree.query(sites_xy, k=1)
    site_inside = np.zeros(n_sites, dtype=bool)
    for i, (lat, lon) in enumerate(zip(sites_lat, sites_lon, strict=True)):
        pt = Point(lon, lat)
        site_inside[i] = bool(nasa_geo.contains(pt).any())

    per_site: list[PerSiteResult] = []
    for i in range(n_sites):
        per_site.append(
            {
                "site_id": str(top_sites["site_id"].iloc[i]),
                "rank": int(top_sites["rank"].iloc[i]),
                "nearest_region": str(nasa_geo["name"].iloc[int(site_to_region_idx[i])]),
                "distance_km": float(site_to_region_dist[i] / 1000.0),
                "inside_region": bool(site_inside[i]),
            }
        )

    # Per-region: nearest site, distance, contains-any-top-site check.
    region_to_site_dist, region_to_site_idx = site_tree.query(region_centroids_xy, k=1)
    region_contains = np.zeros(n_regions, dtype=bool)
    for r in range(n_regions):
        polygon = nasa_geo.geometry.iloc[r]
        for s in range(n_sites):
            if polygon.contains(Point(sites_lon[s], sites_lat[s])):
                region_contains[r] = True
                break

    per_region: list[PerRegionResult] = []
    for r in range(n_regions):
        idx = int(region_to_site_idx[r])
        per_region.append(
            {
                "name": str(nasa_geo["name"].iloc[r]),
                "nearest_site_id": str(top_sites["site_id"].iloc[idx]),
                "nearest_site_rank": int(top_sites["rank"].iloc[idx]),
                "distance_km": float(region_to_site_dist[r] / 1000.0),
                "contains_top_site": bool(region_contains[r]),
            }
        )

    near_threshold_m = near_km * 1000.0
    return {
        "n_top_sites": n_sites,
        "n_nasa_regions": n_regions,
        "sites_within_any_region": int(site_inside.sum()),
        "sites_within_25km_of_region": int(np.sum(site_to_region_dist <= near_threshold_m)),
        "regions_with_a_top_site": int(region_contains.sum()),
        "near_km": near_km,
        "per_region": per_region,
        "per_site": per_site,
    }


def render_summary(result: ProximityResult) -> str:
    """Compact stdout table for ``selene validate``."""
    n_sites = result["n_top_sites"]
    n_regions = result["n_nasa_regions"]
    inside = result["sites_within_any_region"]
    near = result["sites_within_25km_of_region"]
    matched = result["regions_with_a_top_site"]
    lines: list[str] = []
    lines.append(f"top {n_sites} sites vs {n_regions} NASA candidates:")
    lines.append(f"  inside any region:                       {inside:>3} / {n_sites}")
    near_label = f"  within {result['near_km']:.0f} km of any centroid:"
    lines.append(f"{near_label:<43}{near:>3} / {n_sites}")
    lines.append(f"  regions matched (>=1 top site inside):    {matched:>3} / {n_regions}")
    lines.append("")
    lines.append(f"{'NASA region':<22} {'nearest':<10} {'dist (km)':>10}  inside?")
    for row in result["per_region"]:
        flag = "yes" if row["contains_top_site"] else "no"
        lines.append(
            f"{row['name']:<22} {row['nearest_site_id']:<10} {row['distance_km']:>10.1f}  {flag}"
        )
    return "\n".join(lines)
