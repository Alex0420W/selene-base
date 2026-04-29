"""Proximity metrics between selene-base ranked sites and NASA candidates.

Two complementary metric families:

- **Centroid-distance** (legacy, week 4): for each site/region, distance
  from the site to the nearest NASA *centroid*. Cheap, but penalises
  any model that picks rim cells over centroid cells, which is what
  NASA itself does within each region.
- **Polygon-inside** (week 8): for each site/region, whether the site
  falls inside any 15 km NASA disk, plus the signed distance to the
  nearest disk *edge* (negative when the site is inside). Same disks
  as before; different geometric primitive.

Both families ship in a single :class:`ProximityResult` so a single
validation run produces both tables; downstream tooling can read whichever
key it cares about.

Distances are computed in lunar south-polar stereographic metres so a
"km" is a real km on the lunar sphere of R = 1 737 400 m.
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
    distance_km: float  # to centroid (legacy)
    distance_to_edge_km: float  # negative when nearest site is inside the disk
    contains_top_site: bool


class PerSiteResult(TypedDict):
    site_id: str
    rank: int
    nearest_region: str
    distance_km: float  # to centroid (legacy)
    distance_to_edge_km: float  # negative when site is inside the disk
    inside_region: bool  # legacy alias of inside_any_region
    inside_any_region: bool


class ProximityResult(TypedDict):
    n_top_sites: int
    n_nasa_regions: int
    near_km: float
    # Centroid-distance metrics (legacy):
    sites_within_any_region: int  # legacy alias of sites_inside_any_region
    sites_within_25km_of_region: int
    regions_with_a_top_site: int  # legacy alias of regions_containing_top_site
    # Polygon-based metrics (week 8):
    sites_inside_any_region: int
    regions_containing_top_site: int
    regions_with_top_site_within_disk_radius: int
    per_region: list[PerRegionResult]
    per_site: list[PerSiteResult]


def _project_xy(lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    xs, ys = transformer.transform(lons, lats)
    return np.column_stack([xs, ys])


def _signed_edge_distance_m(point: Point, polygon) -> float:
    """Signed distance (metres) from ``point`` to the nearest polygon edge.

    Positive when the point is outside the polygon; negative when inside.
    Both ``point`` and ``polygon`` must be in the same CRS — we use
    polar metres throughout this module.
    """
    edge_distance = polygon.boundary.distance(point)
    return -edge_distance if polygon.contains(point) else edge_distance


def proximity_analysis(
    top_sites: gpd.GeoDataFrame,
    nasa_regions: gpd.GeoDataFrame,
    *,
    near_km: float = DEFAULT_NEAR_KM,
) -> ProximityResult:
    """Compute alignment metrics — both centroid-distance and polygon-inside.

    Args:
        top_sites: GeoDataFrame with at minimum ``site_id``, ``rank``,
            ``lat``, ``lon`` columns.
        nasa_regions: GeoDataFrame from
            :func:`selene_base.validation.nasa_regions.regions_to_geodataframe`.
            Each row's geometry is the 15 km disk approximation in lunar
            geographic CRS.
        near_km: Threshold for the legacy "within X km of centroid"
            metric.

    Returns:
        :class:`ProximityResult` dict with both metric families.
    """
    n_sites = len(top_sites)
    n_regions = len(nasa_regions)

    if n_sites == 0 or n_regions == 0:
        return {
            "n_top_sites": n_sites,
            "n_nasa_regions": n_regions,
            "near_km": near_km,
            "sites_within_any_region": 0,
            "sites_within_25km_of_region": 0,
            "regions_with_a_top_site": 0,
            "sites_inside_any_region": 0,
            "regions_containing_top_site": 0,
            "regions_with_top_site_within_disk_radius": 0,
            "per_region": [],
            "per_site": [],
        }

    nasa_geo = nasa_regions.to_crs(LUNAR_GEOGRAPHIC_CRS)
    nasa_polar = nasa_regions.to_crs(POLAR_PROJ)
    region_radii_m = nasa_geo["radius_km"].astype(float).to_numpy() * 1000.0
    sites_lat = top_sites["lat"].to_numpy()
    sites_lon = top_sites["lon"].to_numpy()
    sites_xy = _project_xy(sites_lon, sites_lat)
    region_centroids_xy = _project_xy(
        nasa_geo["lon"].to_numpy(),
        nasa_geo["lat"].to_numpy(),
    )

    site_tree = cKDTree(sites_xy)
    region_tree = cKDTree(region_centroids_xy)

    # Per-site: nearest region (centroid distance), inside-disk check,
    # and signed distance to the nearest disk edge.
    site_to_region_dist, site_to_region_idx = region_tree.query(sites_xy, k=1)
    site_inside = np.zeros(n_sites, dtype=bool)
    site_to_edge_m = np.zeros(n_sites, dtype=np.float64)
    site_points_polar = [Point(xy[0], xy[1]) for xy in sites_xy]
    region_polygons_polar = [poly for poly in nasa_polar.geometry]
    for i, pt in enumerate(site_points_polar):
        # Inside-any-region: scan polygons for a hit.
        inside = False
        for polygon in region_polygons_polar:
            if polygon.contains(pt):
                inside = True
                break
        site_inside[i] = inside
        # Signed distance to the nearest disk edge:
        # - if the site is inside any polygon, return the *most negative*
        #   (deepest inside) so containment dominates the answer;
        # - if the site is outside every polygon, return the smallest
        #   positive (closest edge).
        # Both cases reduce to ``edge_distances.min()`` because negatives
        # always sort below positives.
        edge_distances = np.array(
            [_signed_edge_distance_m(pt, polygon) for polygon in region_polygons_polar]
        )
        site_to_edge_m[i] = float(edge_distances.min())

    per_site: list[PerSiteResult] = []
    for i in range(n_sites):
        per_site.append(
            {
                "site_id": str(top_sites["site_id"].iloc[i]),
                "rank": int(top_sites["rank"].iloc[i]),
                "nearest_region": str(nasa_geo["name"].iloc[int(site_to_region_idx[i])]),
                "distance_km": float(site_to_region_dist[i] / 1000.0),
                "distance_to_edge_km": float(site_to_edge_m[i] / 1000.0),
                "inside_region": bool(site_inside[i]),
                "inside_any_region": bool(site_inside[i]),
            }
        )

    # Per-region: nearest site (centroid + edge), contains-any-top-site.
    region_to_site_dist, region_to_site_idx = site_tree.query(region_centroids_xy, k=1)
    region_contains = np.zeros(n_regions, dtype=bool)
    region_to_site_edge_m = np.zeros(n_regions, dtype=np.float64)
    for r in range(n_regions):
        polygon = region_polygons_polar[r]
        # Does any top site fall inside this disk?
        for site_pt in site_points_polar:
            if polygon.contains(site_pt):
                region_contains[r] = True
                break
        # Signed distance from this disk's edge to the nearest top site
        # (negative if that nearest top site is inside the disk).
        nearest_pt = site_points_polar[int(region_to_site_idx[r])]
        region_to_site_edge_m[r] = _signed_edge_distance_m(nearest_pt, polygon)

    per_region: list[PerRegionResult] = []
    for r in range(n_regions):
        idx = int(region_to_site_idx[r])
        per_region.append(
            {
                "name": str(nasa_geo["name"].iloc[r]),
                "nearest_site_id": str(top_sites["site_id"].iloc[idx]),
                "nearest_site_rank": int(top_sites["rank"].iloc[idx]),
                "distance_km": float(region_to_site_dist[r] / 1000.0),
                "distance_to_edge_km": float(region_to_site_edge_m[r] / 1000.0),
                "contains_top_site": bool(region_contains[r]),
            }
        )

    near_threshold_m = near_km * 1000.0
    sites_within_centroid = int(np.sum(site_to_region_dist <= near_threshold_m))
    n_inside = int(site_inside.sum())
    n_regions_contain = int(region_contains.sum())

    # "within disk radius of edge" = inside the disk OR within disk_radius
    # of the boundary. Using each region's own radius lets us extend
    # gracefully if NASA ever publishes non-uniform disk sizes.
    nearest_site_to_edge_m = region_to_site_edge_m
    within_disk_radius = (nearest_site_to_edge_m <= 0) | (
        np.abs(nearest_site_to_edge_m) <= region_radii_m
    )
    n_regions_within_disk_radius = int(within_disk_radius.sum())

    return {
        "n_top_sites": n_sites,
        "n_nasa_regions": n_regions,
        "near_km": near_km,
        "sites_within_any_region": n_inside,
        "sites_within_25km_of_region": sites_within_centroid,
        "regions_with_a_top_site": n_regions_contain,
        "sites_inside_any_region": n_inside,
        "regions_containing_top_site": n_regions_contain,
        "regions_with_top_site_within_disk_radius": n_regions_within_disk_radius,
        "per_region": per_region,
        "per_site": per_site,
    }


def render_summary(result: ProximityResult) -> str:
    """Two-table stdout block for ``selene validate``.

    Prints the centroid-distance summary first (legacy headline), then
    the polygon-inside summary (week 8 headline). Per-region table at
    the bottom shows both metrics side by side.
    """
    n_sites = result["n_top_sites"]
    n_regions = result["n_nasa_regions"]
    near_km = result["near_km"]

    near = result["sites_within_25km_of_region"]
    inside = result["sites_inside_any_region"]
    contains = result["regions_containing_top_site"]
    within_disk = result["regions_with_top_site_within_disk_radius"]

    lines: list[str] = []
    lines.append(f"top {n_sites} sites vs {n_regions} NASA candidates:")
    lines.append("")
    lines.append("centroid-distance metrics (legacy, week 4):")
    label = f"  within {near_km:.0f} km of any centroid:"
    lines.append(f"{label:<48}{near:>3} / {n_sites}")
    lines.append("")
    lines.append("polygon-inside metrics (week 8):")
    lines.append(f"  inside any 15 km disk:                          {inside:>3} / {n_sites}")
    lines.append(f"  regions containing a top site:                  {contains:>3} / {n_regions}")
    lines.append(
        f"  regions with a top site within 1 disk radius:   {within_disk:>3} / {n_regions}"
    )
    lines.append("")
    lines.append(
        f"{'NASA region':<22} {'nearest':<10} {'dist-c (km)':>11} {'dist-edge (km)':>14}  inside?"
    )
    for row in result["per_region"]:
        flag = "yes" if row["contains_top_site"] else "no"
        lines.append(
            f"{row['name']:<22} {row['nearest_site_id']:<10} "
            f"{row['distance_km']:>11.1f} {row['distance_to_edge_km']:>14.1f}  {flag}"
        )
    return "\n".join(lines)
