"""Score Wueller 2026 sites against selene-base's per-criterion rasters.

Companion to :mod:`selene_base.validation.wueller_comparison`. The
spatial comparison there answers "are the two catalogs picking the same
*cells*?"; this module answers "do the two catalogs agree on what makes
a *good* cell?" by sampling selene's six active 240 m criterion rasters
plus the aggregate score at each Wueller-site coordinate.

Per v1.5's design (see :mod:`selene_base.pipeline.rank_per_region_tiled`
header), per-criterion rasters are not materialised globally at 20 m —
the v1.5 ranker upsamples the 240 m product onto each tile via
bilinear ``reproject_match``. This module therefore evaluates against
the 240 m criterion rasters that drive selene's aggregate score, which
is the same per-criterion stack v1.5 itself feeds into the score
upsample. v1.5's per-region tiled framing reduces the median
matched-pair distance from 1.88 km to 1.76 km without changing the
80 % → 81 % selene-to-Wueller match rate; the methodology converges at
20 m, so 240 m criterion evaluation is representative of selene's
view of any in-scope cell.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (registers .rio accessor on xarray)
import xarray as xr
from pyproj import Transformer

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS
from selene_base.scoring.ranking import (
    HLS_DTE_VISIBILITY_MIN,
    HLS_ILLUMINATION_MIN,
    HLS_SLOPE_MAX_DEG,
)
from selene_base.validation.wueller_comparison import POLAR_PROJ, load_wueller_sites

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_OUTPUTS_DIR = Path("data/outputs")
DEFAULT_SCORED_SUBDIR = "scored"

CRITERION_RASTER_FILENAMES: dict[str, str] = {
    "slope": "slope_score_southpole_240m.tif",
    "illumination": "illumination_score_southpole_240m.tif",
    "eva_psr_access": "eva_psr_access_score_southpole_240m.tif",
    "thermal": "thermal_score_southpole_240m.tif",
    "multi_volatile": "multi_volatile_score_southpole_240m.tif",
    "los_to_earth": "los_to_earth_score_southpole_240m.tif",
    "pgv_seismic": "pgv_seismic_score_southpole_240m.tif",
}

RAW_RASTER_FILENAMES: dict[str, str] = {
    "slope_deg": "lola_slope_deg_southpole_240m.tif",
    "illumination": "illumination_southpole_240m.tif",
    "los_visibility": "los_visibility_fraction_southpole_240m.tif",
}

AGGREGATE_FILENAME = "score_southpole.tif"


def _open_raster(path: Path) -> xr.DataArray:
    """Open a 240 m COG as a polar-stereographic 2D DataArray."""
    da = rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)
    return da


def _sample_at_xy(da: xr.DataArray, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Nearest-pixel sampling at projected (x, y) coordinates."""
    pts = da.sel(
        x=xr.DataArray(xs, dims="point"),
        y=xr.DataArray(ys, dims="point"),
        method="nearest",
    )
    arr = pts.to_numpy()
    return arr.astype(np.float64)


def score_wueller_sites(
    *,
    wueller_sites: gpd.GeoDataFrame | None = None,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    outputs_dir: Path = DEFAULT_OUTPUTS_DIR,
    scored_subdir: str = DEFAULT_SCORED_SUBDIR,
    in_scope_only: bool = True,
) -> pd.DataFrame:
    """Evaluate selene's per-criterion rasters at every Wueller-site coordinate.

    Args:
        wueller_sites: Pre-loaded Wueller GeoDataFrame, or ``None`` to
            load the bundled real catalog via
            :func:`load_wueller_sites`. The frame must carry
            ``wueller_site_id``, ``region``, ``lat``, ``lon``, and
            ``in_usgs_scope`` columns.
        processed_dir: Directory holding the 240 m raw + scored COGs.
        outputs_dir: Directory holding the aggregate ``score_southpole.tif``.
        scored_subdir: Subdirectory under ``processed_dir`` holding the
            per-criterion score COGs.
        in_scope_only: When ``True`` (default), evaluate only the
            in-scope sites (region in NASA's Oct 2024 down-selected nine).

    Returns:
        DataFrame with one row per Wueller site and columns:
        ``wueller_site_id``, ``region``, ``in_usgs_scope``, ``lat``,
        ``lon``, ``x_m``, ``y_m``, ``aggregate_score``,
        per-criterion ``score_<crit>``, raw ``slope_deg``,
        ``illumination``, ``los_visibility``, plus a derived
        ``hls_compliant`` flag (slope ≤ 8 °, illum ≥ 0.33,
        DTE ≥ 0.50). The 100 m buffer-from-steep constraint is sub-pixel
        at 240 m and is therefore omitted here — it only becomes
        meaningful at v1.5's 20 m per-region tiles.
    """
    if wueller_sites is None:
        wueller_sites = load_wueller_sites()
    if in_scope_only:
        wueller_sites = wueller_sites[wueller_sites["in_usgs_scope"]].reset_index(drop=True)

    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    xs, ys = transformer.transform(wueller_sites["lon"].to_numpy(), wueller_sites["lat"].to_numpy())
    xs = np.asarray(xs, dtype=np.float64)
    ys = np.asarray(ys, dtype=np.float64)

    rows: dict[str, np.ndarray | list[object]] = {
        "wueller_site_id": wueller_sites["wueller_site_id"].astype(str).tolist(),
        "region": wueller_sites["region"].astype(str).tolist(),
        "in_usgs_scope": wueller_sites["in_usgs_scope"].astype(bool).tolist(),
        "lat": wueller_sites["lat"].to_numpy(),
        "lon": wueller_sites["lon"].to_numpy(),
        "x_m": xs,
        "y_m": ys,
    }

    scored_dir = processed_dir / scored_subdir
    for crit, fname in CRITERION_RASTER_FILENAMES.items():
        path = scored_dir / fname
        if not path.exists():
            rows[f"score_{crit}"] = np.full(len(wueller_sites), np.nan, dtype=np.float64)
            continue
        da = _open_raster(path)
        rows[f"score_{crit}"] = _sample_at_xy(da, xs, ys)

    for raw_name, fname in RAW_RASTER_FILENAMES.items():
        path = processed_dir / fname
        if not path.exists():
            rows[raw_name] = np.full(len(wueller_sites), np.nan, dtype=np.float64)
            continue
        da = _open_raster(path)
        rows[raw_name] = _sample_at_xy(da, xs, ys)

    agg_path = outputs_dir / AGGREGATE_FILENAME
    if agg_path.exists():
        da = _open_raster(agg_path)
        rows["aggregate_score"] = _sample_at_xy(da, xs, ys)
    else:
        rows["aggregate_score"] = np.full(len(wueller_sites), np.nan, dtype=np.float64)

    df = pd.DataFrame(rows)
    df["hls_compliant"] = (
        np.isfinite(df["slope_deg"])
        & np.isfinite(df["illumination"])
        & np.isfinite(df["los_visibility"])
        & (df["slope_deg"] <= HLS_SLOPE_MAX_DEG)
        & (df["illumination"] >= HLS_ILLUMINATION_MIN)
        & (df["los_visibility"] >= HLS_DTE_VISIBILITY_MIN)
    )

    column_order = [
        "wueller_site_id",
        "region",
        "in_usgs_scope",
        "lat",
        "lon",
        "x_m",
        "y_m",
        "aggregate_score",
        *[f"score_{c}" for c in CRITERION_RASTER_FILENAMES],
        "slope_deg",
        "illumination",
        "los_visibility",
        "hls_compliant",
    ]
    return df[column_order]
