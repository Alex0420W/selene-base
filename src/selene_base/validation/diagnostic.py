"""Per-criterion diagnostic: where do our top sites differ from NASA's?

Samples each criterion's score grid at our top-N site coordinates and
at the NASA candidate centroids, and reports mean ± std for both site
sets plus their signed difference. The output is the table the README
needs: which criteria favour our top sites, which favour NASA's, by how
much.

Two notes for the reader:

- The two-sample ``t_statistic`` returned alongside is informative,
  not inferential. With ``n=20`` against ``n=9`` and known structural
  differences in how the site sets were generated, a strict t-test
  hypothesis is the wrong frame; we report the statistic as a
  *ranking* signal that orders criteria by how strongly they
  separate the two samples.
- All sampling is nearest-pixel-centre, no interpolation. Score COGs
  are smooth enough that bilinear sampling would not change the
  per-region story.
"""

from __future__ import annotations

from collections.abc import Mapping

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS


def _sample_at_lonlat(
    grid: xr.DataArray,
    lons: np.ndarray,
    lats: np.ndarray,
) -> np.ndarray:
    """Sample ``grid`` at lon/lat points via nearest-pixel-centre lookup."""
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, grid.rio.crs, always_xy=True)
    xs, ys = transformer.transform(lons, lats)
    transform = grid.rio.transform()
    cols = np.round((np.asarray(xs) - transform.c) / transform.a - 0.5).astype(int)
    rows = np.round((np.asarray(ys) - transform.f) / transform.e - 0.5).astype(int)
    h, w = grid.sizes["y"], grid.sizes["x"]
    in_grid = (rows >= 0) & (rows < h) & (cols >= 0) & (cols < w)
    arr = grid.to_numpy()
    out = np.full(rows.shape, np.nan, dtype=np.float64)
    if in_grid.any():
        valid_rows = rows[in_grid]
        valid_cols = cols[in_grid]
        out[in_grid] = arr[valid_rows, valid_cols]
    return out


def per_criterion_comparison(
    top_sites: gpd.GeoDataFrame,
    nasa_regions: gpd.GeoDataFrame,
    score_maps: Mapping[str, xr.DataArray],
) -> pd.DataFrame:
    """Compute mean ± std at top sites vs at NASA centroids, per criterion.

    Args:
        top_sites: GeoDataFrame with at minimum ``lat`` and ``lon``
            columns. Empty input is rejected.
        nasa_regions: GeoDataFrame with ``lat`` and ``lon`` columns
            (centroids from ``regions_to_geodataframe``).
        score_maps: Mapping from criterion name to its [0, 1] score
            grid on the common projected CRS.

    Returns:
        DataFrame indexed by criterion name with columns:

        - ``our_top_n_mean``, ``our_top_n_std``
        - ``nasa_mean``, ``nasa_std``
        - ``delta`` (= ``our_top_n_mean - nasa_mean``)
        - ``abs_delta``
        - ``t_statistic`` (Welch's two-sample, informative-only)

    Raises:
        ValueError: When ``top_sites`` or ``score_maps`` is empty.
    """
    if len(top_sites) == 0:
        raise ValueError("top_sites is empty; nothing to compare")
    if len(score_maps) == 0:
        raise ValueError("score_maps is empty; nothing to compare")

    nasa_geo = nasa_regions.to_crs(LUNAR_GEOGRAPHIC_CRS)
    top_lons = top_sites["lon"].to_numpy()
    top_lats = top_sites["lat"].to_numpy()
    nasa_lons = nasa_geo["lon"].to_numpy()
    nasa_lats = nasa_geo["lat"].to_numpy()

    rows: list[dict[str, object]] = []
    for name, grid in score_maps.items():
        ours = _sample_at_lonlat(grid, top_lons, top_lats)
        theirs = _sample_at_lonlat(grid, nasa_lons, nasa_lats)
        ours_finite = ours[np.isfinite(ours)]
        theirs_finite = theirs[np.isfinite(theirs)]
        n_ours = ours_finite.size
        n_theirs = theirs_finite.size

        ours_mean = float(np.mean(ours_finite)) if n_ours else float("nan")
        ours_std = float(np.std(ours_finite, ddof=1)) if n_ours > 1 else float("nan")
        theirs_mean = float(np.mean(theirs_finite)) if n_theirs else float("nan")
        theirs_std = float(np.std(theirs_finite, ddof=1)) if n_theirs > 1 else float("nan")

        if n_ours > 1 and n_theirs > 1:
            num = ours_mean - theirs_mean
            denom = np.sqrt((ours_std**2) / n_ours + (theirs_std**2) / n_theirs)
            t_stat = float(num / denom) if denom > 0 else float("nan")
        else:
            t_stat = float("nan")

        rows.append(
            {
                "criterion": name,
                "our_top_n_mean": ours_mean,
                "our_top_n_std": ours_std,
                "nasa_mean": theirs_mean,
                "nasa_std": theirs_std,
                "delta": ours_mean - theirs_mean,
                "abs_delta": abs(ours_mean - theirs_mean),
                "t_statistic": t_stat,
                "n_ours": n_ours,
                "n_nasa": n_theirs,
            }
        )

    df = pd.DataFrame(rows).set_index("criterion")
    return df.sort_values("abs_delta", ascending=False)


def render_summary(df: pd.DataFrame) -> str:
    """Compact stdout table suitable for the CLI."""
    header = f"{'criterion':<14} {'our 20':>16} {'nasa 9':>16} {'delta':>8}  {'|t|':>6}"
    lines = [header, "-" * len(header)]
    for crit, row in df.iterrows():
        ours = f"{row['our_top_n_mean']:.3f} +- {row['our_top_n_std']:.3f}"
        theirs = f"{row['nasa_mean']:.3f} +- {row['nasa_std']:.3f}"
        delta = f"{row['delta']:+.3f}"
        t_value = row["t_statistic"]
        t_abs = abs(t_value) if t_value == t_value else float("nan")
        t_str = f"{t_abs:.2f}" if t_abs == t_abs else "n/a"
        lines.append(f"{str(crit):<14} {ours:>16} {theirs:>16} {delta:>8}  {t_str:>6}")
    return "\n".join(lines)
