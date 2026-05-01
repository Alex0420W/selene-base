"""Per-site HTML reports.

For each top-N site we generate a self-contained HTML page with the
site's lat/lon, total score, per-criterion bar chart (embedded base64
PNG), a 50 km × 50 km mini-map cropped from the aggregate score COG,
the nearest NASA Artemis IV (formerly Artemis III) region with distance, and the dominant
criterion. Plus a top-level ``index.html`` linking every report so the
``data/outputs/sites/`` directory is browsable in any browser, offline.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401
from pyproj import Transformer

from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

POLAR_PROJ = (
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)
MINI_MAP_RADIUS_KM = 25.0  # → 50 km × 50 km crop


def _png_data_uri_from_fig(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _is_finite(v: object) -> bool:
    if v is None:
        return False
    try:
        f = float(v)
    except (TypeError, ValueError):
        return False
    return f == f  # excludes NaN


def _bar_chart_uri(row: gpd.GeoSeries) -> str:
    crit_cols = [c for c in row.index if c.startswith("score_")]
    names = [c.removeprefix("score_") for c in crit_cols]
    available = [_is_finite(row[c]) for c in crit_cols]
    values = [float(row[c]) if _is_finite(row[c]) else 0.0 for c in crit_cols]

    fig, ax = plt.subplots(figsize=(5.0, 2.6))
    bars = ax.barh(names, values, color="#06d6a0")
    for bar, ok in zip(bars, available, strict=True):
        if not ok:
            bar.set_color("#bbbbbb")
            bar.set_hatch("///")
    ax.set_xlim(0, 1)
    ax.set_xlabel("score")
    ax.invert_yaxis()
    for i, (v, ok) in enumerate(zip(values, available, strict=True)):
        label = f"{v:.2f}" if ok else "n/a"
        ax.text(min(v + 0.02, 0.95), i, label, va="center", fontsize=9)
    ax.set_title("per-criterion score (grey = data not yet wired)")
    fig.tight_layout()
    return _png_data_uri_from_fig(fig)


def _mini_map_uri(score_cog: Path, x_m: float, y_m: float) -> str:
    da = rioxarray.open_rasterio(score_cog, masked=True).squeeze("band", drop=True)
    transform = da.rio.transform()
    pixel_size = abs(float(transform.a))
    radius_pixels = int(round((MINI_MAP_RADIUS_KM * 1000.0) / pixel_size))
    col = int(round((x_m - transform.c) / transform.a - 0.5))
    row = int(round((y_m - transform.f) / transform.e - 0.5))
    h = da.sizes["y"]
    w = da.sizes["x"]
    r0 = max(0, row - radius_pixels)
    r1 = min(h, row + radius_pixels + 1)
    c0 = max(0, col - radius_pixels)
    c1 = min(w, col + radius_pixels + 1)
    crop = da.isel(y=slice(r0, r1), x=slice(c0, c1)).to_numpy()

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(crop, cmap="plasma", origin="upper", vmin=0, vmax=1)
    cy = max(0, min(crop.shape[0] - 1, row - r0))
    cx = max(0, min(crop.shape[1] - 1, col - c0))
    ax.scatter([cx], [cy], s=120, edgecolor="white", facecolor="cyan", linewidths=2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"±{int(MINI_MAP_RADIUS_KM)} km neighbourhood")
    return _png_data_uri_from_fig(fig)


def _nearest_nasa_region(
    site_lat: float, site_lon: float, nasa_regions: gpd.GeoDataFrame | None
) -> tuple[str, float] | None:
    if nasa_regions is None or len(nasa_regions) == 0:
        return None
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, POLAR_PROJ, always_xy=True)
    sx, sy = transformer.transform(site_lon, site_lat)
    nasa_geo = nasa_regions.to_crs(LUNAR_GEOGRAPHIC_CRS)
    rxs, rys = transformer.transform(
        nasa_geo["lon"].to_numpy(),
        nasa_geo["lat"].to_numpy(),
    )
    dists_km = np.hypot(np.asarray(rxs) - sx, np.asarray(rys) - sy) / 1000.0
    idx = int(np.argmin(dists_km))
    return str(nasa_geo["name"].iloc[idx]), float(dists_km[idx])


def _dominant_criteria(row: gpd.GeoSeries, top_n: int = 2) -> list[str]:
    sub = {
        c.removeprefix("score_"): float(row[c])
        for c in row.index
        if c.startswith("score_") and _is_finite(row[c])
    }
    if not sub:
        return []
    ranked = sorted(sub.items(), key=lambda kv: -kv[1])
    return [name for name, _ in ranked[:top_n]]


def generate_site_report(
    site_row: gpd.GeoSeries,
    score_cog: Path,
    output_dir: Path,
    *,
    nasa_regions: gpd.GeoDataFrame | None = None,
) -> Path:
    """Generate a per-site HTML report and return the path written.

    Args:
        site_row: One row of the ranked sites GeoDataFrame; must carry
            ``site_id``, ``rank``, ``lat``, ``lon``, ``score``,
            ``x_m``, ``y_m``, and any ``score_<criterion>`` columns.
        score_cog: Path to the aggregate score COG (used for the
            mini-map crop).
        output_dir: Where the HTML lands. Created if missing.
        nasa_regions: Optional NASA candidates GeoDataFrame; when
            present, the nearest one + distance is included.

    Returns:
        Path to the written HTML.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    site_id = str(site_row["site_id"])
    rank = int(site_row["rank"])
    out_path = output_dir / f"{site_id}.html"

    bar_uri = _bar_chart_uri(site_row)
    mini_uri = _mini_map_uri(score_cog, float(site_row["x_m"]), float(site_row["y_m"]))
    nearest = _nearest_nasa_region(float(site_row["lat"]), float(site_row["lon"]), nasa_regions)
    dominant = _dominant_criteria(site_row, top_n=3)

    nearest_html = (
        f"<p>Nearest NASA Artemis IV region: <b>{nearest[0]}</b> ({nearest[1]:.1f} km away)</p>"
        if nearest is not None
        else ""
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{site_id} · selene-base site report</title>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 760px;
          margin: 2em auto; padding: 0 1em; color: #222; }}
  header {{ border-bottom: 1px solid #ddd; padding-bottom: 0.6em; margin-bottom: 1em; }}
  h1 {{ margin: 0; font-size: 1.6em; }}
  .score {{ font-size: 1.4em; font-weight: bold; color: #06d6a0; }}
  .meta {{ color: #555; font-size: 0.95em; }}
  .panels {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1em; align-items: start; }}
  img {{ max-width: 100%; height: auto; }}
  footer {{ margin-top: 2em; color: #888; font-size: 0.85em; }}
  a {{ color: #1d4ed8; }}
</style>
</head>
<body>
<header>
  <h1>{site_id} · rank {rank}</h1>
  <p class="meta">lat {float(site_row["lat"]):.3f}°, lon {float(site_row["lon"]):.3f}°</p>
  <p class="score">aggregate score {float(site_row["score"]):.3f}</p>
</header>
<section class="panels">
  <div>
    <img src="{bar_uri}" alt="per-criterion scores">
  </div>
  <div>
    <img src="{mini_uri}" alt="local score map">
  </div>
</section>
{nearest_html}
<p>Dominant criteria: {", ".join(dominant) if dominant else "—"}</p>
<footer>
  <p><a href="index.html">← all sites</a></p>
  <p>Generated by selene-base · MIT licence · data © NASA / LRO</p>
</footer>
</body>
</html>
"""
    out_path.write_text(html, encoding="utf-8")
    return out_path


def generate_site_index(sites: gpd.GeoDataFrame, output_dir: Path) -> Path:
    """Write a top-level ``index.html`` that links every site report."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[str] = []
    for _, row in sites.iterrows():
        rows.append(
            f"<tr>"
            f"<td>{int(row['rank'])}</td>"
            f"<td><a href='{row['site_id']}.html'>{row['site_id']}</a></td>"
            f"<td>{float(row['lat']):.3f}</td>"
            f"<td>{float(row['lon']):.3f}</td>"
            f"<td>{float(row['score']):.3f}</td>"
            f"</tr>"
        )
    table_rows = "\n".join(rows)
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>selene-base · top sites</title>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 720px;
          margin: 2em auto; padding: 0 1em; color: #222; }}
  table {{ border-collapse: collapse; width: 100%; margin-top: 1em; }}
  th, td {{ text-align: left; padding: 0.4em 0.7em; border-bottom: 1px solid #eee; }}
  tr:hover td {{ background: #f5f5f5; }}
  th {{ font-weight: 600; }}
  a {{ color: #1d4ed8; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
</style>
</head>
<body>
<h1>selene-base · top {len(sites)} candidate sites</h1>
<p>Click a site to see its per-criterion breakdown and ±25 km neighbourhood.</p>
<table>
<thead><tr><th>rank</th><th>id</th><th>lat (°)</th><th>lon (°)</th><th>score</th></tr></thead>
<tbody>
{table_rows}
</tbody>
</table>
<p style="margin-top: 2em"><a href="../webmap.html">← interactive web map</a></p>
</body>
</html>
"""
    index_path = output_dir / "index.html"
    index_path.write_text(html, encoding="utf-8")
    return index_path
