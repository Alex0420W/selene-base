"""Tests for the v1.5 per-region tiled ranker.

End-to-end synthetic-data run: a small smooth elevation surface (so
slopes are well below the HLS 8° max), a uniform illumination, a
uniform aggregate score, and the v1.5 horizon-profile NPZ produced by
:func:`run_tiled_per_region`. The ranker should accept top-N sites
inside the polygon, with all four HLS filters satisfied at the
high-resolution grid.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

from selene_base.data.load import LUNAR_SOUTH_POLAR_CRS
from selene_base.pipeline import preprocess_tiled, rank_per_region_tiled
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

POLAR_CRS = str(LUNAR_SOUTH_POLAR_CRS)


def _build_synthetic_inputs(tmp_path: Path, *, resolution_m: float = 240.0) -> dict[str, Path]:
    """Materialise a 240 m global score+illumination grid plus a synthetic LOLA source.

    Returns a dict of paths the ranker expects.
    """
    raw_dir = tmp_path / "raw"
    processed = tmp_path / "processed"
    outputs = tmp_path / "outputs"
    (raw_dir / "lola").mkdir(parents=True)
    processed.mkdir(parents=True)
    outputs.mkdir(parents=True)

    polygons = regions_polygons_to_geodataframe(target_crs=POLAR_CRS)
    bx_min, by_min, bx_max, by_max = polygons.total_bounds
    pad_m = 110_000.0
    xmin = float(np.floor((bx_min - pad_m) / resolution_m) * resolution_m)
    ymin = float(np.floor((by_min - pad_m) / resolution_m) * resolution_m)
    xmax = float(np.ceil((bx_max + pad_m) / resolution_m) * resolution_m)
    ymax = float(np.ceil((by_max + pad_m) / resolution_m) * resolution_m)
    width = int(round((xmax - xmin) / resolution_m))
    height = int(round((ymax - ymin) / resolution_m))
    xs = np.linspace(xmin + resolution_m / 2, xmax - resolution_m / 2, width)
    ys = np.linspace(ymax - resolution_m / 2, ymin + resolution_m / 2, height)
    xx, yy = np.meshgrid(xs, ys)

    # Smooth elevation: gentle bowl, no steep features → slope filter passes.
    z = (xx**2 + yy**2) / 1e10  # 1e10 → max slope < 1° over the grid
    z = z.astype(np.float32)

    elev = xr.DataArray(
        z, dims=("y", "x"), coords={"y": ys, "x": xs}, name="elevation_m"
    ).rio.write_crs(POLAR_CRS, inplace=False)
    src_tif = raw_dir / "lola" / "synthetic.tif"
    elev.rio.to_raster(src_tif, driver="GTiff", compress="DEFLATE")

    # Illumination = 0.7 everywhere (passes >= 0.33 threshold).
    illum = xr.DataArray(
        np.full((height, width), 0.7, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": ys, "x": xs},
    ).rio.write_crs(POLAR_CRS, inplace=False)
    illum_path = processed / "illumination_southpole_240m.tif"
    illum.rio.to_raster(illum_path, driver="GTiff", compress="DEFLATE")

    # Aggregate score with a smooth peak, so the ranker has a preference order.
    score_arr = (0.5 + 0.3 * np.exp(-((xx**2 + yy**2) / (2 * 200_000.0**2)))).astype(np.float32)
    score = xr.DataArray(score_arr, dims=("y", "x"), coords={"y": ys, "x": xs}).rio.write_crs(
        POLAR_CRS, inplace=False
    )
    score_path = outputs / "score_southpole.tif"
    score.rio.to_raster(score_path, driver="GTiff", compress="DEFLATE")

    return {
        "raw": raw_dir,
        "processed": processed,
        "outputs": outputs,
        "lola_source": src_tif,
        "illum": illum_path,
        "score": score_path,
    }


def test_tiled_rank_runs_end_to_end_on_synthetic(tmp_path: Path) -> None:
    paths = _build_synthetic_inputs(tmp_path, resolution_m=240.0)

    # Step 1: horizon NPZ via the v1.5 preprocess driver.
    preprocess_tiled.run_tiled_per_region(
        resolution_m=240.0,
        region_codes=["SP"],
        processed_dir=paths["processed"],
        source_path=paths["lola_source"],
        use_gpu=False,
        echo=lambda _msg: None,
    )

    # Step 2: tiled ranker.
    # The synthetic landscape is intentionally smooth (so slopes pass
    # the 8° HLS filter), but a smooth surface at -88° S means Earth
    # is below the curvature-only horizon for most of the libration
    # cycle — real polar terrain provides the local relief that lifts
    # visibility above 0.5. Drop the DTE threshold for this test so
    # the rest of the pipeline can be exercised end-to-end.
    sites = rank_per_region_tiled.run(
        resolution_m=240.0,
        region_codes=["SP"],
        processed_dir=paths["processed"],
        outputs_dir=paths["outputs"],
        raw_dir=paths["raw"],
        score_map_path=paths["score"],
        source_path=paths["lola_source"],
        hls_dte_visibility_min=0.0,
        echo=lambda _msg: None,
    )

    assert len(sites) >= 1
    assert {"region_name", "region_code", "rank_in_region", "score", "lat", "lon"}.issubset(
        sites.columns
    )
    assert (sites["region_code"] == "SP").all()
    assert sites["resolution_m"].iloc[0] == 240.0
    assert sites["rank_in_region"].min() == 1

    # Output artefacts written.
    assert (paths["outputs"] / "per_region_tiled" / "sites.geojson").exists()
    assert (paths["outputs"] / "per_region_tiled" / "per_region_summary.json").exists()


def test_tiled_rank_handles_missing_horizon_npz(tmp_path: Path) -> None:
    paths = _build_synthetic_inputs(tmp_path, resolution_m=240.0)
    # Skip the preprocess step → no NPZ on disk.
    sites = rank_per_region_tiled.run(
        resolution_m=240.0,
        region_codes=["SP"],
        processed_dir=paths["processed"],
        outputs_dir=paths["outputs"],
        raw_dir=paths["raw"],
        score_map_path=paths["score"],
        source_path=paths["lola_source"],
        echo=lambda _msg: None,
    )
    # Empty result, no crash.
    assert len(sites) == 0
    assert (paths["outputs"] / "per_region_tiled" / "sites.geojson").exists()
