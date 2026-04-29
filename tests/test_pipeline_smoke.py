"""End-to-end smoke test of preprocess → score on synthetic input.

Builds a tiny tilted-plane DEM in a temp directory, runs the pipeline
through both ``selene_base.pipeline.preprocess.run`` and the typer CLI,
then asserts the aggregate score COG was written and contains finite
values. No downloaded data required, runs in CI.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr
import yaml
from pyproj import CRS
from typer.testing import CliRunner

from selene_base.cli import app
from selene_base.data.reproject import is_cog
from selene_base.pipeline import preprocess as _preprocess
from selene_base.pipeline import score as _score

LUNAR_SOUTH_POLAR = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)


@pytest.fixture
def synthetic_lola(tmp_path: Path) -> Path:
    """Write a 1024×1024 tilted-plane DEM as a GeoTIFF on the polar grid."""
    raw_dir = tmp_path / "raw" / "lola"
    raw_dir.mkdir(parents=True, exist_ok=True)
    out = raw_dir / "synthetic_dem.tif"

    pixel_m = 240.0
    h = w = 1024
    grad = math.tan(math.radians(8.0))  # 8° tilt — under our 15° cutoff
    x = np.arange(w) * pixel_m
    z = np.broadcast_to(x * grad, (h, w)).astype(np.float32)

    da = xr.DataArray(
        z,
        dims=("y", "x"),
        coords={
            "y": np.linspace(150_000, -150_000, h),
            "x": np.linspace(-150_000, 150_000, w),
        },
        name="elevation_m",
    ).rio.write_crs(LUNAR_SOUTH_POLAR, inplace=False)
    da.rio.to_raster(out, driver="GTiff", compress="DEFLATE")
    return out


@pytest.fixture
def region_yaml(tmp_path: Path) -> Path:
    cfg = {
        "name": "south_pole_test",
        "crs": str(LUNAR_SOUTH_POLAR.to_proj4()),
        "resolution_m": 240,
        "bounds_m": [-150_000, -150_000, 150_000, 150_000],
        "lat_max": -80.0,
    }
    p = tmp_path / "region.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return p


@pytest.fixture
def weights_yaml(tmp_path: Path) -> Path:
    weights = {
        "illumination": 0.30,
        "ice": 0.25,
        "slope": 0.15,
        "thermal": 0.10,
        "hazard": 0.10,
        "seismic": 0.10,
    }
    p = tmp_path / "weights.yaml"
    p.write_text(yaml.safe_dump(weights))
    return p


def _spec_for(synthetic_dem: Path) -> tuple[_preprocess.DatasetSpec, ...]:
    def loader() -> xr.DataArray:
        return rioxarray.open_rasterio(synthetic_dem, masked=True).squeeze("band", drop=True)

    return (
        _preprocess.DatasetSpec(
            name="lola",
            raw_check=synthetic_dem,
            loader=loader,
            resampling="bilinear",
            note="synthetic test DEM",
        ),
    )


def test_pipeline_run_end_to_end(
    tmp_path: Path,
    synthetic_lola: Path,
    region_yaml: Path,
    weights_yaml: Path,
) -> None:
    processed = tmp_path / "processed"
    outputs = tmp_path / "outputs"

    results = _preprocess.run(
        region_config=region_yaml,
        processed_dir=processed,
        datasets=_spec_for(synthetic_lola),
        echo=lambda _msg: None,
    )
    assert any(r.status == "cached" and r.name == "lola" for r in results)
    cog = processed / "lola_southpole_240m.tif"
    assert cog.exists() and is_cog(cog)

    summary = _score.run(
        weights_path=weights_yaml,
        region_config=region_yaml,
        processed_dir=processed,
        outputs_dir=outputs,
        echo=lambda _msg: None,
    )
    assert summary.output_path.exists()
    assert summary.n_finite > 0
    assert 0.0 <= summary.minimum <= summary.maximum <= 1.0
    # 8° tilt under a 15° threshold should give non-trivial positive scores
    assert summary.mean > 0.0


def test_cli_preprocess_and_score_exit_zero(
    tmp_path: Path,
    synthetic_lola: Path,
    region_yaml: Path,
    weights_yaml: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(_preprocess, "RASTER_DATASETS", _spec_for(synthetic_lola))
    processed = tmp_path / "processed"
    outputs = tmp_path / "outputs"

    runner = CliRunner()
    pre = runner.invoke(
        app,
        [
            "preprocess",
            "--region-config",
            str(region_yaml),
            "--processed-dir",
            str(processed),
        ],
    )
    assert pre.exit_code == 0, pre.output
    assert "lola" in pre.output
    assert (processed / "lola_southpole_240m.tif").exists()

    sc = runner.invoke(
        app,
        [
            "score",
            "--weights",
            str(weights_yaml),
            "--region-config",
            str(region_yaml),
            "--processed-dir",
            str(processed),
            "--outputs-dir",
            str(outputs),
        ],
    )
    assert sc.exit_code == 0, sc.output
    assert (outputs / "score_southpole.tif").exists()


def test_score_without_preprocessed_inputs_raises(
    tmp_path: Path,
    region_yaml: Path,
    weights_yaml: Path,
) -> None:
    with pytest.raises(RuntimeError, match="no criterion score grids"):
        _score.run(
            weights_path=weights_yaml,
            region_config=region_yaml,
            processed_dir=tmp_path / "empty_processed",
            outputs_dir=tmp_path / "outputs",
            echo=lambda _msg: None,
        )


def test_pipeline_rank_runs_after_score(
    tmp_path: Path,
    synthetic_lola: Path,
    region_yaml: Path,
    weights_yaml: Path,
) -> None:
    from selene_base.pipeline import rank as _rank

    processed = tmp_path / "processed"
    outputs = tmp_path / "outputs"

    _preprocess.run(
        region_config=region_yaml,
        processed_dir=processed,
        datasets=_spec_for(synthetic_lola),
        echo=lambda _msg: None,
    )
    _score.run(
        weights_path=weights_yaml,
        region_config=region_yaml,
        processed_dir=processed,
        outputs_dir=outputs,
        echo=lambda _msg: None,
    )
    sites = _rank.run(
        processed_dir=processed,
        outputs_dir=outputs,
        top_n=5,
        min_distance_km=1.0,
        min_score=0.0,  # synthetic plane stays low; let any cell qualify
        echo=lambda _msg: None,
    )

    assert (outputs / "top_sites.geojson").exists()
    assert (outputs / "top_sites.csv").exists()
    assert len(sites) > 0
    for col in ("site_id", "rank", "score", "lat", "lon", "score_slope"):
        assert col in sites.columns


def test_cli_rank_exit_zero(
    tmp_path: Path,
    synthetic_lola: Path,
    region_yaml: Path,
    weights_yaml: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(_preprocess, "RASTER_DATASETS", _spec_for(synthetic_lola))
    processed = tmp_path / "processed"
    outputs = tmp_path / "outputs"

    runner = CliRunner()
    pre = runner.invoke(
        app,
        [
            "preprocess",
            "--region-config",
            str(region_yaml),
            "--processed-dir",
            str(processed),
        ],
    )
    assert pre.exit_code == 0, pre.output
    sc = runner.invoke(
        app,
        [
            "score",
            "--weights",
            str(weights_yaml),
            "--region-config",
            str(region_yaml),
            "--processed-dir",
            str(processed),
            "--outputs-dir",
            str(outputs),
        ],
    )
    assert sc.exit_code == 0, sc.output
    rk = runner.invoke(
        app,
        [
            "rank",
            "--top-n",
            "5",
            "--min-distance-km",
            "1.0",
            "--min-score",
            "0.0",
            "--processed-dir",
            str(processed),
            "--outputs-dir",
            str(outputs),
        ],
    )
    assert rk.exit_code == 0, rk.output
    assert (outputs / "top_sites.geojson").exists()
