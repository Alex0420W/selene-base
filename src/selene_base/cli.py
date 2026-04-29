"""Command-line entry point for selene-base.

Exposes a typer app whose subcommands trace the pipeline:
``download`` (week 1, real) ->``preprocess`` (week 2) ->``score`` (week 3)
→ ``rank`` (week 3) ->``viz`` (week 4). Subcommands beyond download still
raise NotImplementedError until their target week.
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path

import typer

from selene_base.data import download as _download
from selene_base.pipeline import compare as _compare
from selene_base.pipeline import coupling_sweep as _coupling_sweep
from selene_base.pipeline import preprocess as _preprocess
from selene_base.pipeline import rank as _rank
from selene_base.pipeline import score as _score
from selene_base.pipeline import sensitivity as _sensitivity
from selene_base.pipeline import validate as _validate
from selene_base.pipeline import viz as _viz_pipeline

app = typer.Typer(
    name="selene",
    help="Multi-criteria habitat suitability analyzer for the lunar south pole.",
    no_args_is_help=True,
    add_completion=False,
)


class Dataset(StrEnum):
    """Dataset choices for ``selene download``."""

    robbins = "robbins"
    lola = "lola"
    diviner = "diviner"
    illumination = "illumination"
    lend = "lend"
    scarps = "scarps"
    all = "all"


@app.command()
def download(
    dataset: Dataset = typer.Argument(
        Dataset.all,
        help="Dataset to fetch (or 'all' for every dataset).",
    ),
    sample: bool = typer.Option(
        False,
        "--sample",
        help=(
            "Fetch the bundled ~12 MB sample dataset instead of the full "
            "raw products. Lets the pipeline run end-to-end in seconds."
        ),
    ),
) -> None:
    """Fetch one or all source datasets to ``data/raw/``.

    Idempotent: rerunning skips any file already present and large
    enough to pass its sanity-check threshold.
    """
    if sample:
        _download.download_sample_data()
        typer.echo("")
        typer.echo("Next: selene preprocess && selene score && selene rank")
        return
    if dataset is Dataset.all:
        results = _download.download_all()
        for name, path in results.items():
            typer.echo(f"  {name:<14} ->{path}")
        return
    handler = _download.DATASETS[dataset.value]
    path = handler(_download.DEFAULT_RAW_DIR / dataset.value)
    typer.echo(f"  {dataset.value} ->{path}")


@app.command()
def preprocess(
    region_config: Path = typer.Option(
        Path("config/region_southpole.yaml"),
        "--region-config",
        help="Path to the south-polar grid definition YAML.",
        dir_okay=False,
    ),
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory to write cached COGs into.",
        file_okay=False,
    ),
    overwrite: bool = typer.Option(
        False,
        "--overwrite",
        help="Re-cache COGs even if they already exist on disk.",
    ),
    skip_los: bool = typer.Option(
        False,
        "--skip-los",
        help=(
            "Skip the slow Earth line-of-sight horizon-profile derivation "
            "(~5–15 min on full-resolution LOLA). Useful for smoke tests "
            "and runs that don't need the los_to_earth criterion."
        ),
    ),
) -> None:
    """Reproject every available raw raster onto the common 240 m grid.

    Idempotent: rasters already cached as
    ``data/processed/<name>_southpole_240m.tif`` are skipped unless
    ``--overwrite`` is set. Datasets whose raw bytes are missing are
    logged and skipped; Robbins (vector) is rasterised by the hazard
    criterion in week 3 and is not warped here.
    """
    results = _preprocess.run(
        region_config=region_config,
        processed_dir=processed_dir,
        overwrite=overwrite,
        compute_los=not skip_los,
    )
    typer.echo("")
    typer.echo(_preprocess.format_summary(results))


@app.command()
def score(
    weights: Path = typer.Option(
        Path("config/weights_default.yaml"),
        "--weights",
        "-w",
        help="Path to a YAML file mapping criterion name to non-negative weight.",
        exists=False,
        dir_okay=False,
    ),
    region_config: Path = typer.Option(
        Path("config/region_southpole.yaml"),
        "--region-config",
        help="Path to the south-polar grid definition YAML.",
        dir_okay=False,
    ),
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory holding cached input + per-criterion score COGs.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write the aggregate score COG into.",
        file_okay=False,
    ),
    overwrite: bool = typer.Option(
        False,
        "--overwrite",
        help="Re-compute per-criterion score COGs even if cached.",
    ),
) -> None:
    """Run every available criterion and aggregate into one score map.

    Week 2 ships ``slope`` only; missing criteria trigger a warning and
    the remaining weights are renormalised. Output:
    ``<outputs_dir>/score_southpole.tif``.
    """
    _score.run(
        weights_path=weights,
        region_config=region_config,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        overwrite=overwrite,
    )


@app.command()
def rank(
    top_n: int = typer.Option(20, "--top-n", "-n", min=1, help="Number of sites to extract."),
    min_distance_km: float = typer.Option(
        5.0,
        "--min-distance-km",
        help="Minimum pairwise distance between returned sites, in kilometres.",
    ),
    min_score: float = typer.Option(
        0.5,
        "--min-score",
        min=0.0,
        max=1.0,
        help="Floor on candidate score; cells below this never enter the running.",
    ),
    score_map: Path | None = typer.Option(
        None,
        "--score-map",
        help="Aggregate score COG. Defaults to <outputs-dir>/score_southpole.tif.",
        dir_okay=False,
    ),
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory holding per-criterion score COGs under scored/.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write top_sites.geojson and top_sites.csv into.",
        file_okay=False,
    ),
) -> None:
    """Extract the top-N geographically-distinct candidate sites via NMS."""
    _rank.run(
        score_map_path=score_map,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        top_n=top_n,
        min_distance_km=min_distance_km,
        min_score=min_score,
    )


@app.command()
def viz(
    score_map: Path | None = typer.Option(
        None,
        "--score-map",
        help="Aggregate score COG. Defaults to <outputs-dir>/score_southpole.tif.",
        dir_okay=False,
    ),
    sites_path: Path | None = typer.Option(
        None,
        "--sites",
        help="Top sites GeoJSON. Defaults to <outputs-dir>/top_sites.geojson.",
        dir_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write webmap.html and per-site reports into.",
        file_okay=False,
    ),
    skip_reports: bool = typer.Option(
        False,
        "--skip-reports",
        help="Build only the web map, skip per-site HTML reports.",
    ),
) -> None:
    """Generate the folium web map and per-site HTML reports."""
    _viz_pipeline.run(
        score_map_path=score_map,
        sites_path=sites_path,
        outputs_dir=outputs_dir,
        skip_reports=skip_reports,
    )


@app.command()
def validate(
    sites_path: Path | None = typer.Option(
        None,
        "--sites",
        help="Top sites GeoJSON. Defaults to <outputs-dir>/top_sites.geojson.",
        dir_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write validation.json into.",
        file_okay=False,
    ),
    near_km: float = typer.Option(
        25.0,
        "--near-km",
        help="Distance threshold (km) for the 'within X km of any centroid' metric.",
    ),
) -> None:
    """Compare ranked sites against NASA's Artemis III candidate regions."""
    _validate.run(
        sites_path=sites_path,
        outputs_dir=outputs_dir,
        near_km=near_km,
    )


@app.command()
def compare(
    sites_path: Path | None = typer.Option(
        None,
        "--sites",
        help="Top sites GeoJSON. Defaults to <outputs-dir>/top_sites.geojson.",
        dir_okay=False,
    ),
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory holding per-criterion score COGs under scored/.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write comparison.json and comparison.png into.",
        file_okay=False,
    ),
) -> None:
    """Per-criterion diagnostic: where do our top sites differ from NASA's?"""
    _compare.run(
        sites_path=sites_path,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
    )


@app.command()
def sensitivity(
    n_samples: int = typer.Option(
        200,
        "--n-samples",
        min=2,
        help="Number of weight vectors to draw from the simplex.",
    ),
    top_n: int = typer.Option(20, "--top-n", min=1, help="Sites to extract per sample."),
    min_distance_km: float = typer.Option(
        25.0,
        "--min-distance-km",
        help="NMS minimum pairwise separation, in kilometres.",
    ),
    near_km: float = typer.Option(
        25.0,
        "--near-km",
        help="Proximity threshold for the headline 'regions matched' metric.",
    ),
    far_km: float = typer.Option(
        100.0,
        "--far-km",
        help="Wider proximity threshold reported alongside (e.g. 100 km).",
    ),
    seed: int = typer.Option(42, "--seed", help="PRNG seed for the sweep."),
    weights: Path = typer.Option(
        Path("config/weights_default.yaml"),
        "--weights",
        "-w",
        help="Default-weights YAML used to mark the histogram baseline.",
        dir_okay=False,
    ),
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory holding per-criterion score COGs under scored/.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write sensitivity_results.parquet + .png into.",
        file_okay=False,
    ),
) -> None:
    """Sweep weight vectors and report robustness of the validation result."""
    _sensitivity.run(
        n_samples=n_samples,
        top_n=top_n,
        min_distance_km=min_distance_km,
        proximity_threshold_km=near_km,
        far_threshold_km=far_km,
        seed=seed,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        weights_path=weights,
    )


@app.command(name="coupling-sweep")
def coupling_sweep(
    distances_km: list[float] = typer.Option(
        [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0],
        "--distance-km",
        help="Coupling-distance values (km) to sweep over. Repeatable.",
    ),
    top_n: int = typer.Option(20, "--top-n", min=1, help="Sites per sweep."),
    min_distance_km: float = typer.Option(
        25.0,
        "--min-distance-km",
        help="NMS minimum pairwise separation, in kilometres.",
    ),
    near_km: float = typer.Option(
        25.0,
        "--near-km",
        help="Proximity threshold for the 'regions matched' metric.",
    ),
    weights: Path = typer.Option(
        Path("config/weights_default.yaml"),
        "--weights",
        "-w",
        help="Weights YAML; the coupling weight is honoured.",
        dir_okay=False,
    ),
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory holding cached input + per-criterion score COGs.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory to write coupling_sweep.parquet + .png into.",
        file_okay=False,
    ),
) -> None:
    """Sweep the spatial-coupling distance and report NASA-region alignment."""
    _coupling_sweep.run(
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        weights_path=weights,
        distances_km=list(distances_km),
        top_n=top_n,
        min_distance_km=min_distance_km,
        proximity_threshold_km=near_km,
    )


if __name__ == "__main__":
    app()
