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
    all = "all"


@app.command()
def download(
    dataset: Dataset = typer.Argument(
        Dataset.all,
        help="Dataset to fetch (or 'all' for every dataset).",
    ),
) -> None:
    """Fetch one or all source datasets to ``data/raw/``.

    Idempotent: rerunning skips any file already present and large
    enough to pass its sanity-check threshold.
    """
    if dataset is Dataset.all:
        results = _download.download_all()
        for name, path in results.items():
            typer.echo(f"  {name:<14} ->{path}")
        return
    handler = _download.DATASETS[dataset.value]
    path = handler(_download.DEFAULT_RAW_DIR / dataset.value)
    typer.echo(f"  {dataset.value} ->{path}")


@app.command()
def preprocess() -> None:
    """Reproject every raw raster onto the common south-polar grid.

    Reads ``config/region_southpole.yaml`` for the target CRS, bounds,
    and resolution, then writes Cloud-Optimized GeoTIFFs into
    ``data/processed/``.

    Filled in week 2.
    """
    raise NotImplementedError("filled in week 2")


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
) -> None:
    """Run every criterion and aggregate into a single suitability map.

    Each criterion in :mod:`selene_base.criteria` produces a [0, 1] score
    grid; :func:`selene_base.scoring.aggregate.weighted_sum` combines them.

    Filled in week 3.
    """
    raise NotImplementedError("filled in week 3")


@app.command()
def rank(
    top_n: int = typer.Option(20, "--top-n", "-n", min=1, help="Number of sites to extract."),
    min_distance_km: float = typer.Option(
        25.0,
        "--min-distance-km",
        help="Minimum pairwise distance between returned sites, in kilometres.",
    ),
) -> None:
    """Extract the top-N geographically-distinct candidate sites via NMS.

    Filled in week 3.
    """
    raise NotImplementedError("filled in week 3")


@app.command()
def viz() -> None:
    """Generate the folium web map and per-site HTML reports.

    Writes ``data/outputs/webmap.html`` plus one report per ranked site.

    Filled in week 4.
    """
    raise NotImplementedError("filled in week 4")


if __name__ == "__main__":
    app()
