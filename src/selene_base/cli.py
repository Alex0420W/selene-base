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
from selene_base.pipeline import compare_wueller as _compare_wueller
from selene_base.pipeline import coupling_sweep as _coupling_sweep
from selene_base.pipeline import preprocess as _preprocess
from selene_base.pipeline import preprocess_tiled as _preprocess_tiled
from selene_base.pipeline import rank as _rank
from selene_base.pipeline import rank_per_region as _rank_per_region
from selene_base.pipeline import rank_per_region_tiled as _rank_per_region_tiled
from selene_base.pipeline import score as _score
from selene_base.pipeline import score_wueller_sites as _score_wueller_sites
from selene_base.pipeline import sensitivity as _sensitivity
from selene_base.pipeline import validate as _validate
from selene_base.pipeline import validate_per_region as _validate_per_region
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
    resolution: str = typer.Option(
        "80m",
        "--resolution",
        help=(
            "LOLA-only: native grid resolution to fetch. One of '80m' "
            "(v1.4 default, ~115 MB) or '20m' (v1.5 high-res, ~1.85 GB). "
            "Ignored for other datasets."
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
    if dataset is Dataset.lola:
        res_m = _parse_resolution(resolution)
        if res_m not in _download.LOLA_RESOLUTIONS_M:
            raise typer.BadParameter(
                f"--resolution must be one of {_download.LOLA_RESOLUTIONS_M!r} "
                f"(got {resolution!r})",
                param_hint="--resolution",
            )
        path = _download.download_lola(
            _download.DEFAULT_RAW_DIR / dataset.value, resolution_m=res_m
        )
        typer.echo(f"  {dataset.value} ({res_m} m) ->{path}")
        return
    handler = _download.DATASETS[dataset.value]
    path = handler(_download.DEFAULT_RAW_DIR / dataset.value)
    typer.echo(f"  {dataset.value} ->{path}")


def _parse_resolution(value: str) -> int:
    """Parse a CLI resolution string like '20m', '20', or '20 m' to int metres."""
    s = str(value).strip().lower().rstrip("m").strip()
    try:
        return int(s)
    except ValueError as exc:
        raise typer.BadParameter(
            f"could not parse --resolution {value!r} as an integer in metres",
            param_hint="--resolution",
        ) from exc


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
    resolution: str | None = typer.Option(
        None,
        "--resolution",
        help=(
            "Override the analysis-grid resolution (in metres). Currently "
            "only consulted in --tiled-per-region mode; the global driver "
            "always uses the value from --region-config (240 m for v1.4)."
        ),
    ),
    tiled_per_region: bool = typer.Option(
        False,
        "--tiled-per-region",
        help=(
            "v1.5 mode: skip the global preprocess loop and instead run the "
            "Earth-LOS horizon-profile derivation per USGS region tile "
            "(polygon bbox + 100 km buffer) on the GPU. Tiles are cached "
            "as horizon_profile_southpole_<resolution>m_<region_code>.npz. "
            "Pair with --resolution 20 for a Wueller-class high-resolution run."
        ),
    ),
    region_code: list[str] | None = typer.Option(
        None,
        "--region-code",
        help=(
            "Tiled mode only: restrict the run to one or more USGS RegionCodes "
            "(e.g. 'SP' for Slater Plain). Repeatable. Default: all 9 regions."
        ),
    ),
    lola_source: Path | None = typer.Option(
        None,
        "--lola-source",
        help=(
            "Tiled mode only: explicit path to a LOLA elevation source file "
            "(.tif/.lbl/.img). Bypasses the PDS-naming auto-detect, enabling "
            "PGDA mosaics (e.g. ldem_83s_10mpp_adj.tif) and arbitrary "
            "user-supplied DEMs. When omitted, the finest-resolution source "
            "in data/raw/lola/ is picked automatically."
        ),
        exists=False,
        dir_okay=False,
    ),
) -> None:
    """Reproject every available raw raster onto the common 240 m grid.

    Idempotent: rasters already cached as
    ``data/processed/<name>_southpole_240m.tif`` are skipped unless
    ``--overwrite`` is set. Datasets whose raw bytes are missing are
    logged and skipped; Robbins (vector) is rasterised by the hazard
    criterion in week 3 and is not warped here.

    With ``--tiled-per-region`` (v1.5 mode), the global pipeline is
    bypassed and the GPU horizon-profile derivation is run per USGS
    polygon tile at ``--resolution`` metres.
    """
    if tiled_per_region:
        resolution_m = _parse_resolution(resolution) if resolution is not None else 20
        tiled_results = _preprocess_tiled.run_tiled_per_region(
            resolution_m=float(resolution_m),
            region_codes=region_code if region_code else None,
            processed_dir=processed_dir,
            overwrite=overwrite,
            source_path=lola_source,
        )
        typer.echo("")
        typer.echo(_preprocess_tiled.format_summary(tiled_results))
        return

    if resolution is not None:
        raise typer.BadParameter(
            "--resolution is only honoured in --tiled-per-region mode; "
            "the global driver reads resolution from --region-config.",
            param_hint="--resolution",
        )
    if region_code:
        raise typer.BadParameter(
            "--region-code is only honoured in --tiled-per-region mode.",
            param_hint="--region-code",
        )
    if lola_source is not None:
        raise typer.BadParameter(
            "--lola-source is only honoured in --tiled-per-region mode.",
            param_hint="--lola-source",
        )

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
    method: str = typer.Option(
        "weighted_sum",
        "--method",
        help=(
            "Aggregator to combine the per-criterion score grids. "
            "'weighted_sum' (default, v1.4+) is a renormalising weighted "
            "linear sum. 'topsis' (v1.7) is the Hwang & Yoon 1981 "
            "Technique for Order Preference by Similarity to Ideal "
            "Solution — same input + weights, but rewards balanced "
            "profiles over lopsided ones (a cell that saturates one "
            "criterion at the cost of zeroing another scores lower than "
            "an even-but-moderate cell with the same weighted sum)."
        ),
    ),
) -> None:
    """Run every available criterion and aggregate into one score map.

    Week 2 ships ``slope`` only; missing criteria trigger a warning and
    the remaining weights are renormalised. Output:
    ``<outputs_dir>/score_southpole.tif``.
    """
    if method not in ("weighted_sum", "topsis"):
        raise typer.BadParameter(
            f"--method must be one of 'weighted_sum' or 'topsis' (got {method!r})",
            param_hint="--method",
        )
    _score.run(
        weights_path=weights,
        region_config=region_config,
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        overwrite=overwrite,
        method=method,  # type: ignore[arg-type]
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


@app.command(name="rank-per-region")
def rank_per_region(
    n_per_region: int = typer.Option(
        10,
        "--n-per-region",
        min=1,
        help="Maximum number of HLS-compliant sites to extract per USGS region.",
    ),
    min_distance_km: float = typer.Option(
        2.0,
        "--min-distance-km",
        help="Minimum pairwise distance between sites within the same region (km).",
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
        help="Directory holding the slope, illumination, and LOS COGs.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory under which the per_region/ subdirectory is written.",
        file_okay=False,
    ),
    tiled_per_region: bool = typer.Option(
        False,
        "--tiled-per-region",
        help=(
            "v1.5 mode: rank within per-USGS-region tiles at high resolution. "
            "Reads horizon_profile_southpole_<resolution>m_<code>.npz produced "
            "by `selene preprocess --tiled-per-region`, derives 20 m slope and "
            "Earth-LOS visibility on the tile grid, applies HLS filters at the "
            "fine resolution, and writes sites to data/outputs/per_region_tiled/."
        ),
    ),
    resolution: str | None = typer.Option(
        None,
        "--resolution",
        help=(
            "Tiled mode only: analysis-grid resolution (m) of the per-tile "
            "horizon NPZ produced by `selene preprocess --tiled-per-region`. "
            "Default 20 m for v1.5."
        ),
    ),
    region_code: list[str] | None = typer.Option(
        None,
        "--region-code",
        help=(
            "Tiled mode only: restrict the run to one or more USGS RegionCodes "
            "(repeatable). Default: all 9 regions."
        ),
    ),
    method: str = typer.Option(
        "weighted_sum",
        "--method",
        help=(
            "Aggregator the consumed score map was produced with. "
            "'weighted_sum' (default) reads <outputs-dir>/score_southpole.tif "
            "and writes per_region(_tiled)/. 'topsis' (v1.7) reads "
            "<outputs-dir>/topsis/score_southpole.tif and writes "
            "per_region(_tiled)_topsis/, so a topsis catalog can sit "
            "side-by-side with the weighted-sum catalog without "
            "clobbering it. Run `selene score --method topsis "
            "--outputs-dir data/outputs/topsis` first to produce the "
            "topsis aggregate, then this command picks it up. Override "
            "explicitly with --score-map and --outputs-dir if needed."
        ),
    ),
) -> None:
    """Rank top-N HLS-compliant sites *within each USGS region*.

    Applies NASA's published HLS hard filters (slope ≤ 8°, 100 m
    buffer, illumination ≥ 33 %, DTE visibility ≥ 50 %) inside every
    USGS-published Artemis IV (formerly Artemis III) polygon, then ranks the survivors by
    aggregate score. Sites are guaranteed inside their named polygon
    by construction.

    With ``--tiled-per-region`` (v1.5 mode), the HLS filters and NMS run
    on a per-tile high-resolution grid using the v1.5 horizon profile.
    """
    if method not in ("weighted_sum", "topsis"):
        raise typer.BadParameter(
            f"--method must be one of 'weighted_sum' or 'topsis' (got {method!r})",
            param_hint="--method",
        )
    # When --method=topsis, default the score-map and outputs-subdir to
    # the topsis-suffixed location so a topsis catalog can sit alongside
    # the default weighted-sum catalog without clobbering it. Explicit
    # --score-map / --outputs-dir overrides win.
    effective_score_map = score_map
    effective_outputs_dir = outputs_dir
    per_region_subdir = "per_region_tiled" if tiled_per_region else "per_region"
    if method == "topsis":
        if score_map is None:
            effective_score_map = Path("data/outputs/topsis/score_southpole.tif")
        if outputs_dir == Path("data/outputs"):
            effective_outputs_dir = Path("data/outputs/topsis")
        per_region_subdir = "per_region_tiled_topsis" if tiled_per_region else "per_region_topsis"

    if tiled_per_region:
        resolution_m = _parse_resolution(resolution) if resolution is not None else 20
        _rank_per_region_tiled.run(
            resolution_m=float(resolution_m),
            region_codes=region_code if region_code else None,
            processed_dir=processed_dir,
            outputs_dir=effective_outputs_dir,
            score_map_path=effective_score_map,
            n_per_region=n_per_region,
            min_distance_km=min_distance_km,
            per_region_subdir=per_region_subdir,
        )
        return

    if resolution is not None:
        raise typer.BadParameter(
            "--resolution is only honoured in --tiled-per-region mode.",
            param_hint="--resolution",
        )
    if region_code:
        raise typer.BadParameter(
            "--region-code is only honoured in --tiled-per-region mode.",
            param_hint="--region-code",
        )
    _rank_per_region.run(
        score_map_path=effective_score_map,
        processed_dir=processed_dir,
        outputs_dir=effective_outputs_dir,
        n_per_region=n_per_region,
        min_distance_km=min_distance_km,
        per_region_subdir=per_region_subdir,
    )


@app.command(name="score-wueller-sites")
def score_wueller_sites(
    processed_dir: Path = typer.Option(
        Path("data/processed"),
        "--processed-dir",
        help="Directory holding the 240 m raw + scored COGs.",
        file_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help=(
            "Directory holding score_southpole.tif and where the "
            "v1.5_catalog_wueller_evaluation.csv + "
            "wueller_sites_scored_by_selene.json deliverables are written."
        ),
        file_okay=False,
    ),
    in_scope_only: bool = typer.Option(
        True,
        "--in-scope-only/--no-in-scope-only",
        help=(
            "When on (default), evaluate only the 73 Wueller sites whose "
            "region is in NASA's October 2024 down-selected nine. Off "
            "evaluates all 130 Wueller sites."
        ),
    ),
) -> None:
    """Score Wueller 2026 sites against selene's per-criterion rasters.

    Companion to ``selene compare-wueller``. The latter answers "are the
    two catalogs picking the same cells?" via spatial nearest-neighbour
    matching; this subcommand answers "do the two catalogs agree on
    what makes a *good* cell?" by sampling selene's six active 240 m
    criterion rasters plus the aggregate score at every Wueller-site
    coordinate.

    Per v1.5's design, per-criterion rasters are not materialised
    globally at 20 m. The 240 m criterion stack is what drives selene's
    aggregate score (which v1.5 upsamples bilinearly onto each
    per-region 20 m tile); the methodology converges between
    resolutions (median matched-pair distance 1.88 km at 240 m → 1.76
    km at 20 m), so 240 m criterion evaluation is representative.
    """
    _score_wueller_sites.run(
        processed_dir=processed_dir,
        outputs_dir=outputs_dir,
        in_scope_only=in_scope_only,
    )


@app.command(name="compare-wueller")
def compare_wueller(
    match_threshold_km: float = typer.Option(
        5.0,
        "--match-threshold-km",
        help=(
            "Distance below which a (selene, Wueller) pair is marked "
            "'matched'. Default 5 km — the upper end of typical regional "
            "candidate-site granularity. Not tuned from the agreement "
            "result."
        ),
    ),
    sites: Path | None = typer.Option(
        None,
        "--sites",
        help="selene per-region sites GeoJSON. Defaults to data/outputs/per_region/sites.geojson.",
        dir_okay=False,
    ),
    wueller_source: Path | None = typer.Option(
        None,
        "--wueller-source",
        "--wueller-csv",  # backward-compat alias for v1.4.0 callers
        help=(
            "Wueller 2026 site source. Defaults to the bundled real "
            "shapefile at src/selene_base/validation/data/wueller_2026/"
            "LandingSites.shp (130 sites, doi:10.5281/zenodo.17084058, "
            "CC-BY 4.0). Pass a .shp or .csv path to override; "
            "CSV mode is retained for backward compatibility with the "
            "v1.4.0 synthetic placeholder."
        ),
        dir_okay=False,
    ),
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory under which wueller_comparison.{json,csv} are written.",
        file_okay=False,
    ),
    filter_to_usgs_scope: bool = typer.Option(
        True,
        "--filter-to-usgs-scope/--no-filter-to-usgs-scope",
        help=(
            "When on (default), drop the ~57 Wueller sites whose region "
            "is not in NASA's October 2024 down-selected nine. selene-base "
            "only ranks within USGS-scope regions, so this is the apples-"
            "to-apples comparison. Pass --no-filter-to-usgs-scope to "
            "compare against all 130 Wueller sites."
        ),
    ),
) -> None:
    """Compare selene-base per-region sites against Wueller et al. 2026.

    Quantitative comparison against the 130-site Wueller 2026 catalog
    (doi:10.5281/zenodo.17084058, CC-BY 4.0). Default mode restricts
    to NASA's October 2024 down-selected nine regions; pass
    ``--no-filter-to-usgs-scope`` to compare the full 130-site catalog.
    """
    _compare_wueller.run(
        selene_sites_path=sites,
        wueller_source=wueller_source,
        outputs_dir=outputs_dir,
        match_threshold_km=match_threshold_km,
        filter_to_usgs_scope=filter_to_usgs_scope,
    )


@app.command(name="validate-per-region")
def validate_per_region(
    outputs_dir: Path = typer.Option(
        Path("data/outputs"),
        "--outputs-dir",
        help="Directory containing the per_region/ artefacts.",
        file_okay=False,
    ),
) -> None:
    """Summarise the per-region HLS-compliant ranking against USGS polygons.

    Reads ``data/outputs/per_region/sites.geojson`` and the cached
    summary JSON from ``selene rank-per-region``, prints a per-region
    table, and writes ``per_region_validation.json``.
    """
    _validate_per_region.run(outputs_dir=outputs_dir)


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
    """Compare ranked sites against NASA's Artemis IV (formerly Artemis III) candidate regions."""
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
