"""Weight-vector sensitivity analysis.

Sweeps the criterion weight simplex via Latin hypercube sampling, runs
``aggregate → top_n_sites → proximity_analysis`` for each sample, and
returns a tidy DataFrame with one row per sample plus the headline
proximity metrics. The output answers the natural question that
follows the validation result: *is this robust to the weight choice?*

Implementation notes:

- Score COGs are loaded once and passed in as ``xr.DataArray``. The
  loop never touches disk again.
- Aggregation is the cheap step (linear combination of pre-computed
  rasters); top-N + proximity are the bottleneck. We accept the cost
  and keep the loop straightforward.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from selene_base.scoring.aggregate import weighted_sum
from selene_base.scoring.ranking import top_n_sites
from selene_base.validation.comparison import proximity_analysis

DEFAULT_NEAR_KM = 25.0
DEFAULT_FAR_KM = 100.0


def latin_hypercube_weights(
    n_samples: int,
    criterion_names: list[str],
    *,
    seed: int = 42,
) -> np.ndarray:
    """Generate ``n_samples`` weight vectors uniformly over the simplex.

    Uses the classic order-statistics construction: draw
    ``n_criteria - 1`` quantile cuts from a Latin hypercube, sort each
    row, prepend 0 / append 1, and take the differences. The resulting
    rows sum to 1 and are uniformly distributed on the simplex.

    Args:
        n_samples: Number of weight vectors to produce. Must be > 0.
        criterion_names: Order of weights in each row (only the length
            is read here; the names are echoed back by ``sweep_weights``).
        seed: PRNG seed for the underlying scipy QMC generator.

    Returns:
        Array of shape ``(n_samples, len(criterion_names))``; each row
        sums to 1, every entry is in ``[0, 1]``.

    Raises:
        ValueError: If ``n_samples <= 0`` or ``criterion_names`` is empty.
    """
    if n_samples <= 0:
        raise ValueError(f"n_samples must be positive, got {n_samples!r}")
    n_dim = len(criterion_names)
    if n_dim == 0:
        raise ValueError("criterion_names must be non-empty")
    if n_dim == 1:
        return np.ones((n_samples, 1), dtype=np.float64)

    from scipy.stats import qmc

    sampler = qmc.LatinHypercube(d=n_dim - 1, seed=seed)
    cuts = sampler.random(n=n_samples)  # shape (n_samples, n_dim - 1)
    cuts = np.sort(cuts, axis=1)
    padded = np.concatenate(
        [np.zeros((n_samples, 1)), cuts, np.ones((n_samples, 1))],
        axis=1,
    )
    weights = np.diff(padded, axis=1)
    return weights


def sweep_weights(
    score_maps: Mapping[str, xr.DataArray],
    weight_samples: np.ndarray,
    nasa_regions: gpd.GeoDataFrame,
    *,
    top_n: int = 20,
    min_distance_km: float = 25.0,
    proximity_threshold_km: float = DEFAULT_NEAR_KM,
    far_threshold_km: float = DEFAULT_FAR_KM,
    min_score: float = 0.0,
    criterion_order: list[str] | None = None,
    nasa_regions_polygons: gpd.GeoDataFrame | None = None,
) -> pd.DataFrame:
    """Run aggregate → top-N → proximity for each weight sample.

    Args:
        score_maps: Mapping from criterion name to its [0, 1] score
            DataArray on the common grid. Must be non-empty.
        weight_samples: 2-D array of shape ``(n_samples, n_criteria)``;
            each row a weight vector. Column order is given by
            ``criterion_order`` (or ``sorted(score_maps)`` when ``None``).
        nasa_regions: NASA candidate-region GeoDataFrame from
            :func:`selene_base.validation.nasa_regions.regions_to_geodataframe`.
        top_n: Number of sites the NMS extracts each iteration.
        min_distance_km: NMS minimum pairwise separation, in km.
        proximity_threshold_km: Distance threshold for the
            ``n_within_proximity_km`` column.
        far_threshold_km: Wider distance threshold reported alongside.
        min_score: Floor on candidate site score. Default 0.0 because
            different weight regimes shift the score distribution and
            a fixed 0.5 floor would silently discard sweeps where every
            cell scores below 0.5.
        criterion_order: Override the column → criterion mapping.

    Returns:
        DataFrame with one row per sample. Columns:

        - ``w_<criterion>`` for each criterion.
        - ``n_inside_region`` --top sites inside any NASA disk.
        - ``n_within_proximity_km`` --top sites within
          ``proximity_threshold_km`` of any centroid.
        - ``n_within_far_km`` --top sites within ``far_threshold_km`` of
          any centroid.
        - ``n_regions_with_top_site`` --NASA regions containing ≥1 top.
        - ``n_regions_within_proximity_km`` --NASA regions whose nearest
          top site is within the proximity threshold.
        - ``mean_score_at_nasa_centroids`` --mean aggregate score sampled
          at the NASA centroids.
        - ``mean_top_n_score`` --mean aggregate score of the top-N sites.

    Raises:
        ValueError: On shape mismatch or empty inputs.
    """
    if not score_maps:
        raise ValueError("score_maps must be non-empty")
    names = list(criterion_order) if criterion_order is not None else sorted(score_maps.keys())
    if any(name not in score_maps for name in names):
        missing = [name for name in names if name not in score_maps]
        raise ValueError(f"criterion_order references missing maps: {missing}")
    n_dim = len(names)
    if weight_samples.ndim != 2 or weight_samples.shape[1] != n_dim:
        raise ValueError(f"weight_samples must be (n_samples, {n_dim}); got {weight_samples.shape}")

    rows: list[dict[str, float]] = []
    centroid_xs: np.ndarray | None = None
    centroid_ys: np.ndarray | None = None

    # Pre-compute the NASA centroid pixel indices once --every sample
    # samples the same locations.
    sample_map = score_maps[names[0]]
    transform = sample_map.rio.transform()
    if transform is None:
        raise ValueError("score_maps lack a transform; cannot sample NASA centroids")

    from pyproj import Transformer

    from selene_base.data.load import LUNAR_GEOGRAPHIC_CRS

    nasa_geo = nasa_regions.to_crs(LUNAR_GEOGRAPHIC_CRS)
    transformer = Transformer.from_crs(LUNAR_GEOGRAPHIC_CRS, sample_map.rio.crs, always_xy=True)
    centroid_xs, centroid_ys = transformer.transform(
        nasa_geo["lon"].to_numpy(),
        nasa_geo["lat"].to_numpy(),
    )
    centroid_cols = np.round((centroid_xs - transform.c) / transform.a - 0.5).astype(int)
    centroid_rows = np.round((centroid_ys - transform.f) / transform.e - 0.5).astype(int)
    h, w = sample_map.sizes["y"], sample_map.sizes["x"]
    in_grid = (
        (centroid_rows >= 0) & (centroid_rows < h) & (centroid_cols >= 0) & (centroid_cols < w)
    )

    for sample_idx in range(weight_samples.shape[0]):
        weights = {name: float(weight_samples[sample_idx, k]) for k, name in enumerate(names)}
        # weighted_sum warns when there are extra weights beyond scores;
        # not relevant here (we always sweep over all available criteria),
        # but silence the warning anyway in case a caller passes a
        # subset.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            aggregate = weighted_sum(score_maps, weights)
        agg_arr = aggregate.to_numpy()

        sites = top_n_sites(
            aggregate,
            n=top_n,
            min_distance_m=min_distance_km * 1000.0,
            min_score=min_score,
        )
        if len(sites) > 0:
            prox = proximity_analysis(
                sites,
                nasa_regions,
                near_km=proximity_threshold_km,
                nasa_regions_polygons=nasa_regions_polygons,
            )
            far_prox = proximity_analysis(sites, nasa_regions, near_km=far_threshold_km)
            # The dict field name carries "25km" historically; the analysis
            # is parameterised by near_km so it really means "within
            # whichever threshold we passed in".
            n_within_far = far_prox["sites_within_25km_of_region"]
            n_regions_within = int(
                np.sum(
                    np.array([r["distance_km"] for r in prox["per_region"]])
                    <= proximity_threshold_km
                )
            )
            mean_top_n = float(sites["score"].mean())
        else:
            prox = {
                "sites_within_any_region": 0,
                "sites_within_25km_of_region": 0,
                "regions_with_a_top_site": 0,
            }
            n_within_far = 0
            n_regions_within = 0
            mean_top_n = float("nan")

        if in_grid.any():
            valid_rows = centroid_rows[in_grid]
            valid_cols = centroid_cols[in_grid]
            scores_at_centroids = agg_arr[valid_rows, valid_cols]
            mean_centroid_score = float(np.nanmean(scores_at_centroids))
        else:
            mean_centroid_score = float("nan")

        row: dict[str, float] = {f"w_{name}": weights[name] for name in names}
        row.update(
            {
                "n_inside_region": int(prox["sites_within_any_region"]),
                "n_within_proximity_km": int(prox["sites_within_25km_of_region"]),
                "n_within_far_km": int(n_within_far),
                "n_regions_with_top_site": int(prox["regions_with_a_top_site"]),
                "n_regions_within_proximity_km": int(n_regions_within),
                "mean_score_at_nasa_centroids": mean_centroid_score,
                "mean_top_n_score": mean_top_n,
                "n_inside_usgs_polygon": int(prox.get("sites_inside_any_usgs_polygon", 0)),
                "n_usgs_regions_with_top_site": int(
                    prox.get("regions_with_top_site_inside_usgs_polygon", 0)
                ),
                "median_distance_to_nearest_usgs_polygon_km": float(
                    prox.get("median_distance_to_nearest_usgs_polygon_km", float("nan"))
                ),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def sweep_coupling_distance(
    score_maps_without_coupling: Mapping[str, xr.DataArray],
    distance_to_psr: xr.DataArray,
    distance_to_ridge: xr.DataArray,
    weights: Mapping[str, float],
    nasa_regions: gpd.GeoDataFrame,
    *,
    distances_km: list[float] | None = None,
    top_n: int = 20,
    min_distance_km: float = 25.0,
    proximity_threshold_km: float = 25.0,
    far_threshold_km: float = 100.0,
    min_score: float = 0.0,
) -> pd.DataFrame:
    """Sweep ``coupling_distance_km``; recompute coupling and aggregate each.

    Args:
        score_maps_without_coupling: The five non-coupling per-criterion
            score grids (slope, illumination, hazard, thermal, ice).
        distance_to_psr: Pre-computed distance grid (metres) from
            :func:`selene_base.criteria.coupling.derive_distance_to_psr`.
        distance_to_ridge: Pre-computed distance grid (metres) from
            :func:`selene_base.criteria.coupling.derive_distance_to_sunlit_ridge`.
        weights: Weight vector (must include ``coupling`` plus every
            entry in ``score_maps_without_coupling``).
        nasa_regions: NASA Artemis III GeoDataFrame from
            :func:`selene_base.validation.nasa_regions.regions_to_geodataframe`.
        distances_km: List of coupling-distance values to test.
            Defaults to ``[1, 2, 3, 5, 7, 10, 15, 20]`` km.
        top_n: NMS site count.
        min_distance_km: NMS minimum pairwise separation, in km.
        proximity_threshold_km: Distance for the headline "regions
            matched" metric.
        far_threshold_km: Wider distance reported alongside.
        min_score: Floor on candidate site score.

    Returns:
        DataFrame with one row per ``coupling_distance_km`` tested:

        - ``coupling_distance_km``
        - ``n_inside_region``, ``n_within_proximity_km``,
          ``n_within_far_km``
        - ``n_regions_with_top_site``,
          ``n_regions_within_proximity_km``
        - ``mean_score_at_nasa_centroids``, ``mean_top_n_score``
    """
    if distances_km is None:
        distances_km = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

    from selene_base.criteria import coupling as coupling_criterion

    rows: list[dict[str, float]] = []
    for d_km in distances_km:
        coupling_score = coupling_criterion.compute(
            distance_to_psr,
            distance_to_ridge,
            coupling_distance_km=d_km,
        )
        score_maps = dict(score_maps_without_coupling)
        score_maps["coupling"] = coupling_score
        names = sorted(score_maps.keys())
        weight_vec = np.array([float(weights[n]) for n in names], dtype=np.float64)
        weight_vec /= weight_vec.sum()
        weight_samples = weight_vec.reshape(1, -1)
        result = sweep_weights(
            score_maps,
            weight_samples,
            nasa_regions,
            top_n=top_n,
            min_distance_km=min_distance_km,
            proximity_threshold_km=proximity_threshold_km,
            far_threshold_km=far_threshold_km,
            min_score=min_score,
            criterion_order=names,
        )
        row = {"coupling_distance_km": float(d_km)}
        for col in (
            "n_inside_region",
            "n_within_proximity_km",
            "n_within_far_km",
            "n_regions_with_top_site",
            "n_regions_within_proximity_km",
            "mean_score_at_nasa_centroids",
            "mean_top_n_score",
        ):
            row[col] = float(result[col].iloc[0])
        rows.append(row)
    return pd.DataFrame(rows)


def best_weights(results: pd.DataFrame, criterion_names: list[str]) -> dict[str, float]:
    """Pick the sample with the most NASA regions matched within proximity.

    Ties are broken by the next-most-discriminating column (number of
    sites inside any region, then mean score at NASA centroids).
    """
    if results.empty:
        raise ValueError("results DataFrame is empty")
    sort_cols = [
        "n_regions_within_proximity_km",
        "n_regions_with_top_site",
        "n_inside_region",
        "mean_score_at_nasa_centroids",
    ]
    ranked = results.sort_values(sort_cols, ascending=False)
    best = ranked.iloc[0]
    return {name: float(best[f"w_{name}"]) for name in criterion_names}


def save_results(results: pd.DataFrame, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(path)
    return path


def render_summary(
    results: pd.DataFrame,
    *,
    default_weights: dict[str, float] | None = None,
    proximity_threshold_km: float = DEFAULT_NEAR_KM,
) -> str:
    """Compact stdout block summarising the sweep distribution."""
    n = len(results)
    lines = [f"sensitivity sweep: {n} weight samples"]
    target_col = "n_regions_within_proximity_km"
    metric = results[target_col]
    lines.append(
        f"  regions matched within {proximity_threshold_km:.0f} km --"
        f"min/median/max: {int(metric.min())} / {int(metric.median())} / {int(metric.max())}"
    )
    sites_col = results["n_within_proximity_km"]
    s_min, s_med, s_max = int(sites_col.min()), int(sites_col.median()), int(sites_col.max())
    lines.append(
        f"  top sites within {proximity_threshold_km:.0f} km --"
        f"min/median/max: {s_min} / {s_med} / {s_max}"
    )
    if default_weights is not None:
        diff = pd.Series(0.0, index=results.index)
        for k, v in default_weights.items():
            col = f"w_{k}"
            if col in results.columns:
                diff = diff + (results[col] - v).abs()
        nearest_idx = diff.idxmin()
        nearest = results.loc[nearest_idx]
        lines.append(
            f"  nearest sample to default weights --regions matched: {int(nearest[target_col])}"
        )
    return "\n".join(lines)
