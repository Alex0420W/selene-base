# selene-base

A multi-criteria habitat suitability engine for the lunar south pole. `selene-base`
fuses LRO-era datasets — LOLA topography, Diviner thermal climatology, Mazarico
illumination maps, LEND hydrogen abundance, the Robbins crater catalog, and the
Watters lobate-scarp catalog — into a single ranked list of candidate Artemis base
sites. It is unusual in that it grounds the modern remote-sensing pipeline against
the historical Apollo seismic record (re-localised to active scarps by Civilini et
al., 2023), so a site's "no-go" zones reflect not just slope and shadow but also
where the Moon is still tectonically alive.

[![CI](https://img.shields.io/badge/ci-pending-lightgrey)](.github/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)](pyproject.toml)
[![Status](https://img.shields.io/badge/status-active%20development-yellow)](#roadmap)

> **Status.** Week 2 complete. Reprojection pipeline and slope criterion
> shipping. Beginning week 3: remaining five criteria.
>
> Today the pipeline can: download Robbins / LOLA / Mazarico illumination,
> warp every available raster onto the common 240 m south-polar grid,
> derive slope from LOLA via central differences, score it with the
> slope criterion, and aggregate via the (renormalising) weighted sum
> in `scoring/aggregate.py`. Diviner / LEND URLs remain TODO-flagged in
> [`data/download.py`](src/selene_base/data/download.py).

## Pipeline (today)

```
data/raw/<dataset>/        ──load──▶  xr.DataArray (native CRS)
                              │
                              ▼ reproject_to_grid(target_crs, bounds, 240 m)
                              │
data/processed/<name>_southpole_240m.tif   (cached COG)
                              │
                              ▼ criterion.compute(...)  [week 2: slope; week 3: rest]
                              │
data/processed/scored/<name>_score_southpole_240m.tif   (per-criterion COG)
                              │
                              ▼ scoring.aggregate.weighted_sum(...)
                              │
data/outputs/score_southpole.tif   (final aggregate COG)
```

## Architecture

```
                ┌───────────────┐
   raw rasters  │ data.download │   (LOLA, Diviner, illumination, LEND,
   + catalogs   │  data.load    │    Robbins craters, Watters scarps)
                └──────┬────────┘
                       │
                       ▼
                ┌───────────────┐
                │ data.reproject│   south-polar stereographic, 240 m,
                │  → 240 m grid │   bounds ±304 km, lat ≤ -80°
                └──────┬────────┘
                       │
        ┌──────────────┼──────────────┬──────────────┬──────────────┐
        ▼              ▼              ▼              ▼              ▼
   ┌─────────┐   ┌──────────┐   ┌─────────┐   ┌────────┐    ┌──────────┐
   │ slope   │   │ illum.   │   │ thermal │   │  ice   │    │ hazard   │
   │criterion│   │criterion │   │criterion│   │criterion│    │ + seismic│
   └────┬────┘   └────┬─────┘   └────┬────┘   └───┬────┘    └────┬─────┘
        └────────────┴───────────────┴────────────┴──────────────┘
                                  │
                                  ▼
                        ┌───────────────────┐
                        │ scoring.aggregate │   weighted_sum →
                        │  (weighted sum)   │   suitability map ∈ [0,1]
                        └─────────┬─────────┘
                                  │
                                  ▼
                        ┌───────────────────┐
                        │ scoring.ranking   │   non-maximum suppression,
                        │  (top-N + NMS)    │   min separation 25 km
                        └─────────┬─────────┘
                                  │
                                  ▼
                        ┌───────────────────┐
                        │ viz.webmap        │   folium HTML +
                        │ viz.site_report   │   per-site HTML reports
                        └───────────────────┘
```

## Quickstart

```bash
git clone https://github.com/<you>/selene-base.git
cd selene-base
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e ".[dev]"

# Week 1 — data acquisition
selene download robbins         # ~92 MB raw, ~400 KB filtered south-polar slice
selene download lola            # ~115 MB south-polar DEM (PDS3 IMG + LBL)
selene download illumination    # ~82 MB Mazarico avgvisib raster
selene download diviner         # URL TODO — verify before running
selene download lend            # URL TODO — verify before running
selene download all             # convenience: every dataset in turn (idempotent)
python notebooks/01_data_inventory.py   # writes sanity plots to data/outputs/sanity/

# Week 2 — reproject + slope criterion (real today)
selene preprocess                                  # warps every available raster to 240 m COGs
selene score --weights config/weights_default.yaml # week 2 ships slope only; warns on missing criteria
python notebooks/02_slope_first_pass.py            # elevation / slope / score side-by-side

# Week 3+ — not yet implemented (still NotImplementedError)
selene rank --top-n 20
selene viz
```

`selene --help` lists every subcommand; `selene <cmd> --help` shows its options.
Every download is idempotent: rerunning skips files already on disk that
pass their minimum-size sanity check.

## Data sources

| Dataset | Product | Resolution | URL status | Role in scoring |
| --- | --- | --- | --- | --- |
| LOLA south-polar DEM | `ldem_80s_80m.{img,lbl}` (PDS3) | 80 m / pixel | [verified](https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/polar/img/) | slope, elevation context |
| Mazarico illumination | `avgvisib_65s_240m_201608.{img,lbl}` | 240 m / pixel | [verified](https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/extras/illumination/release_2016/img/) | illumination fraction |
| Robbins crater catalog | PDS4 bundle CSV (Robbins 2018) | vector, Ø ≥ 1 km | [verified](https://astrogeology.usgs.gov/search/map/Moon/Research/Craters/lunar_crater_database_robbins_2018) | impact / ejecta hazard |
| Diviner Tbol max/min | UCLA `level4_polar/` mosaics | ~240 m / pixel | **TODO** — filename unverified | thermal stability |
| LEND epithermal flux | polar neutron-flux map | ~5–10 km / pixel | **TODO** — product not located | water-ice proxy |
| Watters lobate scarps | catalog (Watters et al. 2015) | vector | not yet wired | seismic exclusion |

**TODO URLs** are flagged in [`src/selene_base/data/download.py`](src/selene_base/data/download.py)
with the original PDS / archive starting points; verify before running
`selene download diviner` or `selene download lend`. The Watters scarp
download is still on the backlog and will land alongside the seismic
criterion in week 3.

All datasets are reprojected onto a single south-polar stereographic grid
(`+proj=stere +lat_0=-90 +lat_ts=-90 +R=1737400`) at 240 m / pixel over
±304 km, defined in [`config/region_southpole.yaml`](config/region_southpole.yaml).

### Resampling choices per dataset

| Dataset | Resampling | Why |
| --- | --- | --- |
| LOLA elevation | bilinear | smooth continuous surface; bilinear is standard |
| Illumination | bilinear | continuous percentage; bilinear preserves the dynamic range |
| Diviner Tmax/Tmin (when wired) | bilinear | continuous Tbol field |
| LEND (when wired) | bilinear | already coarse, smoothing OK at 240 m |
| Robbins | n/a | vector, rasterised by the hazard criterion in week 3 |

The slope criterion derives its gradient from the **already-downsampled
240 m DEM**, not from the native 80 m LOLA DEM averaged after-the-fact.
Computing slope on the high-res DEM and then averaging slope-degrees
double-smooths and biases towards lower values; computing slope on the
target-resolution DEM keeps the result self-consistent with the rest of
the analysis grid.

## Scoring methodology

Each criterion produces a [0, 1] score grid, where 1 is best and 0 is unusable.
Three normalisation primitives in [`scoring/normalize.py`](src/selene_base/scoring/normalize.py)
cover every criterion. The aggregate
[`scoring/aggregate.py`](src/selene_base/scoring/aggregate.py) tolerates
missing criteria: weights for criteria that aren't yet implemented are
silently dropped (with a warning) and the remaining weights are
renormalised to sum to 1, so weeks 2 and 3 can ship a partial pipeline
without rebalancing the weights file.

**Slope** (lower is better, hard cutoff at 15°):

$$ s_\text{slope}(x) = \max\!\left(0,\; 1 - \frac{x}{\theta_\text{max}}\right),\quad \theta_\text{max} = 15° $$

Slope itself is computed via :func:`numpy.gradient` with explicit metric
spacing (Zevenbergen & Thorne 1987 convention), which is the most common
choice in planetary GIS and gives values within ~5% of the Sobel-weighted
Horn (1981) kernel on smooth surfaces. Edge pixels and any cell whose
3×3 stencil touched a NaN are explicitly NaN.

**Illumination** (linear in average sunlit fraction):

$$ s_\text{illum}(x) = \mathrm{clip}\!\left(\frac{x - x_\text{lo}}{x_\text{hi} - x_\text{lo}},\; 0,\; 1\right) $$

**Thermal** (Gaussian peak at a moderate target temperature):

$$ s_\text{therm}(x) = \exp\!\left(-\frac{(x - T^\star)^2}{2\sigma^2}\right),\quad T^\star \approx 230\,\mathrm{K} $$

**Ice / hydrogen** (linear in inferred water-equivalent wt%):

$$ s_\text{ice}(x) = \mathrm{clip}\!\left(\frac{x}{x_\text{ref}},\; 0,\; 1\right) $$

**Hazard** (distance from craters with $D \ge D_\text{min}$):

$$ s_\text{haz}(p) = \mathrm{clip}\!\left(\frac{d(p,\,\text{craters})}{r_\text{eject}},\; 0,\; 1\right) $$

**Seismic** (distance from active lobate scarps):

$$ s_\text{seis}(p) = \mathrm{clip}\!\left(\frac{d(p,\,\text{scarps})}{d_\text{safe}},\; 0,\; 1\right),\quad d_\text{safe} \approx 50\,\mathrm{km} $$

The aggregate suitability is a weighted linear sum:

$$ S(p) = \sum_{c} w_c \cdot s_c(p),\qquad \sum_c w_c = 1 $$

Default weights from [`config/weights_default.yaml`](config/weights_default.yaml):

| Criterion | Weight |
| --- | ---: |
| illumination | 0.30 |
| ice | 0.25 |
| slope | 0.15 |
| thermal | 0.10 |
| hazard | 0.10 |
| seismic | 0.10 |

**Planned upgrade — TOPSIS.** A weighted linear sum lets a strong score on one
criterion mask a near-disqualifying score on another (e.g. excellent illumination
right next to an active scarp). TOPSIS — Technique for Order of Preference by
Similarity to Ideal Solution — ranks each cell by its Euclidean distance to a
synthetic "ideal" and "anti-ideal" point in criterion-score space, which
penalises lop-sided profiles. It is on the roadmap as an alternate aggregator
that callers can opt into via `--method topsis`.

## Validation plan

Top-ranked sites are compared against NASA's nine announced Artemis III candidate
regions: Cabeus B, Haworth, Malapert Massif, Mons Mouton, Mons Mouton Plateau,
Nobile Rim 1, Nobile Rim 2, de Gerlache Rim 2, and Slater Plain. Two metrics:

1. **Coverage** — fraction of the nine NASA candidates that fall within
   `min_distance_km` of any top-N selene-base site.
2. **Per-region rank** — for each NASA candidate, the rank of the closest
   selene-base site.

Results table (filled in week 4):

| NASA candidate | Closest selene-base rank | Distance (km) | Aggregate score |
| --- | ---: | ---: | ---: |
| Cabeus B | \<results pending\> | \<results pending\> | \<results pending\> |
| Haworth | \<results pending\> | \<results pending\> | \<results pending\> |
| Malapert Massif | \<results pending\> | \<results pending\> | \<results pending\> |
| Mons Mouton | \<results pending\> | \<results pending\> | \<results pending\> |
| Mons Mouton Plateau | \<results pending\> | \<results pending\> | \<results pending\> |
| Nobile Rim 1 | \<results pending\> | \<results pending\> | \<results pending\> |
| Nobile Rim 2 | \<results pending\> | \<results pending\> | \<results pending\> |
| de Gerlache Rim 2 | \<results pending\> | \<results pending\> | \<results pending\> |
| Slater Plain | \<results pending\> | \<results pending\> | \<results pending\> |

A high-coverage outcome would be evidence that the chosen criteria and weights
are at least *consistent* with NASA's site-selection process; sensitivity sweeps
over the weights file are planned to characterise how robust the agreement is.

## Project structure

```
selene-base/
├── README.md
├── pyproject.toml
├── LICENSE
├── .python-version
├── .gitignore
├── config/
│   ├── region_southpole.yaml
│   └── weights_default.yaml
├── data/
│   ├── raw/         (gitignored — populated by `selene download`)
│   ├── processed/   (gitignored — populated by `selene preprocess`)
│   └── outputs/     (gitignored — populated by `selene viz`)
├── src/selene_base/
│   ├── cli.py
│   ├── data/        (download, load, reproject)
│   ├── criteria/    (slope, illumination, thermal, ice, hazard, seismic)
│   ├── scoring/     (normalize, aggregate, ranking)
│   └── viz/         (webmap, site_report)
├── notebooks/
├── tests/
└── .github/workflows/ci.yml
```

## Roadmap

A four-week plan, with each module's docstring tagged to its target week.

- **Week 1 — data ingestion.** Implement `data.download.*`, `data.load.*`,
  and a CLI `selene download` that pulls every product into `data/raw/`. The
  three normalisation primitives in `scoring.normalize` are already real and
  tested (they have no upstream dependency), as is `scoring.aggregate.weighted_sum`.
- **Week 2 — common grid + first criterion.** Implement
  `data.reproject.reproject_to_grid` and the slope criterion against the LOLA
  DEM. End-to-end `selene preprocess` should run on the real raw data.
- **Week 3 — full scoring + ranking.** Implement the remaining five criteria
  (illumination, thermal, ice, hazard, seismic) and `scoring.ranking.top_n_sites`
  with non-maximum suppression. `selene score` and `selene rank` work end-to-end.
- **Week 4 — visualisation + validation.** Implement `viz.webmap.build_webmap`
  and `viz.site_report.render_site_report`. Run the NASA Artemis III candidate
  comparison, fill in the results table, and freeze v0.1.

## References

- Robbins, S. J. (2019). *A new global database of lunar impact craters >1–2 km:
  1. Crater locations and sizes, comparisons with published databases, and global
  analysis.* Journal of Geophysical Research: Planets, 124, 871–892.
- Mazarico, E., Neumann, G. A., Smith, D. E., Zuber, M. T., & Torrence, M. H.
  (2011). *Illumination conditions of the lunar polar regions using LOLA
  topography.* Icarus, 211(2), 1066–1081.
- Watters, T. R., Robinson, M. S., Banks, M. E., Tran, T., & Denevi, B. W.
  (2015). *Global thrust faulting on the Moon and the influence of tidal
  stresses.* (Watters lobate-scarp catalog used here.)
- Civilini, F., Weber, R. C., Jiang, Z., Phillips, D., & Pan, W. (2023).
  *Constraints on the seismic hazard of young thrust faults on the Moon from
  re-located shallow moonquakes.* (Used for the seismic exclusion criterion.)
- NASA (2022). *Artemis III candidate landing regions* — public press materials
  identifying the nine candidate regions used for validation.

## License

MIT — see [LICENSE](LICENSE).
