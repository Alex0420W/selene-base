# selene-base

> Multi-criteria habitat suitability for the lunar south pole, validated against NASA's nine announced Artemis III candidate landing regions.

[![CI](https://img.shields.io/badge/ci-pending-lightgrey)](.github/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)](pyproject.toml)
[![Status](https://img.shields.io/badge/status-v0.1-brightgreen)](#roadmap)

NASA's Artemis III mission will land humans near the lunar south pole around 2027. Selecting a base site there is a multi-criteria optimisation problem: the south pole is a maze of crater rims that catch grazing sunlight, deep permanently-shadowed cold-traps that may host water ice, and active thrust faults that re-localised Apollo-era shallow moonquakes have placed within tens of kilometres of candidate sites. **`selene-base`** fuses the modern LRO-era remote-sensing record (LOLA topography, Diviner thermal climatology, Mazarico illumination maps, LEND hydrogen abundance, the Robbins crater catalog, the Watters lobate-scarp catalog) with the historical Apollo seismic context to score every 240 m pixel of the polar cap and rank top candidate sites. The pipeline is end-to-end reproducible — `selene download && selene preprocess && selene score && selene rank && selene validate && selene viz` produces a ranked GeoJSON of sites, a per-site HTML report, and an interactive web map, on a developer laptop, in minutes, from public data.

## Headline result

Run on the three criteria with verified source data today (slope from LOLA, illumination from Mazarico, impact-hazard from the Robbins catalog rasterised at 3 km), with the remaining three criteria (thermal, ice, seismic) wired against the expected interfaces but skipping cleanly while their NASA archives remain TODO-flagged:

> **0 of selene-base's top 20 candidate sites land inside any of NASA's nine Artemis III candidate regions; 0 of 20 fall within 25 km of any centroid.** Two NASA regions — Slater Plain and de Gerlache Rim 2 — have a top-20 site between 25 and 30 km away; the rest are 50–157 km from the nearest top site. NASA's centroids score 0.20 to 0.63 on our aggregate map (max 0.97, 95th percentile 0.74) — comfortably below the 20th-best site at 0.88.

This is a **real disagreement**, not a bug, and the explanation is structural: NASA's selection emphasises proximity to permanently-shadowed water-ice deposits and Earth line-of-sight corridors. Both depend on criteria selene-base implements but cannot yet feed (`criteria/ice.py` needs a south-polar LEND product; `criteria/thermal.py` needs Diviner Tbol mosaics; `criteria/seismic.py` needs the Watters lobate-scarp catalog). The renormalised three-criterion regime that runs today instead favours flat, low-crater-density plains, which exist at every longitude — pulling the top-N off NASA's preferred meridian band. As the missing inputs land, the same pipeline will re-evaluate against the same nine NASA regions; the validation harness is in place.

![selene-base top 20 candidates vs NASA Artemis III regions](docs/img/webmap_screenshot.png)

| metric | value |
| --- | --- |
| top sites inside any NASA region (15 km disk) | 0 / 20 |
| top sites within 25 km of any NASA centroid | 0 / 20 |
| NASA regions with a top site within 30 km | 2 / 9 (Slater Plain, de Gerlache Rim 2) |
| NASA regions with a top site within 100 km | 8 / 9 (all but Cabeus B) |
| top-20 score range | 0.880–0.971 |
| score range across NASA centroids | 0.205 (Malapert Massif) – 0.629 (Cabeus B) |

For the per-region distance table see `data/outputs/validation.json`, or run `selene validate` on a fresh checkout to regenerate it. The interactive map lives at [`data/outputs/webmap.html`](data/outputs/webmap.html) after `selene viz`; per-site HTML reports under [`data/outputs/sites/`](data/outputs/sites/).

## Pipeline

```
data/raw/<dataset>/        ──load──▶  xr.DataArray (native CRS)
                              │
                              ▼ reproject_to_grid(target_crs, bounds, 240 m)
                              │
data/processed/<name>_southpole_240m.tif        (cached COG)
                              │
                              ▼ criterion.compute(...)            [six criteria]
                              │
data/processed/scored/<name>_score_southpole_240m.tif
                              │
                              ▼ scoring.aggregate.weighted_sum()  [renormalises]
                              │
data/outputs/score_southpole.tif                (final aggregate COG)
                              │
                              ▼ scoring.ranking.top_n_sites()     [NMS at 25 km]
                              │
data/outputs/top_sites.{geojson,csv}            (ranked sites + per-criterion sub-scores)
                              │
                              ▼ validation.comparison + viz
                              │
data/outputs/validation.json + webmap.html + sites/
```

### Quickstart

```bash
git clone https://github.com/Alex0420W/selene-base.git
cd selene-base
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e .

# Five-line clone-to-webmap path on the bundled ~12 MB sample dataset:
selene download --sample        # downloads + extracts data/raw/<sample>
selene preprocess               # warps + crater-density rasterisation -> data/processed/
selene score                    # six criteria; missing ones renormalise out cleanly
selene rank --top-n 20          # NMS + per-criterion sub-scores -> top_sites.{geojson,csv}
selene viz                      # webmap.html + per-site HTML reports

# Diagnostic & robustness:
selene validate                 # alignment metrics vs NASA's nine candidates
selene compare                  # per-criterion delta our top-20 vs NASA centroids
selene sensitivity --n-samples 200   # 200-sample weight-vector simplex sweep

# Full-resolution analysis (~290 MB raw, ~3 verified URLs):
selene download robbins
selene download lola
selene download illumination
# selene download diviner / lend / scarps remain TODO-flagged
selene preprocess && selene score && selene rank --top-n 20 --min-distance-km 25
```

`selene --help` lists every subcommand; `selene <cmd> --help` shows its options.

## Methodology

Every criterion produces a `[0, 1]` score grid where 1 is "best" and 0 is "unusable", aligned to the common 240 m south-polar stereographic grid (`+proj=stere +lat_0=-90 +lat_ts=-90 +R=1737400`, ±304 km, defined in [`config/region_southpole.yaml`](config/region_southpole.yaml)). Three normalisation primitives in [`scoring/normalize.py`](src/selene_base/scoring/normalize.py) — `min_max`, `optimal_range` (Gaussian), `inverse_threshold` — cover every criterion. The aggregate is a weighted linear sum that **renormalises across whichever criteria are present at score-time**, so a partial pipeline (today: slope, illumination, hazard) produces a comparable score to a complete one — only the absolute meaning of "0.97" shifts.

| Criterion | Score function | Source dataset | Resolution | Resampling | Default knobs |
| --- | --- | --- | --- | --- | --- |
| **Slope** | $s = \max(0,\,1-x/\theta_{\max})$ | LOLA LDEM 80 m (PDS3) | 80 m → 240 m | bilinear | $\theta_{\max} = 15°$ |
| **Illumination** | $s = \min(x/x_t,\,1)$ | Mazarico avgvisib 65°S 240 m | 240 m | bilinear | $x_t = 0.70$ |
| **Thermal** | $s = e^{-(\bar T - T^\star)^2/(2\sigma^2)} \cdot \max(0,\,1-(T_{\max}-T_{\min})/\Delta_{\max})$ | Diviner Tbol max/min (TODO) | ~240 m | bilinear | $T^\star=180\,$K, $\sigma=50\,$K, $\Delta_{\max}=200\,$K |
| **Ice** | $s = \mathrm{clip}(1-\mathrm{minmax}(\phi) + b\cdot\mathbb{1}[d(p, \mathcal{P})\le R],\,0,\,1)$ | LEND CSETN flux + PSR mask (TODO) | ~5–10 km | bilinear | $R=5\,$km, $b=0.3$ |
| **Hazard** | $s = \mathrm{clip}(1-d/d_{\mathrm{sat}},\,0,\,1)$ | Robbins 2018 catalog | vector → 240 m density | KDTree, 3 km radius | $d_{\mathrm{sat}}=50$ |
| **Seismic** | $s = \mathrm{clip}(\delta/\delta_{\mathrm{safe}},\,0,\,1)$ | Watters scarp catalog (TODO) | vector → 240 m distance | KDTree, 1 km densified vertices | $\delta_{\mathrm{safe}}=50\,$km |

Slope is computed at the 240 m target resolution from the already-downsampled LOLA DEM via `numpy.gradient` with explicit metric spacing (Zevenbergen & Thorne 1987 convention; ~5 % off Horn 1981 on smooth surfaces). Computing slope on the high-res 80 m DEM and then averaging slope-degrees double-smooths and biases low; computing on the target-resolution DEM keeps everything self-consistent.

The PSR mask used by the ice criterion is derived from the Mazarico illumination raster (`illumination < 0.001`); this gives the ice criterion something useful to do as soon as the LEND product is wired even if a more authoritative PSR catalog is not.

Default weights from [`config/weights_default.yaml`](config/weights_default.yaml): illumination 0.30, ice 0.25, slope 0.15, thermal 0.10, hazard 0.10, seismic 0.10. With only slope, illumination, and hazard available today, the renormalised effective weights are 0.27 (slope), 0.55 (illumination), 0.18 (hazard).

**Planned upgrade — TOPSIS.** A weighted linear sum lets a strong score on one criterion mask a near-disqualifying score on another (e.g. excellent illumination next to an active scarp). TOPSIS ranks each cell by its Euclidean distance to a synthetic "ideal" and "anti-ideal" point in criterion-score space, which penalises lop-sided profiles. It is on the roadmap as an alternate aggregator behind a `--method topsis` flag.

## Validation

`selene validate` compares the top-N ranked sites (from `data/outputs/top_sites.geojson`) against the disk-approximation polygons of NASA's nine announced Artemis III candidate regions in [`src/selene_base/validation/nasa_regions.py`](src/selene_base/validation/nasa_regions.py). Centroids are public information from NASA's October 2024 Artemis III site-selection announcement; we approximate each region as a 15 km disk around its centroid because NASA's actual polygons are not openly published in machine-readable form. **The disks are not authoritative geometry** — they're a defensible proximity proxy for this comparison.

Two metrics for each top site:

1. **Inside any region** — does the site fall inside any of the nine 15 km disks?
2. **Within X km of any centroid** — distance from the site to the nearest NASA centroid.

And two for each NASA region:

1. **Distance to nearest top-N site** — how far away is the closest selene-base candidate?
2. **Contains a top-N site** — is at least one selene-base candidate inside this region's disk?

### Per-region results (today)

![Distance from each NASA Artemis III candidate to the nearest selene-base top site](docs/img/validation_table.png)

| NASA candidate | nearest site | distance (km) | inside region? |
| --- | --- | ---: | --- |
| Cabeus B | site_16 | 157.3 | no |
| Haworth | site_08 | 52.1 | no |
| Malapert Massif | site_08 | 65.4 | no |
| Mons Mouton | site_16 | 98.6 | no |
| Mons Mouton Plateau | site_12 | 79.3 | no |
| Nobile Rim 1 | site_18 | 64.3 | no |
| Nobile Rim 2 | site_18 | 79.8 | no |
| de Gerlache Rim 2 | site_17 | 27.8 | no |
| Slater Plain | site_12 | 25.8 | no |

The honest read: with three of six criteria active and the operationally-dominant two (water-ice proximity, thermal stability) absent, selene-base ranks far-side and inter-crater plains highly because they're flat, low-crater-density, and well-illuminated, but they are not where NASA wants to land. The pipeline behaves correctly given the inputs it has; the result is informative about the limits of partial-criterion analysis. As Diviner / LEND / Watters land, this table updates.

## Robustness

Anyone reading 0/20 fairly asks: *is that just a function of the default weights?* Run `selene sensitivity --n-samples 200` to find out: it draws 200 weight vectors via Latin hypercube on the simplex over the criteria currently available, runs `aggregate → top_n_sites → proximity_analysis` for each, and reports the distribution of "NASA regions matched within 25 km" alongside the default-weight result.

![sensitivity over 200 weight samples](docs/img/sensitivity_distribution.png)

The result is bimodal at 0 and 2 region matches:

- **143 / 200 samples (71.5 %) match 0 regions within 25 km** — the modal outcome. The default-weights result is in this bucket.
- **56 / 200 samples (28 %) match 2 regions within 25 km** — Slater Plain (25.8 km) and de Gerlache Rim 2 (27.8 km), the two NASA candidates that already sit just beyond the 25 km threshold under the default weights.
- **0 / 200 samples match more than 2 regions.**

Achieving 2/9 region matches requires an extreme weight regime: the best sample uses `slope = 0.01, illumination = 0.04, hazard = 0.95` — essentially "rank by lowest crater density alone." That regime doesn't *find* additional NASA-aligned sites, it just renames the 25.8 km / 27.8 km near-misses as inside-threshold by virtue of weighting heavily a criterion that happens to score Slater and de Gerlache slightly higher than other plains. **No realistic three-criterion weight regime substantially improves alignment with NASA's selection.** This is the strongest statement we can make about the methodology: the 0/20 result is robust, and the disagreement is structural, not a function of weight choice.

## Diagnostic comparison

Run `selene compare` to ask a sharper question: *at NASA's centroids vs at our top-20, which criteria favour which set, by how much?* The output is a single table that explains the disagreement quantitatively.

![Per-criterion score: where we differ from NASA](docs/img/comparison.png)

| criterion | our top-20 | NASA 9 centroids | delta | \|t\| |
| --- | --- | --- | ---: | ---: |
| slope | 0.880 ± 0.089 | 0.285 ± 0.288 | +0.595 | 6.08 |
| illumination | 0.902 ± 0.063 | 0.321 ± 0.274 | +0.580 | 6.27 |
| hazard | 0.971 ± 0.031 | 0.969 ± 0.023 | +0.002 | 0.21 |

Three things to notice:

1. **Hazard is essentially identical.** Both NASA's nine and our top-20 sit in low-crater-density terrain (0.97 ± 0.03 vs 0.97 ± 0.02). On the criterion both site sets share, both site sets agree. **Our hazard layer is doing what NASA does.**
2. **Slope and illumination differ by ~0.6 in our favour.** NASA's regions span steep terrain (slope 0.285 ± 0.288 — note the high variance: Malapert Massif is on a literal massif) and accept low illumination (0.321 ± 0.274). Our pipeline treats both as disqualifying.
3. **The variance ratio matters as much as the means.** Our top-20 is tightly clustered (std ≈ 0.06–0.09) — the ranking is consistent. NASA's nine span the full range (std ≈ 0.27–0.29) — they're choosing site by site for reasons beyond what slope+illumination+hazard captures.

Reading the table together: NASA's selection accepts low slope and illumination scores in exchange for criteria selene-base does not yet consume — proximity to permanently-shadowed water-ice deposits and stable thermal regimes. When `criteria/ice.py` and `criteria/thermal.py` get their source data, the same `selene compare` rerun will produce different numbers; the diagnostic harness is in place.

The ``|t|`` column is a Welch two-sample t-statistic, reported informationally only. With ``n = 20`` against ``n = 9`` and structurally different sampling, a strict inferential frame is the wrong tool — the values are useful as a *ranking signal* (slope and illumination separate the site sets ~30× more strongly than hazard does), not as p-values.

## Architecture

```
selene-base/
├── src/selene_base/
│   ├── data/                # download + load + reproject + rasterize
│   ├── criteria/            # six [0,1] scoring functions
│   ├── scoring/             # normalize, aggregate (renormalising), ranking (NMS)
│   ├── validation/          # NASA candidate regions + proximity_analysis
│   ├── viz/                 # folium webmap + per-site HTML reports
│   ├── pipeline/            # one orchestrator module per CLI subcommand
│   └── cli.py               # typer CLI: download, preprocess, score, rank, validate, viz
├── config/                  # region_southpole.yaml, weights_default.yaml
├── data/                    # raw/ processed/ outputs/ (all gitignored)
├── notebooks/               # jupytext .py scripts; one per week
├── tests/                   # synthetic-data unit tests + skipif-guarded data tests
└── .github/workflows/ci.yml
```

The dependency graph is one-way: `data/` is the foundation; `criteria/` reads loaded rasters; `scoring/` aggregates criterion outputs; `validation/` and `viz/` consume scoring outputs; `pipeline/` orchestrates; `cli.py` exposes the orchestrators. Tests follow the same layering.

223 tests, ~80 % combined branch coverage, all running synthetically in CI on Python 3.11 and 3.12. Real-data tests are guarded with `pytest.mark.skipif(not Path(...).exists())` so the suite stays green without 200 MB of cached LRO data. CI runs a separate `pipeline-smoke` job on push to `main` that downloads the bundled ~12 MB sample tarball, runs `preprocess → score → rank → validate → compare`, and asserts every output file is on disk and schema-valid.

## Roadmap

- **Week 1 — data acquisition.** ✅ `selene download` for Robbins, LOLA, Mazarico illumination; Diviner / LEND / scarps URLs flagged.
- **Week 2 — common grid + slope criterion.** ✅ `reproject_to_grid`, COG cache, slope criterion end-to-end on real data.
- **Week 3 — full scoring + ranking.** ✅ All six criteria (3 on real data, 3 skip cleanly), KDTree crater density, NMS top-N extraction.
- **Week 4 — validation + visualisation.** ✅ NASA Artemis III proximity comparison, interactive folium web map, per-site HTML reports, validated v0.1.
- **Week 5 — robustness, diagnostic, sample data.** ✅ Latin-hypercube weight-sensitivity sweep, per-criterion `selene compare` diagnostic, bundled ~12 MB sample tarball, CI pipeline smoke test on the sample.
- **Future work.**
  - Resolve TODO-flagged URLs (Diviner Tbol max/min, LEND CSETN south-polar map, Watters scarp catalog) and rerun the validation against the now-six-criterion top-N.
  - TOPSIS aggregator behind `--method topsis`.
  - Earth line-of-sight criterion derived from LOLA elevation horizon checks.
  - ML-based criterion inputs (planned as a separate project, `selene-vision`).

## References

- Robbins, S. J. (2019). *A new global database of lunar impact craters >1–2 km: 1. Crater locations and sizes, comparisons with published databases, and global analysis.* Journal of Geophysical Research: Planets, 124, 871–892. [doi:10.1029/2018JE005592](https://doi.org/10.1029/2018JE005592)
- Mazarico, E., Neumann, G. A., Smith, D. E., Zuber, M. T., & Torrence, M. H. (2011). *Illumination conditions of the lunar polar regions using LOLA topography.* Icarus, 211(2), 1066–1081. [doi:10.1016/j.icarus.2010.10.030](https://doi.org/10.1016/j.icarus.2010.10.030)
- Smith, D. E., et al. (2010). *The Lunar Orbiter Laser Altimeter investigation on the Lunar Reconnaissance Orbiter mission.* Space Science Reviews, 150(1–4), 209–241. [doi:10.1007/s11214-009-9512-y](https://doi.org/10.1007/s11214-009-9512-y)
- Paige, D. A., et al. (2010). *The Lunar Reconnaissance Orbiter Diviner Lunar Radiometer Experiment.* Space Science Reviews, 150(1–4), 125–160. [doi:10.1007/s11214-009-9529-2](https://doi.org/10.1007/s11214-009-9529-2)
- Mitrofanov, I. G., et al. (2010). *Hydrogen mapping of the lunar south pole using the LRO Neutron Detector Experiment LEND.* Science, 330(6003), 483–486.
- Watters, T. R., Robinson, M. S., Banks, M. E., Tran, T., & Denevi, B. W. (2015). *Global thrust faulting on the Moon and the influence of tidal stresses.* Geology, 43(10), 851–854. [doi:10.1130/G37120.1](https://doi.org/10.1130/G37120.1)
- Civilini, F., Weber, R. C., Jiang, Z., Phillips, D., & Pan, W. (2023). *Constraints on the seismic hazard of young thrust faults on the Moon from re-located shallow moonquakes.* (Used as motivation for the seismic exclusion criterion.)
- NASA (October 2024). *Artemis III candidate landing regions.* [https://www.nasa.gov/feature/artemis-iii](https://www.nasa.gov/feature/artemis-iii)

## Notes for the reader

The interactive web map (`data/outputs/webmap.html` after `selene viz`) is built with folium / Leaflet and pulls Leaflet's JS and CSS from a CDN — so it needs an internet connection on first open. Everything else (the score raster, polygons, popups, per-site reports) is inlined and works offline. Per-site HTML reports under `data/outputs/sites/` are fully self-contained.

## License

MIT — see [LICENSE](LICENSE).
