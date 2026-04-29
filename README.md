# selene-base

> Multi-criteria habitat suitability for the lunar south pole, validated against NASA's nine announced Artemis III candidate landing regions.

[![CI](https://img.shields.io/badge/ci-pending-lightgrey)](.github/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)](pyproject.toml)
[![Status](https://img.shields.io/badge/status-v0.1-brightgreen)](#roadmap)

NASA's Artemis III mission will land humans near the lunar south pole around 2027. Selecting a base site there is a multi-criteria optimisation problem: the south pole is a maze of crater rims that catch grazing sunlight, deep permanently-shadowed cold-traps that may host water ice, and active thrust faults that re-localised Apollo-era shallow moonquakes have placed within tens of kilometres of candidate sites. **`selene-base`** fuses the modern LRO-era remote-sensing record (LOLA topography, Diviner thermal climatology, Mazarico illumination maps, LEND hydrogen abundance, the Robbins crater catalog, the Watters lobate-scarp catalog) with the historical Apollo seismic context to score every 240 m pixel of the polar cap and rank top candidate sites. The pipeline is end-to-end reproducible — `selene download && selene preprocess && selene score && selene rank && selene validate && selene viz` produces a ranked GeoJSON of sites, a per-site HTML report, and an interactive web map, on a developer laptop, in minutes, from public data.

## Headline finding

> **A six-week engineering arc — identify a methodology limit, hypothesise a structural fix, implement and measure it.** Adding two well-validated Diviner-PRP-derived criteria (thermal, ice) to a 3-criterion baseline *worsened* NASA-region alignment (median distance 65 km → 102 km) because weighted-sum MCDA cannot model NASA's coupled spatial constraint. We then implemented a derived **spatial-coupling criterion** that scores cells by joint proximity to a PSR *and* a sunlit ridge — a multiplicative product that encodes the AND. Result: median distance 102 km → 115 km on default weights (the new criterion makes the aggregate harder), but the *200-sample sensitivity ceiling* improved from 2/9 → 3/9 NASA regions matched, and the criterion's diagnostic shows it correctly identifies the rim-band geometry NASA targets. The fix is real but partial. The remaining gap is in the validation metric, not the criterion.

This is the project's strongest result: a methodological story arc with measurements at each step, including the negative ones.

![Coupling score with NASA candidates and our top-20 overlaid](docs/img/coupling_overlay.png)

### Three-stage validation history

| Validation metric | 3-criteria | 5-criteria (+ Diviner PRP) | 6-criteria (+ spatial coupling) |
| --- | --- | --- | --- |
| top sites inside any NASA region (15 km disk) | 0 / 20 | 0 / 20 | 0 / 20 |
| top sites within 25 km of any centroid | 0 / 20 | 0 / 20 | 0 / 20 |
| closest NASA region | Slater Plain @ 25.8 km | de Gerlache Rim 2 @ 64.8 km | de Gerlache Rim 2 @ **47.8 km** |
| 2nd closest | de Gerlache Rim 2 @ 27.8 km | Slater Plain @ 71.6 km | Slater Plain @ **55.4 km** |
| median NASA region distance | 65 km | 102 km | 115 km |
| 200-sample sensitivity, *best* regions matched within 25 km | 2 / 9 | 2 / 9 | **3 / 9** |
| 200-sample sensitivity, modal outcome | 0/9 (71.5 %) | 0/9 (92.5 %) | 0/9 (87 %) |
| coupling score at NASA centroids (mean ± std) | n/a | n/a | **0.000 ± 0.000** |
| coupling score at our top-20 (mean ± std) | n/a | n/a | **0.053 ± 0.165** |

### What the coupling criterion did and did not fix

Two things shifted in the right direction:

1. **The closest-distance metric improved by a third.** The closest NASA region went from 64.8 km (5-crit) to 47.8 km (6-crit) — adding the coupling criterion pulled the highest-ranked sites back toward the south polar rim where NASA's candidates cluster. site_01 at (-89.7°, +17.7°) is now within 50 km of de Gerlache Rim 2.
2. **The sensitivity ceiling lifted from 2/9 → 3/9.** The 5-criterion weight simplex could not produce more than 2 NASA-region matches at any weight vector. The 6-criterion simplex *can* produce 3 matches — and the best-found regime weights coupling at 0.27, not zero. The new criterion is doing real work in the search space, not being ignored.

Two things did not shift:

3. **The headline 0/20 number did not move.** The minimum gap to a NASA centroid is still 47.8 km vs the 25 km threshold — closer, but not inside.
4. **NASA centroids score 0.000 on coupling.** This is the most informative single number in the whole project. NASA's nine centroids land *inside* their candidate disks — usually inside a crater, on a massif, or in a plain — not on the rim band that maximises the coupling score. The criterion is correctly identifying the rim geometry NASA *targets within* each region, but the validation metric measures distance to the *centroid* of each region, which sits 5–15 km off-rim by construction. The coupling score at our top-20 is small but non-zero (0.053); the coupling score at NASA's centroids is exactly zero. Both site sets disagree with the centroid as a landing site.

### So is this a methodology win or a validation-metric failure?

Both, and that's the finding. The coupling criterion *does* what the week 6 diagnostic predicted it would: it encodes the AND that linear-sum aggregation can't, and it lifts the sensitivity ceiling. But the validation metric we've been measuring against — distance to NASA centroid — penalises *any* model that picks rim cells over centroid cells, which is what NASA itself does within its regions. A more honest validation metric would compare against NASA's actual published landing-site footprints inside each region, not the disk centroid. Building that comparison is the next step (see [Roadmap](#roadmap)).

For the per-region distance table see `data/outputs/validation.json`, or run `selene validate` on a fresh checkout. The interactive map lives at [`data/outputs/webmap.html`](data/outputs/webmap.html) after `selene viz`; per-site reports under [`data/outputs/sites/`](data/outputs/sites/).

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

# Full-resolution analysis (~900 MB raw, 4 verified URLs):
selene download robbins         # ~92 MB
selene download lola            # ~115 MB
selene download illumination    # ~82 MB
selene download diviner         # ~605 MB Diviner Polar Resource Product (PRP)
# selene download lend / scarps remain TODO-flagged
selene preprocess && selene score && selene rank --top-n 20 --min-distance-km 25
```

`selene --help` lists every subcommand; `selene <cmd> --help` shows its options.

## Methodology

Every criterion produces a `[0, 1]` score grid where 1 is "best" and 0 is "unusable", aligned to the common 240 m south-polar stereographic grid (`+proj=stere +lat_0=-90 +lat_ts=-90 +R=1737400`, ±304 km, defined in [`config/region_southpole.yaml`](config/region_southpole.yaml)). Three normalisation primitives in [`scoring/normalize.py`](src/selene_base/scoring/normalize.py) — `min_max`, `optimal_range` (Gaussian), `inverse_threshold` — cover every criterion. The aggregate is a weighted linear sum that **renormalises across whichever criteria are present at score-time**, so a partial pipeline (today: slope, illumination, hazard, thermal, ice) produces a comparable score to a complete one — only the absolute meaning of "0.97" shifts.

| Criterion | Score function | Source dataset | Resolution | Resampling | Default knobs |
| --- | --- | --- | --- | --- | --- |
| **Slope** | $s = \max(0,\,1-x/\theta_{\max})$ | LOLA LDEM 80 m (PDS3) | 80 m → 240 m | bilinear | $\theta_{\max} = 15°$ |
| **Illumination** | $s = \min(x/x_t,\,1)$ | Mazarico avgvisib 65°S 240 m | 240 m | bilinear | $x_t = 0.70$ |
| **Thermal** | $s = e^{-(\bar T - T^\star)^2/(2\sigma^2)}$ on annual-mean Tavg | **Diviner PRP** `temp_avg` (PDS4) | triangle mesh → 240 m | linear griddata | $T^\star=230\,$K, $\sigma=50\,$K |
| **Ice** | $s = \mathrm{clip}(1-d/d_{\max} + \text{bonuses},\,0,\,1)$ on PRP ice-stability depth | **Diviner PRP** `ice_depth` (PDS4) | triangle mesh → 240 m | nearest griddata | $d_{\max}=2.87\,$m, surface bonus 0.5, near-PSR bonus 0.2 |
| **Hazard** | $s = \mathrm{clip}(1-d/d_{\mathrm{sat}},\,0,\,1)$ | Robbins 2018 catalog | vector → 240 m density | KDTree, 3 km radius | $d_{\mathrm{sat}}=50$ |
| **Seismic** | $s = \mathrm{clip}(\delta/\delta_{\mathrm{safe}},\,0,\,1)$ | Watters scarp catalog (TODO) | vector → 240 m distance | KDTree, 1 km densified vertices | $\delta_{\mathrm{safe}}=50\,$km |
| **Coupling** | $s = \max(0,\,1-d_{\text{PSR}}/d_c) \cdot \max(0,\,1-d_{\text{ridge}}/d_c)$ | derived: illumination + slope | 240 m | distance transform | $d_c = 5\,$km |

Slope is computed at the 240 m target resolution from the already-downsampled LOLA DEM via `numpy.gradient` with explicit metric spacing (Zevenbergen & Thorne 1987 convention; ~5 % off Horn 1981 on smooth surfaces). Computing slope on the high-res 80 m DEM and then averaging slope-degrees double-smooths and biases low; computing on the target-resolution DEM keeps everything self-consistent.

The thermal and ice criteria are both fed by the **Diviner Polar Resource Product** ([`dlre_prp_south.tab`](https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_south.tab)) — a single PDS4 character table of 2.88 M triangular-mesh facets, ~605 MB raw. [`data/pds4_table.py`](src/selene_base/data/pds4_table.py) parses it via the matching XML label; [`data/triangle_to_grid.py`](src/selene_base/data/triangle_to_grid.py) interpolates each scalar field onto the project's 240 m polar stereographic grid using `scipy.interpolate.griddata` (linear for temperatures, nearest for the discontinuous ice-depth field). Outputs are cached as three GeoTIFFs in `data/processed/` so the slow ~30 s parse step happens once.

The PSR mask used by the ice criterion is still derived from the Mazarico illumination raster (`illumination < 0.001`); the PRP is a thermal-stability calculation, not an ice-existence map, so PSR proximity adds an orthogonal signal.

Default weights from [`config/weights_default.yaml`](config/weights_default.yaml): illumination 0.20, ice 0.20, coupling 0.20, slope 0.15, thermal 0.10, hazard 0.10, seismic 0.05. The pre-week-7 weight vector (no coupling criterion) is preserved as [`config/weights_legacy.yaml`](config/weights_legacy.yaml); pass `--weights config/weights_legacy.yaml` to reproduce the 5-criterion baseline.

### The spatial-coupling criterion (week 7)

`criteria/coupling.py` is the structural fix the week 6 diagnostic identified. It scores cells by *joint proximity* to two distinct features:

1. **Distance to the nearest PSR**: derived from the Mazarico illumination raster as `illumination < 0.001`, then `scipy.ndimage.distance_transform_edt` with explicit pixel sampling.
2. **Distance to the nearest sunlit ridge**: a cell qualifies if `illumination ≥ 0.70` AND `5° ≤ slope ≤ 25°` — the geometry of polar crater rims (steeper than plains, not cliff-like, well-sunlit). Same distance transform.

The score is the **product** of two linear distance falloffs:

$$ s = \max\!\left(0,\, 1-\tfrac{d_{\text{PSR}}}{d_c}\right)\cdot \max\!\left(0,\, 1-\tfrac{d_{\text{ridge}}}{d_c}\right),\quad d_c = 5\,\text{km}. $$

The product (not sum) is the conjunction. Failing either falloff drives the score to zero — exactly the structural property a linear weighted sum cannot encode. A cell deep inside a far-side PSR (PSR distance 0, ridge distance 200 km) scores 0; a cell on a sunlit rim 30 km from the nearest PSR also scores 0; the rim cell adjacent to a PSR-floor (both distances < 5 km) scores high.

The criterion produces a **sparse mask**: only 0.12 % of finite cells exceed 0.0; only 0.07 % exceed 0.1. The polar rim band shows up clearly; the rest of the cap is essentially black:

![Coupling score, log-scaled — sparse rim band](docs/img/coupling_score.png)

The single tuning knob is `coupling_distance_km`. `selene coupling-sweep` runs the validation pipeline at 1–20 km in 8 steps and plots NASA-region alignment as a function of that knob:

![Spatial-coupling distance sweep vs alignment](docs/img/coupling_distance_sweep.png)

| coupling_distance_km | regions matched within 25 km |
| ---: | ---: |
| 1 – 10 | 0 / 9 |
| 15 – 20 | 2 / 9 |

The curve is monotone-but-flat: tightening the cap below 10 km matches no NASA regions; loosening it to 15 km picks up the Slater Plain / de Gerlache Rim 2 pair (the two NASA candidates already nearest the high-coupling band) and stays at 2/9 through 20 km. At default 5 km the criterion *is* doing what its math says (sparse rim band selection), but the band selene-base identifies is geometrically distinct from NASA's centroid points.

**Planned upgrade — TOPSIS.** A weighted linear sum lets a strong score on one criterion mask a near-disqualifying score on another. TOPSIS ranks each cell by its Euclidean distance to a synthetic "ideal" and "anti-ideal" point in criterion-score space, which penalises lop-sided profiles globally rather than at one specific spatial coupling. It's the second candidate fix from week 6 and is still on the roadmap behind a `--method topsis` flag.

## Validation

`selene validate` compares the top-N ranked sites (from `data/outputs/top_sites.geojson`) against the disk-approximation polygons of NASA's nine announced Artemis III candidate regions in [`src/selene_base/validation/nasa_regions.py`](src/selene_base/validation/nasa_regions.py). Centroids are public information from NASA's October 2024 Artemis III site-selection announcement; we approximate each region as a 15 km disk around its centroid because NASA's actual polygons are not openly published in machine-readable form. **The disks are not authoritative geometry** — they're a defensible proximity proxy for this comparison.

Two metrics for each top site:

1. **Inside any region** — does the site fall inside any of the nine 15 km disks?
2. **Within X km of any centroid** — distance from the site to the nearest NASA centroid.

And two for each NASA region:

1. **Distance to nearest top-N site** — how far away is the closest selene-base candidate?
2. **Contains a top-N site** — is at least one selene-base candidate inside this region's disk?

### Per-region results (today, 6-criterion run with coupling)

![Distance from each NASA Artemis III candidate to the nearest selene-base top site](docs/img/validation_table.png)

| NASA candidate | nearest site | distance (km) | inside region? |
| --- | --- | ---: | --- |
| Cabeus B | site_17 | 127.1 | no |
| Haworth | site_01 | 97.6 | no |
| Malapert Massif | site_01 | 115.4 | no |
| Mons Mouton | site_09 | 126.0 | no |
| Mons Mouton Plateau | site_09 | 131.0 | no |
| Nobile Rim 1 | site_02 | 64.9 | no |
| Nobile Rim 2 | site_02 | 80.8 | no |
| de Gerlache Rim 2 | site_01 | **47.8** | no |
| Slater Plain | site_01 | 55.4 | no |

The closest two pairs (de Gerlache Rim 2, Slater Plain) collapsed onto a single near-pole site — `site_01` at (-89.7°, +17.7°), score 0.707, dominant criterion *illumination*. site_01 sits in the small region of the polar cap where coupling is non-zero AND illumination is high — the kind of compromise cell the linear-sum aggregator could not previously surface as the rank-1 result. **As the Watters scarp catalog lands and the seismic criterion lights up, this table updates.**

## Robustness

Anyone reading 0/20 fairly asks: *is that just a function of the default weights?* Run `selene sensitivity --n-samples 200` to find out: it draws 200 weight vectors via Latin hypercube over the active-criteria simplex, runs `aggregate → top_n_sites → proximity_analysis` for each, and reports the distribution of "NASA regions matched within 25 km" alongside the default-weight result.

![sensitivity over 200 weight samples](docs/img/sensitivity_distribution.png)

The 6-criterion sensitivity (with coupling) **lifted the ceiling**:

- **174 / 200 samples (87 %) match 0 regions within 25 km** — the modal outcome.
- **25 / 200 samples (12.5 %) match 2 regions** — Slater Plain + de Gerlache Rim 2, the consistent close-but-not-inside pair.
- **1 / 200 sample matches 3 regions.** This bucket did not exist in the 5-criterion sweep.

The best weight regime found uses `slope = 0.00, illumination = 0.17, coupling = 0.27, thermal = 0.02, ice = 0.31, hazard = 0.23` — coupling at 0.27, *not* zero. Compare to the 5-criterion best (`hazard = 0.74, illumination = 0.22, slope = 0.03, thermal = 0.01, ice = 0.00`), which down-weighted the Diviner criteria to near-zero. **The 6-criterion best regime *uses* every new criterion**; the 5-criterion best didn't. The structural fix changed the weight-space landscape, not just the modal outcome.

The new criterion's contribution is informationally measurable: the sensitivity ceiling moved from 2/9 to 3/9, the best regime gives coupling 27 % weight, and the close-distance metric improved by a third (64.8 km → 47.8 km). It did not break through to 25 km — see [the diagnostic comparison](#diagnostic-comparison) for why.

## Diagnostic comparison

Run `selene compare` to ask: *at NASA's centroids vs at our top-20, which criteria agree and which disagree?* The table is reordered below by signed delta — agreement first, disagreement last — because that's the actual structure of the result.

![Per-criterion score: where we differ from NASA](docs/img/comparison.png)

| criterion | our top-20 | NASA 9 centroids | delta | \|t\| | reads as |
| --- | --- | --- | ---: | ---: | --- |
| **hazard** | 0.976 ± 0.026 | 0.969 ± 0.023 | **+0.007** | 0.74 | strong agreement |
| **coupling** | 0.053 ± 0.165 | 0.000 ± 0.000 | **+0.053** | 1.45 | both near zero — see below |
| **ice** | 0.993 ± 0.024 | 0.916 ± 0.096 | **+0.076** | 2.34 | strong agreement |
| **thermal** | 0.325 ± 0.198 | 0.113 ± 0.149 | **+0.212** | 3.19 | weak signal (tuning issue) |
| illumination | 0.767 ± 0.106 | 0.321 ± 0.274 | +0.446 | 4.72 | major disagreement |
| slope | 0.929 ± 0.121 | 0.285 ± 0.288 | +0.645 | 6.47 | major disagreement |

**Of the six criteria, three agree closely.** Hazard agrees almost identically (delta +0.007). Ice agrees within 0.08. Thermal scores low at both site sets (a tuning issue: the Gaussian target sits at 230 K but polar means are 100–150 K, leaving the criterion in its tail at every cell — bringing the target into the 130–150 K band would let it discriminate again).

**Coupling is the most informative single number in the project.** *NASA centroids score exactly 0.000 ± 0.000 on coupling.* Our top-20 score 0.053 — small, but non-zero. The criterion is correctly identifying the rim band where PSR meets sunlit ridge at a 5 km coupling distance — and *neither site set sits on that band*. NASA's centroids are *inside* their candidate regions (a centroid of a 15 km disk is, by construction, in the middle of the disk, which means inside a crater for Cabeus B / Haworth, on a massif for Malapert, on a plain for Slater Plain), 5–15 km off the rim that the actual NASA landing footprints would target. The 0.000 number is a discovery about the validation metric: distance to centroid is the wrong proxy for "match a NASA candidate region," because NASA's preferred landing sites within each region are off-centroid by construction.

**Slope and illumination remain the major disagreement** because the linear-sum aggregator cannot represent the *spatial coupling* between them — the methodology finding from week 6 still holds. The week 7 coupling criterion captures the conjunction directly (the product is the AND), but with only 0.20 default weight it cannot dominate the linear sum across the other five criteria. Pushing coupling weight higher would require either changing the default or — more architecturally honest — switching the aggregator itself to TOPSIS, which penalises lop-sided profiles globally rather than at one specific spatial coupling.

The ``|t|`` column is a Welch two-sample t-statistic, informational only.

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

264 tests, ~80 % combined branch coverage, all running synthetically in CI on Python 3.11 and 3.12. Real-data tests are guarded with `pytest.mark.skipif(not Path(...).exists())` so the suite stays green without ~900 MB of cached LRO data. CI runs a separate `pipeline-smoke` job on push to `main` that downloads the bundled ~12 MB sample tarball, runs `preprocess → score → rank → validate → compare`, and asserts every output file is on disk and schema-valid.

## Roadmap

- **Week 1 — data acquisition.** ✅ `selene download` for Robbins, LOLA, Mazarico illumination; Diviner / LEND / scarps URLs flagged.
- **Week 2 — common grid + slope criterion.** ✅ `reproject_to_grid`, COG cache, slope criterion end-to-end on real data.
- **Week 3 — full scoring + ranking.** ✅ All six criteria (3 on real data, 3 skip cleanly), KDTree crater density, NMS top-N extraction.
- **Week 4 — validation + visualisation.** ✅ NASA Artemis III proximity comparison, interactive folium web map, per-site HTML reports, validated v0.1.
- **Week 5 — robustness, diagnostic, sample data.** ✅ Latin-hypercube weight-sensitivity sweep, per-criterion `selene compare` diagnostic, bundled ~12 MB sample tarball, CI pipeline smoke test on the sample.
- **Week 6 — Diviner Polar Resource Product integration.** ✅ PDS4 character-table parser, triangle-mesh-to-grid rasteriser, three new score grids (`temp_avg`, `temp_max`, `ice_depth`) on the common 240 m grid, thermal+ice criteria switched to PRP defaults. Five of six criteria now run on real data. Validation rerun *surfaced the methodology finding*: linear-sum MCDA can't model NASA's spatial-coupling constraint regardless of which criteria are summed.
- **Week 7 — spatial-coupling criterion.** ✅ `criteria/coupling.py` scores cells by joint proximity to a PSR and a sunlit ridge via a multiplicative product of two distance falloffs, encoding the AND that linear-sum aggregation cannot. `selene coupling-sweep` tunes the single `coupling_distance_km` knob across 1–20 km. The fix lifted the sensitivity ceiling from 2/9 → 3/9 region matches and improved closest-distance from 64.8 km → 47.8 km, but did not break the 25 km threshold. The diagnostic surfaced a sharper finding: **NASA centroids score 0.000 on coupling**, meaning the *validation metric* (distance to disk centroid) is itself misaligned — NASA's preferred landing footprints within each region are off-centroid by 5–15 km by construction.

### Where this goes next

The week 7 result reframes the next step. Two candidate fixes:

**(a) TOPSIS aggregator behind `--method topsis`.** TOPSIS (Technique for Order of Preference by Similarity to Ideal Solution) ranks each cell by Euclidean distance to a synthetic "ideal" and "anti-ideal" point in criterion-score space. It penalises lop-sided profiles globally rather than at one specific spatial coupling — complementary to the coupling criterion, not a substitute for it. Easy to drop in behind a CLI flag.

**(b) A better validation metric.** The week 7 diagnostic showed the 0.000 coupling score at NASA centroids is itself a discovery: the *centroid* is the wrong proxy for "is this NASA candidate well-served by selene-base." NASA's published Artemis III selection materials describe each region as a ~15 km disk *containing* preferred landing footprints — the actual operational targets. Building a polygon-based proximity metric (does any selene-base top-20 site fall *anywhere* within the disk, not just within 25 km of the centroid?) would change `0/20` to a meaningful number on the same outputs we already produce. This is the more honest fix, and it doesn't touch the criterion stack — the change is in `validation/comparison.py`.

The intent is both: TOPSIS as a `--method topsis` flag, plus a polygon-based "inside any region" metric refinement. Validation reruns against the same nine NASA candidates after each.

### Smaller follow-ups

- Resolve the last TODO URL — the Watters lobate-scarp catalog — and light up the seismic criterion (six of six live).
- **Thermal target re-tune**: the PRP `temp_avg` peaks at 211 K against a 230 K Gaussian peak, leaving the criterion in its tail at every cell. Bringing the target into the 130–150 K band would let thermal contribute discriminative signal again.
- Earth line-of-sight criterion derived from LOLA elevation horizon checks.
- ML-based criterion inputs (planned as a separate project, `selene-vision`).

## References

- Robbins, S. J. (2019). *A new global database of lunar impact craters >1–2 km: 1. Crater locations and sizes, comparisons with published databases, and global analysis.* Journal of Geophysical Research: Planets, 124, 871–892. [doi:10.1029/2018JE005592](https://doi.org/10.1029/2018JE005592)
- Mazarico, E., Neumann, G. A., Smith, D. E., Zuber, M. T., & Torrence, M. H. (2011). *Illumination conditions of the lunar polar regions using LOLA topography.* Icarus, 211(2), 1066–1081. [doi:10.1016/j.icarus.2010.10.030](https://doi.org/10.1016/j.icarus.2010.10.030)
- Smith, D. E., et al. (2010). *The Lunar Orbiter Laser Altimeter investigation on the Lunar Reconnaissance Orbiter mission.* Space Science Reviews, 150(1–4), 209–241. [doi:10.1007/s11214-009-9512-y](https://doi.org/10.1007/s11214-009-9512-y)
- Paige, D. A., et al. (2010). *The Lunar Reconnaissance Orbiter Diviner Lunar Radiometer Experiment.* Space Science Reviews, 150(1–4), 125–160. [doi:10.1007/s11214-009-9529-2](https://doi.org/10.1007/s11214-009-9529-2)
- Williams, J.-P., et al. (2017). *The global surface temperatures of the Moon as measured by the Diviner Lunar Radiometer Experiment.* Icarus, 283, 300–325. (PRP modeled-ice-stability methodology.)
- Diviner Polar Resource Product (PRP), south pole, version 1.0. PDS Geosciences Node Diviner derived bundle: [`dlre_prp_south.tab`](https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_south.tab).
- Mitrofanov, I. G., et al. (2010). *Hydrogen mapping of the lunar south pole using the LRO Neutron Detector Experiment LEND.* Science, 330(6003), 483–486.
- Watters, T. R., Robinson, M. S., Banks, M. E., Tran, T., & Denevi, B. W. (2015). *Global thrust faulting on the Moon and the influence of tidal stresses.* Geology, 43(10), 851–854. [doi:10.1130/G37120.1](https://doi.org/10.1130/G37120.1)
- Civilini, F., Weber, R. C., Jiang, Z., Phillips, D., & Pan, W. (2023). *Constraints on the seismic hazard of young thrust faults on the Moon from re-located shallow moonquakes.* (Used as motivation for the seismic exclusion criterion.)
- NASA (October 2024). *Artemis III candidate landing regions.* [https://www.nasa.gov/feature/artemis-iii](https://www.nasa.gov/feature/artemis-iii)

## Notes for the reader

The interactive web map (`data/outputs/webmap.html` after `selene viz`) is built with folium / Leaflet and pulls Leaflet's JS and CSS from a CDN — so it needs an internet connection on first open. Everything else (the score raster, polygons, popups, per-site reports) is inlined and works offline. Per-site HTML reports under `data/outputs/sites/` are fully self-contained.

## License

MIT — see [LICENSE](LICENSE).
