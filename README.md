# selene-base

> Multi-criteria habitat suitability for the lunar south pole, validated against authoritative USGS-published Artemis IV (formerly Artemis III) candidate region polygons (DOI 10.5066/P1MEQ6UK). v1.5 produces 69 HLS-compliant candidate landing sites across 8 of NASA's 9 Artemis IV candidate regions at 20 m resolution (Wueller-class); v1.8 activates the eighth (seismic) criterion on the same 69-site catalog, of which **56 (81.2 %) agree within 5 km of a peer-reviewed Wueller et al. 2026 site (JGR Planets, 130 sites); median match distance 1.69 km** — quantitative agreement with peer-reviewed published methodology, computed against the authors' Zenodo data deposit (CC-BY 4.0). v1.5 ships GPU-accelerated horizon-profile derivation (CuPy on Blackwell) enabling per-region high-resolution analysis; v1.7 adds TOPSIS as an opt-in alternative aggregator; v1.8 lights up the eighth (seismic) criterion on a bundled south-polar scarp catalog. The release history below documents the engineering arc across eleven versioned releases.

[![CI](https://github.com/Alex0420W/selene-base/actions/workflows/ci.yml/badge.svg)](https://github.com/Alex0420W/selene-base/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)](pyproject.toml)
[![Release](https://img.shields.io/github/v/release/Alex0420W/selene-base)](https://github.com/Alex0420W/selene-base/releases)

NASA's Artemis IV mission will land humans near the lunar south pole in early 2028 (the program was restructured in February 2026, moving the first crewed landing from Artemis III to Artemis IV; the candidate landing regions identified in NASA's October 2024 down-selection are unchanged). Selecting a base site there is a multi-criteria optimisation problem: the south pole is a maze of crater rims that catch grazing sunlight, deep permanently-shadowed cold-traps that may host water ice, and active thrust faults that re-localised Apollo-era shallow moonquakes have placed within tens of kilometres of candidate sites. **`selene-base`** fuses the modern LRO-era remote-sensing record (LOLA topography, Diviner Polar Resource Product thermal and ice-stability data, Mazarico illumination maps, the Robbins crater catalog) with derived geometric criteria (Earth line-of-sight from horizon ray-marching, spatial coupling between PSRs and sunlit ridges) to score every 240 m pixel of the polar cap and rank top candidate sites. The pipeline is end-to-end reproducible — `selene download && selene preprocess && selene score && selene rank-per-region && selene validate-per-region && selene viz` produces the NASA-aligned per-region catalog, per-site HTML reports, and an interactive web map, on a developer laptop, in minutes, from public data.

## Headline finding

> **56 of 69 selene-base candidate landing sites (81.2 %) match within 5 km of a peer-reviewed Wueller et al. 2026 site, with median match distance 1.69 km (v1.8, all eight criteria active; tightens from 1.76 km at v1.5).** The 69 sites are the per-region HLS-compliant landing catalog at 20 m resolution: candidates across 8 of NASA's 9 Artemis IV candidate regions, all guaranteed inside their published USGS polygon and satisfying NASA's published HLS hard-constraint filters (slope ≤ 8°, 100 m buffer to steeper terrain, illumination ≥ 33 %, and DTE visibility ≥ 50 %) by construction. The sole region with zero compliant cells is **Malapert Massif** — a real terrain-driven finding that persists at 20 m, confirming the v1.4.2 result was not a 240 m sampling artefact (no cell inside the Malapert polygon simultaneously satisfies all four HLS thresholds at either resolution). The best-scoring region is **Mons Mouton Plateau** at score 0.732 with the largest HLS-eligible area (171.95 km², 3.86 % of polygon); the most constrained region with sites is **de Gerlache Rim 2**, where only 0.01 % of polygon cells pass the HLS filters at 20 m and only 2 sites fit at the 2 km NMS separation — also unchanged from 240 m, confirming the dGR2 disagreement with Wueller's catalog is genuine terrain divergence rather than a resolution artefact.

[![selene-base v1.5 (20 m, per-region tiled): 69 candidate sites, 81.2 % within 5 km of Wueller 2026 (median 1.76 km), with USGS polygon outlines and matched-pair connectors](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/selene_vs_wueller_20m.png)](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/selene_vs_wueller_20m.png)

This headline is a **reframing**, not a refinement. The previous versions (v1.0.0–v1.2.0) measured a different question — global ranking of habitat suitability followed by validation against NASA polygons — and reported 0/20 because globally-selected top sites pick the polar rim band where the coupling criterion is non-zero rather than the *interior* of NASA-published regions. The Wueller et al. 2026 (JGR Planets, [doi:10.1029/2025JE009434](https://doi.org/10.1029/2025JE009434)) parallel — which found 130 candidate sites with similar within-region HLS-filtered methodology — confirms that **per-region HLS-filtered ranking is the right framing for a NASA-aligned site catalog**.

The earlier releases (v1.0.0–v1.2.0) are preserved as the project's diagnostic history — they document *why* the global-ranking framing produced 0 inside-polygon results before the v1.3 reframing. The summary below tracks the releases in order; the v1.3 reframing supersedes the global-ranking interpretation but does not invalidate the diagnostic findings.

### Diagnostic findings (v1.0.0 — v1.2.0): why global ranking gave 0 / 20

Across the v1.0.0 — v1.2.0 releases — through five revisions of the model and validation reference against NASA's nine announced Artemis III candidate regions — the **global**-ranking inside-polygon count was 0 / 20:

1. **Adding well-validated criteria can degrade alignment.** Integrating Diviner PRP thermal and ice criteria — both showing strong individual agreement with NASA centroids — *worsened* median distance from 65 km to 102 km. Weighted-sum MCDA cannot model NASA's coupled spatial constraint (near a PSR *and* near a sunlit ridge); it lets a high score on one axis compensate for a near-zero on another, producing top sites that maximize independent criteria but rarely satisfy the conjunction.

2. **A spatial-coupling criterion is the structural fix.** Modeling the conjunction directly (the multiplicative product of two distance falloffs, encoding the AND that linear-sum cannot) lifted the 200-sample sensitivity ceiling from 2/9 to 3/9 NASA regions matched, improved closest-distance from 64.8 km to 47.8 km, and the best-found weight regime now uses coupling at 0.27 weight rather than ignoring it.

3. **A polygon-based validation metric and a thermal-target correction (v1.0.0) tightened the methodology, not the headline.** Replacing the centroid-distance proxy with an "inside any 15 km disk" check and correcting the thermal-criterion target from the data-out-of-support 230 K to the data-median 140 K both produced *informational* shifts: the inside-any-disk count is also 0/20, but the closest-edge distance is now 32.8 km (about disk-radius below the centroid distance, exactly as the geometry predicts); the corrected thermal moved out of the Gaussian's tail (NASA centroids 0.113 → 0.526; our top-20 0.325 → 0.965), and 11/200 sensitivity samples now reach 3/9 region matches (was 1/200). Two methodological loose ends closed; the headline didn't move because the rim band our model identifies is geometrically distinct from NASA's centroid disks at any reasonable proximity threshold.

4. **An Earth line-of-sight criterion (v1.1.0) narrows the gap without closing it.** Adding the most physically prominent missing criterion — a per-pixel Earth-visibility fraction derived from a 36-azimuth horizon ray-march on the LOLA elevation grid plus libration-cycle sampling of the sub-Earth point — shifted the closest-distance from 47.8 km to **27.3 km** (centroid) and from 32.8 km to **12.3 km** (disk edge), put **1/9 NASA regions within one disk radius of a top site** for the first time (Cabeus B), lifted the sensitivity ceiling from 3/9 to **4/9** regions matched, and produced **the first non-zero polygon-inside count in any sensitivity run (21/200 samples now show ≥1 site inside a NASA disk)** — though the default-weights polygon-inside count is still **0/20**. The geometric gap is no longer absolute; it is, however, still *robust* under the operations-driven defaults.

5. **Replacing the 15 km disk approximations with USGS's authoritative published polygons (v1.2.0) confirms the geometric separation is not a disk-approximation artefact.** USGS Data Release 10.5066/P1MEQ6UK ships simplified envelopes (4-vertex quadrilaterals) for all nine Artemis III regions, derived from NASA's LROC QuickMap definitions. They differ substantively from the disk approximations: most regions are ~400 km² quadrilaterals (vs the disk's 707 km²); Mons Mouton Plateau alone is **4,452 km² — over 6× larger than the disk**; one region the legacy code called "Cabeus B" is published as "Peak Near Cabeus B" centred on the rim, not the crater floor; and one disk centroid (the legacy "Slater Plain" at lon -54.3°) sits ~180° away from where USGS publishes Slater Plain (lon +125°). Against these authoritative polygons, the default-weights result is **0/20 top sites inside any USGS polygon, 0/9 USGS regions containing a top site, median distance to the nearest USGS polygon 135.1 km, closest 41.5 km (de Gerlache Rim 2)**. The 200-sample sensitivity sweep produces ≥1 site inside a USGS polygon in **6 / 200 samples** (max 2/20). **The geometric separation is real even against the right validation reference** — the disk approximations were systematically misrepresenting the regions, but the model's rim-band optimum is still geometrically distinct from NASA's authoritative regions.

6. **The right *framing* — per-region ranking with HLS hard filters (v1.3.0) — produces a per-region landing-site catalog across 8/9 NASA regions** (23 sites at the v1.3 default of 3 per region, 70 sites at the v1.4.2 default of 10 per region). The first five stages of the arc all rank globally and ask "did our top-20 fall inside NASA's regions?" The right question, mirroring NASA's own selection process, is "within each NASA region, which cells satisfy the HLS landing requirements and rank highest by suitability?" Reframing the search this way produces a complete per-region landing-site catalog. **The 0/20 result through v1.2.0 reflected the global-ranking framing, not a flaw in the model**: globally-best cells cluster on the polar rim band where the coupling criterion is non-zero; per-region-best HLS-compliant cells cluster inside NASA polygons by construction. Both findings are valid — the diagnostic arc is preserved below for context.

The combined picture: the criteria are tuned, the validation primitive is authoritative, and **the framing is now NASA-aligned**. The previously-open scientific question — "do the v1.3 sites identify the same cells NASA's process identifies?" — is now answered quantitatively in v1.4.2 against Wueller et al. 2026's 130 published sites: **56/70 selene sites (80 %) match within 5 km of an in-scope Wueller site; median match distance 1.88 km** (see "Quantitative comparison against Wueller et al. 2026" below).

![Diagnostic from v1.0–v1.2: model's top-20 (red) cluster on the rim band where the coupling criterion is non-zero; NASA centroids (cyan) sit in the disk interiors](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/coupling_overlay.png)

### Validation history across releases (global ranking, v1.0.0 — v1.2.0)

The table below tracks the *global*-ranking inside-polygon count across the v1.0.0–v1.2.0 releases. v1.3.0 reframes the analysis to per-region ranking; the global-ranking inside-polygon count remains 0/20 by construction of the global framing, but is no longer the project's headline.

| Validation metric | v0.5 baseline | v0.6 (+PRP) | v0.7 (+coupling) | v1.0 (+poly) | v1.1 (+LOS) | **v1.2 (USGS)** |
| --- | --- | --- | --- | --- | --- | --- |
| top sites within 25 km of any centroid | 0 / 20 | 0 / 20 | 0 / 20 | 0 / 20 | 0 / 20 | **0 / 20** |
| top sites inside any 15 km disk | n/a | n/a | n/a | 0 / 20 | 0 / 20 | 0 / 20 |
| **top sites inside any USGS polygon** | n/a | n/a | n/a | n/a | n/a | **0 / 20** |
| **USGS regions containing a top site** | n/a | n/a | n/a | n/a | n/a | **0 / 9** |
| **closest USGS polygon (km)** | n/a | n/a | n/a | n/a | n/a | **41.5** |
| **median distance to nearest USGS polygon (km)** | n/a | n/a | n/a | n/a | n/a | **135.1** |
| regions with a top site within 1 disk radius of edge | n/a | n/a | n/a | 0 / 9 | 1 / 9 | 1 / 9 |
| closest NASA region (centroid) | 25.8 km | 64.8 km | 47.8 km | 47.8 km | **27.3 km** | 27.3 km |
| closest NASA region (disk edge) | n/a | n/a | n/a | 32.8 km | **12.3 km** | 12.3 km |
| 200-sample sensitivity, *best* regions matched (centroid≤25 km) | 2 / 9 | 2 / 9 | 3 / 9 (1 sample) | 3 / 9 (11 samples) | 4 / 9 (5 samples) | **4 / 9 (5 samples)** |
| 200-sample sensitivity, modal outcome | 0/9 (71.5 %) | 0/9 (92.5 %) | 0/9 (87 %) | 0/9 (86.5 %) | 0/9 (78.5 %) | 0/9 (78.5 %) |
| 200-sample sensitivity, samples with ≥1 site inside any disk | n/a | n/a | n/a | 0 / 200 | 21 / 200 | 21 / 200 |
| **200-sample sensitivity, samples with ≥1 site inside any USGS polygon** | n/a | n/a | n/a | n/a | n/a | **6 / 200** |
| coupling score at our top-20 vs NASA centroids | n/a | n/a | 0.053 vs 0.000 | 0.025 vs 0.000 | 0.025 vs 0.000 | 0.025 vs 0.000 |
| thermal score at our top-20 vs NASA centroids | n/a | 0.239 vs 0.113 (in tail) | 0.325 vs 0.113 (in tail) | 0.965 vs 0.526 | 0.962 vs 0.526 | 0.962 vs 0.526 |
| LOS-to-Earth score at our top-20 vs NASA centroids | n/a | n/a | n/a | n/a | 1.000 vs 0.525 (bimodal) | 1.000 vs 0.525 |

The default-weights inside-region number is *0 across every release in the table, against every validation primitive tried*. What moved is the structure behind that number — the sensitivity ceiling lifted twice (2 → 3 → 4 regions), the closest-edge metric collapsed (47.8 → 12.3 km against disks; the v1.2 USGS polygons reset the closest-distance to 41.5 km because USGS's "Slater Plain" sits 180° away in longitude from where the legacy disk centroid placed it), **the polygon-inside count became sample-non-zero against disks** (21/200 weight regimes hit 1/20 or more), and one NASA region (Cabeus B) came within disk-radius of a top site. v1.2.0 against the authoritative USGS polygons confirms the geometric separation is *not* a disk-approximation artefact: even with the right validation reference, only 6/200 weight regimes produce ≥1 site inside any USGS polygon, and the default-weights answer is 0/20. **Under the global-ranking framing, the model picks the polar rim band where the coupling criterion is non-zero; NASA's published landing-region polygons sit 41–135 km away from that band; no setting of physics-driven defaults collapses those two geometries onto each other.** v1.3.0 reframes the analysis to per-region ranking and resolves this — see "Per-region HLS-compliant catalog" below.

For the per-region distance table see `data/outputs/validation.json` (now contains both centroid-distance and polygon-inside metrics), or run `selene validate` on a fresh checkout. The interactive map lives at [`data/outputs/webmap.html`](data/outputs/webmap.html) after `selene viz`; per-site reports under [`data/outputs/sites/`](data/outputs/sites/).

## Pipeline

```
data/raw/<dataset>/        --load-->  xr.DataArray (native CRS)
                              |
                              v reproject_to_grid(target_crs, bounds, 240 m)
                              |
data/processed/<name>_southpole_240m.tif        (cached COG)
                              |
                              v criterion.compute(...)            [eight criteria]
                              |
data/processed/scored/<name>_score_southpole_240m.tif
                              |
                              v scoring.aggregate.weighted_sum()  [renormalises]
                              |
data/outputs/score_southpole.tif                (final aggregate COG)
                              |
                              +--> scoring.ranking.top_n_sites_per_region()  [v1.3 primary]
                              |    + HLS hard filters per USGS polygon
                              |    -> data/outputs/per_region/sites.{geojson,csv}
                              |
                              +--> scoring.ranking.top_n_sites()             [legacy global]
                              |    -> data/outputs/top_sites.{geojson,csv}
                              |
                              v validation.comparison + viz
                              |
data/outputs/validation.json + webmap.html + sites/
```

### Quickstart

```bash
git clone https://github.com/Alex0420W/selene-base.git
cd selene-base
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e .

# Five-line clone-to-webmap path on the bundled ~12 MB sample dataset
# (a single LOLA tile, the south-pole subset of the Mazarico illumination
# raster, and a 1000-row crater catalog excerpt — enough to exercise the
# full preprocess -> score -> rank pipeline end-to-end):
selene download --sample        # downloads + extracts data/raw/<sample>
selene preprocess               # warps + crater-density rasterisation -> data/processed/
selene score                    # seven criteria; missing ones renormalise out cleanly
selene rank-per-region          # NASA-aligned per-USGS-polygon catalog (default n_per_region=10)
selene viz                      # webmap.html + per-site HTML reports
                                # (note: webmap loads Leaflet from a CDN on
                                #  first open; everything else is offline-safe)

# Diagnostic & robustness (legacy global-ranking framing):
selene rank --top-n 20          # global NMS top-N
selene validate                 # alignment metrics vs NASA's nine candidates
selene validate-per-region      # per-region summary: n sites, eligible-area %, best score
selene compare                  # per-criterion delta our top-20 vs NASA centroids
selene sensitivity --n-samples 200       # 200-sample weight-vector simplex sweep
selene coupling-sweep                    # tune coupling_distance_km against alignment
selene compare-wueller          # quantitative comparison vs Wueller 2026 (130 real sites, CC-BY 4.0)

# Full-resolution analysis (~900 MB raw, 4 verified URLs):
selene download robbins         # ~92 MB
selene download lola            # ~115 MB
selene download illumination    # ~82 MB
selene download diviner         # ~605 MB Diviner Polar Resource Product (PRP)
# selene download lend remains TODO-flagged (no verified URL); the seismic-
# criterion scarp catalog (Mishra & Kumar 2022) ships bundled in-repo as of
# v1.8 — no separate download step needed.
selene preprocess && selene score && selene rank-per-region
```

`selene --help` lists every subcommand; `selene <cmd> --help` shows its options.

## Methodology

Every criterion produces a `[0, 1]` score grid where 1 is "best" and 0 is "unusable", aligned to the common 240 m south-polar stereographic grid (`+proj=stere +lat_0=-90 +lat_ts=-90 +R=1737400`, ±304 km, defined in [`config/region_southpole.yaml`](config/region_southpole.yaml)). Three normalisation primitives in [`scoring/normalize.py`](src/selene_base/scoring/normalize.py) — `min_max`, `optimal_range` (Gaussian), `inverse_threshold` — cover every criterion. The aggregate is a weighted linear sum that **renormalises across whichever criteria are present at score-time**, so a partial pipeline produces a comparable score to a complete one — only the absolute meaning of "0.97" shifts. As of v1.8 the methodology is fully realised: all eight criteria (slope, illumination, hazard, thermal, ice, coupling, los_to_earth, seismic) contribute to the aggregate score by default.

| Criterion | Score function | Source dataset | Resolution | Resampling | Default knobs |
| --- | --- | --- | --- | --- | --- |
| **Slope** | $s = \max(0,\,1-x/\theta_{\max})$ | LOLA LDEM 80 m (PDS3) | 80 m -> 240 m | bilinear | $\theta_{\max} = 15°$ |
| **Illumination** | $s = \min(x/x_t,\,1)$ | Mazarico avgvisib 65°S 240 m | 240 m | bilinear | $x_t = 0.70$ |
| **Thermal** | $s = e^{-(\bar T - T^\star)^2/(2\sigma^2)}$ on annual-mean Tavg | **Diviner PRP** `temp_avg` (PDS4) | triangle mesh -> 240 m | linear griddata | $T^\star=140\,$K, $\sigma=30\,$K (v1.0.0 correction; was 230 / 50, outside data support) |
| **Ice** | $s = \mathrm{clip}(1-d/d_{\max} + \text{bonuses},\,0,\,1)$ on PRP ice-stability depth | **Diviner PRP** `ice_depth` (PDS4) | triangle mesh -> 240 m | nearest griddata | $d_{\max}=2.87\,$m, surface bonus 0.5, near-PSR bonus 0.2 |
| **Hazard** | $s = \mathrm{clip}(1-d/d_{\mathrm{sat}},\,0,\,1)$ | Robbins 2019 catalog | vector -> 240 m density | KDTree, 3 km radius | $d_{\mathrm{sat}}=50$ |
| **Coupling** | $s = \max(0,\,1-d_{\text{PSR}}/d_c) \cdot \max(0,\,1-d_{\text{ridge}}/d_c)$ | derived: illumination + slope | 240 m | distance transform | $d_c = 5\,$km |
| **LOS-to-Earth** | linear ramp on per-pixel Earth visibility fraction over libration | derived: LOLA elevation (horizon profile) | 240 m | bilinear ray sampling | $\text{vis}_{\min}=0.20,\,\text{vis}_{\text{target}}=0.50$ |
| **Seismic** (v1.8) | $s = 1/(1+\exp(-(d-d_0)/k))$ on distance to nearest scarp | Mishra & Kumar 2022 (bundled, primary) + Watters 2015 polar (bundled, attribution anchor) | vector -> 240 m distance | KDTree, 1 km densified vertices | $d_0 = 25\,$km, $k = 8\,$km |

Slope is computed at the 240 m target resolution from the already-downsampled LOLA DEM via `numpy.gradient` with explicit metric spacing (Zevenbergen & Thorne 1987 convention; ~5 % off Horn 1981 on smooth surfaces). Computing slope on the high-res 80 m DEM and then averaging slope-degrees double-smooths and biases low; computing on the target-resolution DEM keeps everything self-consistent.

The thermal and ice criteria are both fed by the **Diviner Polar Resource Product** ([`dlre_prp_south.tab`](https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_south.tab)) — a single PDS4 character table of 2.88 M triangular-mesh facets, ~605 MB raw. [`data/pds4_table.py`](src/selene_base/data/pds4_table.py) parses it via the matching XML label; [`data/triangle_to_grid.py`](src/selene_base/data/triangle_to_grid.py) interpolates each scalar field onto the project's 240 m polar stereographic grid using `scipy.interpolate.griddata` (linear for temperatures, nearest for the discontinuous ice-depth field). Outputs are cached as three GeoTIFFs in `data/processed/` so the slow ~30 s parse step happens once.

The PSR mask used by the ice criterion is still derived from the Mazarico illumination raster (`illumination < 0.001`); the PRP is a thermal-stability calculation, not an ice-existence map, so PSR proximity adds an orthogonal signal.

Default weights from [`config/weights_default.yaml`](config/weights_default.yaml) (v1.8 — eighth criterion live): illumination 0.17, ice 0.17, coupling 0.17, los_to_earth 0.14, slope 0.12, thermal 0.07, hazard 0.06, seismic 0.10. **The 0.10 seismic weight reflects that scarp distance is one of multiple safety factors and is most relevant on long habitat-occupation timescales**; the other seven weights are scaled by 0.928 from the v1.5 vector so the simplex sums to 1.00 without changing relative ratios. **The 0.15 LOS weight (now 0.14 after seismic activation) was chosen *before* the v1.1.0 validation rerun on physics-and-operations grounds**: LOS is a real physical criterion NASA prioritises in mission planning, but more binary-ish than continuous, so weighted slightly less than illumination/ice/coupling. Older weight vectors are preserved: [`config/weights_legacy.yaml`](config/weights_legacy.yaml) (5-criterion baseline, no coupling/LOS), [`config/weights_legacy_v6.yaml`](config/weights_legacy_v6.yaml) (v1.0.0 6-criterion, no LOS).

### Mapping to NASA's Figures of Merit framework

NASA's Artemis site-selection process ([Lawrence 2025](https://ntrs.nasa.gov/citations/20250008952), NTRS 20250008952) applies "weighted Figures of Merit" across multiple multi-directorate assessment categories. The specific FOM values are not publicly disclosed; what is documented is the conceptual framework. selene-base's seven physics-driven criteria can be mapped to NASA's FOM categories as follows:

| selene-base criterion | NASA FOM category | data source agreement | NASA values |
|---|---|---|---|
| Slope (≤ 8 °, 100 m buffer) | Site Availability + Deorbit/Descent/Landing | LOLA-derived (consistent) | published threshold (HLS) |
| Illumination | Mission Availability | Mazarico raster (consistent) | published threshold (HLS) |
| Thermal | Science Objectives | Diviner PRP (consistent) | not public |
| Ice (PSR access) | Science Objectives | Diviner PRP + Mazarico (consistent) | not public |
| Hazard (crater density) | Site Availability | Robbins 2019 catalog (consistent) | not public |
| Coupling (PSR ↔ ridge) | Mission Availability + Science Objectives | derived (Mazarico + slope) | not public |
| LOS-to-Earth | Mission Availability | derived (LOLA + libration) | published threshold (HLS) |
| Seismic (Watters scarps) | Site Availability | Mishra & Kumar 2022 + Watters 2015 (consistent; bundled in-repo since v1.8) | not public |

selene-base implements *physics-driven analogs* of NASA's FOM categories using public data and operations-driven thresholds, not direct implementations of NASA's actual weighting (which is not public). The HLS hard-constraint thresholds (slope, slope-buffer, illumination, DTE visibility) are published and are matched verbatim. The soft-criterion weights are physics-and-operations-driven defaults set before any validation rerun, not tuned to match a NASA reference. This positions selene-base honestly: a complete public-data implementation of the conceptual framework documented by Lawrence 2025, calibrated to physics-driven defaults rather than fitted to undisclosed NASA-internal weights.

### Seismic criterion (v1.8)

[Civilini et al. (2023)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019GL083580) re-located several Apollo-era shallow moonquakes and showed several of them cluster within tens of kilometres of mapped lunar lobate scarps — the same fault-scarp landforms catalogued globally by Watters et al. (2015) and extended for the south polar region by Mishra & Kumar (2022). v1.8 activates the eighth criterion as **distance-to-nearest-scarp**, with a logistic mapping centred on Civilini's documented moonquake-to-scarp clustering distance:

- cells within **5 km** of a mapped scarp score ~0.08 (high seismic risk),
- cells at **25 km** score 0.50 (transition; matches typical Apollo-era epicentral-uncertainty scale),
- cells beyond **50 km** score ~0.96 (effectively safe under the Civilini regime).

The default weight is 0.10 (slightly below slope, reflecting that scarp distance is one of multiple safety factors most relevant on long habitat-occupation timescales). The logistic shape is preferred over a hard linear ramp because Civilini's moonquake-to-scarp distance estimates are themselves probabilistic with epicentral uncertainty of order 10 km.

The primary data source is **[Mishra & Kumar (2022)](https://doi.org/10.1029/2022GL098505)** (*Geophysical Research Letters*, [doi:10.5281/zenodo.6624114](https://doi.org/10.5281/zenodo.6624114), CC-BY 4.0): 704 south-polar lobate-scarp line segments combining Watters' published mapping with 75 new scarps from independent LROC NAC analysis. The 704 segments cover all 9 USGS Artemis IV candidate regions plus the 100 km horizon-march buffer. The original LROC SOC `POLAR_SCARP_LOCATIONS` shapefile (citing Watters et al. 2015 verbatim, NASA-hosted) is bundled alongside as an attribution anchor and sanity reference at [`src/selene_base/criteria/data/scarps_watters_2015_polar/`](src/selene_base/criteria/data/scarps_watters_2015_polar/).

**Result against Wueller 2026 with all eight criteria active (20 m per-region tiled run on the v1.5-cached preprocessing):** 56 / 69 selene sites (81.2 %) match within 5 km of an in-scope Wueller site (identical to v1.5's 7-criterion result); 46 / 73 (63.0 %) of in-scope Wueller sites match a selene site within 5 km (**up from 44 / 73 in v1.5**, +2.7 pp); median matched-pair distance **1.69 km vs 1.76 km in v1.5** (-0.07 km, tightens). Same 69-site catalog by site count and same top-3 sites in every region; 13 sites in HW/MMP/N1/N2 reshuffle ranks 4–10 within their polygons because the seismic component lifts the aggregate of slightly-further-from-scarp cells just enough to flip a few interior orderings.

### Aggregator: weighted_sum vs TOPSIS (v1.7)

The default aggregator (v1.4+) is a weighted linear sum:
$\text{score}(c) = \sum_i w_i \cdot s_i(c) / \sum_i w_i$,
where $s_i(c) \in [0, 1]$ is the score of cell $c$ on criterion $i$
and $w_i$ is its weight (renormalised across criteria with score
grids — see [`scoring/aggregate.weighted_sum`](src/selene_base/scoring/aggregate.py)).
The linear sum is simple, matches Wueller et al. 2026's framework,
and is *indifferent to balance*: a cell that saturates one criterion
at 1.0 and zeroes another scores the same as a cell that scores 0.5
on both, given equal weights.

v1.7 adds [`scoring/aggregate.topsis`](src/selene_base/scoring/aggregate.py)
as an alternative aggregator. TOPSIS (Hwang & Yoon 1981) vector-
normalises each criterion to L2 unit length, applies the weights,
identifies the per-criterion ideal (max) and anti-ideal (min) of the
weighted-normalised grid, and scores each cell by closeness:
$\text{score}(c) = d_\text{anti}(c) / (d_\text{ideal}(c) + d_\text{anti}(c))$.
Cells at the per-criterion ideal score 1.0; cells at the anti-ideal
score 0.0; balanced "good across the board" cells outscore lopsided
ones because Euclidean distance to a multidimensional anti-ideal
penalises near-zero scores more than the linear sum does.

**Result on the v1.5 20 m catalog (same HLS filters, same per-region
NMS, only the aggregator changes):**

| metric | v1.5 (weighted_sum) | v1.7 (TOPSIS) | Δ |
| --- | ---: | ---: | ---: |
| selene sites | 69 | 69 | 0 (HLS filters unchanged) |
| selene matched within 5 km of Wueller in-scope | 56 / 69 (81.2 %) | 55 / 69 (79.7 %) | -1 site, -1.5 pp |
| Wueller matched within 5 km of selene | 44 / 73 (60.3 %) | 44 / 73 (60.3 %) | 0 |
| **median matched-pair distance** | 1.76 km | 1.84 km | +0.08 km |
| max matched-pair distance | 4.81 km | 4.59 km | **-0.22 km** |
| top-5-by-score sites | 23, 24, 25, 26, 27 (Mons Mouton Plateau) | 23, 24, 25, 26, 27 (Mons Mouton Plateau) | identical |

Per-region median pair distance:

| region | weighted_sum (v1.5) | TOPSIS (v1.7) | Δ |
| --- | ---: | ---: | ---: |
| Haworth | 1.77 km | **1.26 km** | **-0.51 km (tightens)** |
| Mons Mouton | 1.90 km | 1.95 km | +0.05 km |
| Mons Mouton Plateau | 2.03 km | **1.63 km** | **-0.40 km (tightens)** |
| Nobile Rim 1 | 1.32 km | **0.88 km** | **-0.43 km (tightens)** |
| Nobile Rim 2 | 1.13 km | 1.30 km | +0.17 km |
| Peak Near Cabeus B | 2.04 km | 2.04 km | 0 |
| Slater Plain | 2.03 km | 2.01 km | -0.02 km |
| de Gerlache Rim 2 | — (0/2) | — (0/2) | — |

**Headline reading:** TOPSIS produces an **identical site catalog** (same 69 cells, same top-5) but a **different per-region distance distribution**. Three of seven regions with matched pairs tighten meaningfully (HW, MMP, N1, all by -0.40 to -0.51 km); two widen slightly (MM, N2, both ≤ +0.17 km). The headline 81.2 % → 79.7 % match rate change is one site flipping just past the 5 km threshold; max-pair distance *improves* under TOPSIS. **Both aggregators agree on which cells are top-tier**; they disagree at the margin where 1.5–4 km matter.

Run TOPSIS via `--method topsis` on `selene score` and `selene rank-per-region`:

```
selene score --method topsis --outputs-dir data/outputs/topsis
selene rank-per-region --method topsis --tiled-per-region --resolution 20
selene compare-wueller --sites data/outputs/topsis/per_region_tiled_topsis/sites.geojson \
    --outputs-dir data/outputs/topsis/v17
```

Outputs land under `data/outputs/topsis/` so the TOPSIS catalog sits side-by-side with the weighted-sum default at `data/outputs/per_region_tiled/` without clobbering it. Default stays `weighted_sum` for backward compatibility — v1.5's headline number is unchanged.

### Per-region vs global ranking (v1.3.0)

Through v1.2.0, the pipeline ran **global ranking**: the aggregate score grid is NMS-extracted across the entire ±304 km polar grid, producing a top-20 list that *may or may not* fall inside NASA polygons. Validation then asked "did our top-20 hit the polygons?" — and reported 0/20 across every revision through v1.2.0.

**v1.3.0 reframes the analysis to match NASA's actual selection process.** `selene rank-per-region` searches *within* each USGS polygon, applying NASA's published Human Landing System (HLS) hard-constraint filters as a precondition before ranking by suitability score. The result is a per-region landing-site catalog: every site is inside its named polygon by construction; every site satisfies the published HLS thresholds; the soft-criterion aggregate score is used only to rank within the surviving HLS-compliant cell set.

This mirrors the methodology in **Wueller et al. 2026** (JGR Planets, [doi:10.1029/2025JE009434](https://doi.org/10.1029/2025JE009434)), which catalogued 130 candidate Artemis-III sites with the same per-region + HLS-filter approach. The two approaches use different soft-criterion stacks but the same outer framing.

The legacy `selene rank` command is preserved for continuity with the v1.0.0–v1.2.0 release history; both ship in v1.3.0+. Per-region ranking is the primary path.

### HLS hard-constraint filters (v1.3.0)

The four NASA HLS thresholds applied as a multiplicative AND inside each USGS polygon:

| Filter | Threshold | Source data |
| --- | --- | --- |
| Slope at the landing pad | $\le 8°$ | LOLA-derived slope grid |
| Distance to nearest cell with slope $> 8°$ | $\ge 100\,$m | Slope grid + `scipy.ndimage.distance_transform_edt` |
| Direct solar illumination | $\ge 33\,\%$ | Mazarico average-illumination raster |
| Direct-to-Earth (DTE) visibility over the libration cycle | $\ge 50\,\%$ | The v1.1.0 LOS-to-Earth visibility raster |

These thresholds are NASA-published values. Sources:

- NASA HLS specification (NASA 2019).
- Gracy & Lee 2024, *Update on the Artemis III Reference Mission*, LPSC Abstract #1695.
- Wueller, F., et al. 2026, JGR Planets, doi:10.1029/2025JE009434.

They are not tuneable from the validation result; the ranker accepts them as input parameters but the defaults reflect the published NASA values verbatim. A site that fails any one of the four filters is disqualified regardless of how well it scores on the soft criteria.

The 100 m buffer is computed *globally* on the slope grid (the distance from a cell to the nearest cell with slope $> 8°$ doesn't depend on which polygon owns the cell) and then masked to each polygon. Within each polygon, the surviving HLS-compliant cells are ranked by aggregate score, and up to `n_per_region` sites are NMS-extracted at a default 2 km separation.

### The spatial-coupling criterion (v0.7)

`criteria/coupling.py` is the structural fix the v0.6 diagnostic identified. It scores cells by *joint proximity* to two distinct features:

1. **Distance to the nearest PSR**: derived from the Mazarico illumination raster as `illumination < 0.001`, then `scipy.ndimage.distance_transform_edt` with explicit pixel sampling.
2. **Distance to the nearest sunlit ridge**: a cell qualifies if `illumination >= 0.70` AND `5° <= slope <= 25°` — the geometry of polar crater rims (steeper than plains, not cliff-like, well-sunlit). Same distance transform.

The score is the **product** of two linear distance falloffs:

$$ s = \max\!\left(0,\, 1-\tfrac{d_{\text{PSR}}}{d_c}\right)\cdot \max\!\left(0,\, 1-\tfrac{d_{\text{ridge}}}{d_c}\right),\quad d_c = 5\,\text{km}. $$

The product (not sum) is the conjunction. Failing either falloff drives the score to zero — exactly the structural property a linear weighted sum cannot encode. A cell deep inside a far-side PSR (PSR distance 0, ridge distance 200 km) scores 0; a cell on a sunlit rim 30 km from the nearest PSR also scores 0; the rim cell adjacent to a PSR-floor (both distances < 5 km) scores high.

The criterion produces a **sparse mask**: only 0.12 % of finite cells exceed 0.0; only 0.07 % exceed 0.1. The polar rim band shows up clearly; the rest of the cap is essentially black:

![Coupling score, log-scaled — sparse rim band](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/coupling_score.png)

The single tuning knob is `coupling_distance_km`. `selene coupling-sweep` runs the validation pipeline at 1–20 km in 8 steps and plots NASA-region alignment as a function of that knob:

![Spatial-coupling distance sweep vs alignment](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/coupling_distance_sweep.png)

| coupling_distance_km | regions matched within 25 km |
| ---: | ---: |
| 1 – 10 | 0 / 9 |
| 15 – 20 | 2 / 9 |

The curve is monotone-but-flat: tightening the cap below 10 km matches no NASA regions; loosening it to 15 km picks up the Slater Plain / de Gerlache Rim 2 pair (the two NASA candidates already nearest the high-coupling band) and stays at 2/9 through 20 km. **The 15 km threshold equals the validation disk radius** — that's the geometric fingerprint of the validation-metric finding: the criterion can only "match" a NASA region when the matching tolerance equals the disk approximation, because the sites it picks fall on the rim band, not inside the centroid disks. The default 5 km is kept; the criterion is doing what its math says.

### The Earth line-of-sight criterion (v1.1.0)

`criteria/los_to_earth.py` is the seventh criterion, derived from the already-cached LOLA elevation grid in two passes:

1. **`derive_horizon_profile`** ray-marches outward from each pixel in 36 azimuthal directions (10° resolution) up to 100 km, tracking the maximum elevation angle of obstructing terrain (with curvature correction for the lunar sphere of $R = 1737.4\,$km). Sampling is log-spaced (50 distances from 240 m to 100 km) and the result is a 3D ``(azimuth, y, x)`` field cached as a compressed numpy archive (`lola_horizon_profile_southpole_240m.npz`, ~840 MB; netCDF would be more idiomatic but the only backend available without extra system deps doesn't support compression). One-time cost: ~6 minutes single-threaded on the full 2533×2533 grid.

2. **`compute_earth_visibility_fraction`** samples 24 parametric points on the libration ellipse (Earth's sub-Earth point cycles within $\pm 6.5°$ in latitude and $\pm 7.9°$ in longitude over ~27 days), computes Earth's apparent elevation and azimuth at every pixel for each sample, rotates the geographic azimuth into the grid frame using the local grid convergence $\gamma = \mathrm{atan2}(x_p, y_p)$ for south polar stereographic, and counts the fraction of samples for which Earth's elevation exceeds the horizon angle in the matching azimuth bucket. The result ships as a 2D COG (`los_visibility_fraction_southpole_240m.tif`).

The score is a linear ramp anchored to operational comms thresholds: visibility below `min_visibility = 0.20` scores 0.0 (Apollo-era surface ops baselined ``>20%`` direct-comms duty cycle as a crew-safety floor); visibility at and above `target_visibility = 0.50` scores 1.0 (sustained habitat with redundant relay backup); linear in between. **These thresholds are physics-and-operations driven and were chosen *before* the validation rerun, so they are not validation chasing.** The libration sampling at 24 parametric points is a coarse approximation of the Lissajous trajectory traced by physical (period 27.55 d) and optical (27.32 d) libration; documented as such in the function docstrings.

The geometric pattern in the visibility raster is the expected one: deep crater floors near the pole score near 0 (the rim blocks Earth at every libration phase); high ridge tops score near 1 (clear horizon); the geometric south pole pixel scores ~0.5 (Earth oscillates symmetrically above/below local horizontal as the libration cycle traces $\pm 6.5°$ in sub-Earth latitude). Mean visibility at bottom-decile elevation cells (crater floors) is 0.077; mean at top-decile elevation cells (ridge tops) is 0.625 — an 8× spread that the criterion converts directly into a score gradient.

## Engineering decisions

A few choices in the pipeline are worth surfacing because they materially affect what runs, when, and how reliably:

- **Verified URLs over guessed URLs.** Every dataset's download path was manually verified by browsing the PDS Geosciences directory listing rather than synthesised from a documentation pattern. Two URLs (LEND, Watters scarps) remain TODO-flagged because no verified path exists, rather than silently 404'ing.
- **PDS3 + PDS4 in the same loader namespace.** LOLA is PDS3 (open via the detached `.lbl` label, GDAL's PDS driver requirement); Diviner PRP is PDS4 (parse `.xml` schema, fixed-width ASCII table). Both feed the same `xr.DataArray` interface downstream so callers never see the format difference.
- **Triangle mesh → raster rasterisation.** Diviner PRP's 2.88 M triangular-mesh facets are interpolated to the 240 m polar stereographic grid via `scipy.interpolate.griddata` (linear for temperatures, nearest for the discontinuous ice-depth field) — not rasterised by polygon fill, which would be both slower and physically meaningless for sparse triangle centres.
- **Cloud-Optimized GeoTIFF caching with overviews.** Every reprojected raster is written as a COG with internal tiling and DEFLATE compression, so re-running `selene score` after a fresh clone is essentially free.
- **Two independent sensitivity sweeps.** Latin-hypercube weight-vector sweep (200 samples on the simplex) characterises robustness to weight choice; coupling-distance sweep (8 points across 1–20 km) characterises robustness to the single tunable knob in the spatial-coupling criterion. Both ship as CLI subcommands.
- **Three complementary validation metrics under one command.** `selene validate` reports the legacy centroid-distance number, the v1.0 polygon-inside number against 15 km disks, and the v1.2 polygon-inside number against authoritative USGS polygons side by side, so a single run produces the full geometric story rather than forcing a choice between primitives.
- **CI smoke test on bundled sample data.** A separate CI job downloads the ~12 MB sample tarball and runs `preprocess → score → rank → validate → compare` end-to-end on every push to `main`. Catches integration regressions unit tests miss.
- **Cache horizon-profile as compressed numpy, not GeoTIFF.** The v1.1 LOS criterion needs a 3D ``(azimuth, y, x)`` horizon field that doesn't fit a single-band GeoTIFF cleanly; netCDF is the idiomatic xarray choice but the only netCDF backend available without compiled `netCDF4`/`h5netcdf` system deps is `scipy`, and that backend doesn't accept zlib compression — leaving an uncompressed ~930 MB file. Switching to `np.savez_compressed` keeps the cache pure-numpy, compresses the same float32 grid to ~840 MB without new dependencies, and the consumer (`compute_earth_visibility_fraction`) doesn't need the rio metadata anyway. The 2D consumer-facing artifact (`los_visibility_fraction_southpole_240m.tif`) stays as a normal COG.
- **Authoritative validation reference, bundled in-repo.** Validation against NASA's Artemis IV (formerly Artemis III) candidate regions uses USGS's officially-published simplified region envelopes ([DOI 10.5066/P1MEQ6UK](https://doi.org/10.5066/P1MEQ6UK)), not synthesised disk approximations. The dataset shipped late in development; v1.0.0 and v1.1.0 used 15 km-radius disks, which we found systematically misrepresent the actual region geometries (most are ~400 km² quadrilaterals; Mons Mouton Plateau alone is 4452 km², 6× larger than the disk; one disk centroid for "Slater Plain" sits ~180° in longitude away from the USGS polygon's actual location). v1.2.0 ships the USGS GeoJSON in [`src/selene_base/validation/data/`](src/selene_base/validation/data/) so the validation primitive is reproducible without re-downloading external data, and the disk metrics are kept in parallel for continuity with the v1.0 / v1.1 validation history.
- **Per-region search-space framing.** v1.3.0's `selene rank-per-region` searches *within* each USGS polygon rather than globally, with NASA's published HLS hard-constraint filters applied as a precondition. The two CLI subcommands ship side-by-side: the legacy `selene rank` produces a global top-N for continuity; `selene rank-per-region` produces the NASA-aligned per-polygon catalog. The choice is deliberately surfaced rather than hidden behind a flag — they answer two genuinely different questions ("globally most habitat-suitable cells" vs "best HLS-compliant cells *within each NASA candidate region*"), and the v1.0–v1.2 release history is the diagnostic of why those questions diverge.
- **CC-BY-4.0 third-party data, attributed inline.** The Wueller 2026 130-site catalog ships in-repo at [`src/selene_base/validation/data/wueller_2026/LandingSites.shp`](src/selene_base/validation/data/wueller_2026/) under the original CC-BY 4.0 license, with a [README](src/selene_base/validation/data/wueller_2026/README.md) carrying the full author list, both the paper DOI ([10.1029/2025JE009434](https://doi.org/10.1029/2025JE009434)) and the Zenodo deposit DOI ([10.5281/zenodo.17084058](https://doi.org/10.5281/zenodo.17084058)), and the licence text. The bundle is the small (~60 KB) point-catalog sidecar of the deposit only; the 884 MB HLS slope raster from the same deposit is not bundled. The legacy synthetic-placeholder CSV from v1.4.0 is retained for backward compatibility behind a `DeprecationWarning`.

## Validation

> **Detailed catalog (v1.5):** every site, every per-criterion score, every Wueller comparison — [`docs/v1.5_catalog_report.md`](docs/v1.5_catalog_report.md). Companion to this README: per-site breakdown of the 69-site v1.5 catalog and cross-evaluation of Wueller's 73 in-scope sites against selene's seven criteria.

### Per-region HLS-compliant catalog (v1.3.0)

`selene rank-per-region` followed by `selene validate-per-region` produces the NASA-aligned per-polygon catalog: for each USGS polygon, the cells passing every HLS hard filter are ranked by aggregate suitability score, and up to ``n_per_region`` sites (default 10) are NMS-extracted at 2 km separation.

![Per-region HLS-eligibility maps for all nine USGS polygons](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/per_region_eligibility.png)

**Result (default weights, default HLS thresholds, n_per_region = 10):**

| USGS region | code | n sites | best score | mean score | HLS-eligible area | eligible % |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| **Mons Mouton Plateau** | MP | 10 | **0.746** | 0.734 | ~671 km² | **15.07 %** |
| Nobile Rim 2 | N2 | 10 | 0.720 | 0.698 | ~33.8 km² | 8.46 % |
| Mons Mouton | MM | 10 | 0.710 | 0.675 | ~17.7 km² | 6.92 % |
| Haworth | HW | 10 | 0.703 | 0.682 | ~65.7 km² | 7.40 % |
| de Gerlache Rim 2 | G2 | 2 | 0.703 | 0.680 | ~0.4 km² | **0.10 %** |
| Slater Plain | SP | 9 | 0.707 | 0.669 | ~10.5 km² | 2.62 % |
| Peak Near Cabeus B | CB | 9 | 0.705 | 0.666 | ~9.0 km² | 2.25 % |
| Nobile Rim 1 | N1 | 10 | 0.697 | 0.676 | ~25.9 km² | 6.47 % |
| **Malapert Massif** | MA | **0** | — | — | 0 km² | **0.00 %** |

[![Per-region HLS-compliant ranking: 70 sites color-coded by USGS region across 8 of 9 Artemis IV candidate regions](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/per_region_ranking.png)](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/per_region_ranking.png)

**Headline numbers:**

- **70 total sites across 8 / 9 USGS regions** (NMS-at-2km caps the small regions: Slater Plain and Peak Near Cabeus B saturate at 9; de Gerlache Rim 2 caps at 2).
- **Malapert Massif has zero HLS-compliant cells.** The polygon's terrain (per the LOLA + Mazarico + Diviner stack) does not contain a single 240 m cell that simultaneously satisfies slope ≤ 8°, distance ≥ 100 m to steeper cells, illumination ≥ 33 %, and DTE visibility ≥ 50 %. This is a real terrain-driven finding, consistent with NASA's selection presumably relying on higher-resolution NAC stereo DEMs to characterise this region.
- **Mons Mouton Plateau** is the most "easy" region: 15.07 % of its polygon-cells satisfy every HLS filter (~671 km² of HLS-eligible area inside the 4452 km² polygon), and the best HLS-compliant site there scores **0.746** — the highest in the catalog.
- **de Gerlache Rim 2** is the most constrained of the regions with sites: only 0.10 % of polygon cells (~0.4 km²) are HLS-eligible, and only 2 sites fit at the 2 km NMS separation — independent of `n_per_region`.
- Globally, **9.44 %** of the 240 m polar grid (~605 600 of 6.4 M cells) satisfies every HLS hard filter — most of that is outside any USGS polygon, on the polar rim band the v1.0–v1.2 diagnostic releases identified.

The v1.3 sites are guaranteed inside USGS polygons and HLS-compliant by construction. The relevant scientific question — "do they identify the same cells NASA's process identifies?" — is answered quantitatively in v1.4.2 against **Wueller et al. 2026's 130 published sites**: 56/70 selene sites (80 %) match within 5 km of an in-scope Wueller site, median 1.88 km. See "Quantitative comparison against Wueller et al. 2026" below.

### Global-ranking validation (legacy, v1.0.0 — v1.2.0)

`selene validate` compares the *global* top-N ranked sites (from `data/outputs/top_sites.geojson`) against three references in [`src/selene_base/validation/nasa_regions.py`](src/selene_base/validation/nasa_regions.py):

1. **NASA centroids** as 15 km-radius disks (legacy, v0.4–v1.1.0). Centroids from NASA's October 2024 Artemis III site-selection announcement; the disk radius is the publicly cited "operational region" scale. **Not authoritative geometry** — used through v1.1.0 only because the actual polygons were not openly published.
2. **Disk inside/outside** of those same 15 km disks (v1.0.0 polygon-inside metric).
3. **USGS published polygons** (v1.2.0 headline). The official simplified region envelopes from USGS Data Release [10.5066/P1MEQ6UK](https://doi.org/10.5066/P1MEQ6UK) (McClernan 2024), bundled at [`src/selene_base/validation/data/nasa_regions_polygons_usgs.geojson`](src/selene_base/validation/data/nasa_regions_polygons_usgs.geojson). 4-vertex quadrilaterals in lunar planetocentric lon/lat space, sourced from NASA's LROC QuickMap region definitions. These are simplified envelopes, not the full operational landing footprints, but they are the **authoritative public approximation** of NASA's selected region geometries.

The USGS polygons differ substantively from the disk approximations:

- **Names differ**: USGS calls one region "Peak Near Cabeus B" (centred on the rim peak at lat -83.7°, lon -68.7°), not "Cabeus B" (centred on the crater floor at lat -82.3°, lon -53.3° in the legacy disk list — about 150 km away from the USGS polygon).
- **Sizes vary**: the simplified envelopes total ~8000 km² (vs the disks' uniform 707 km² × 9 = 6362 km²), but Mons Mouton Plateau alone is **4452 km² — over 6× the disk area**, while seven other regions are ~400 km² (smaller than the disks).
- **Locations are not always close**: the legacy disk centroid for "Slater Plain" sits at lon -54.3°, ~180° away from where USGS publishes Slater Plain (lon +125°). The legacy list inherited a press-release centroid that doesn't match the authoritative USGS polygon.

Three metrics for each top site:

1. **Within X km of any centroid** (legacy v0.4)
2. **Inside any 15 km disk** (v1.0.0)
3. **Inside any USGS polygon** (v1.2.0 — *headline through v1.2.0*)

And two for each USGS region:

1. **Distance to nearest top-N site** — distance from the USGS polygon boundary to the closest selene-base candidate.
2. **Contains a top-N site** — is at least one selene-base candidate inside the USGS polygon?

![Score map with USGS polygons (red) and 15 km disks (grey, dashed) overlaid](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/usgs_polygon_validation.png)

### Per-USGS-region results (v1.2.0: 7-criterion + LOS-to-Earth, USGS polygons)

| USGS region | code | area (km²) | nearest top site | distance to polygon (km) | contains a top site? |
| --- | --- | ---: | --- | ---: | --- |
| de Gerlache Rim 2 | G2 | 400.0 | site_01 | **41.5** | no |
| Peak Near Cabeus B | CB | 399.8 | site_18 | 67.3 | no |
| Haworth | HW | 888.0 | site_01 | 71.5 | no |
| Mons Mouton Plateau | MP | 4451.8 | site_17 | 72.4 | no |
| Slater Plain | SP | 399.8 | site_01 | 75.2 | no |
| Malapert Massif | MA | 441.0 | site_01 | 102.1 | no |
| Nobile Rim 1 | N1 | 400.0 | site_01 | 115.6 | no |
| Mons Mouton | MM | 256.0 | site_01 | 118.7 | no |
| Nobile Rim 2 | N2 | 400.0 | site_07 | 134.8 | no |

The closest USGS polygon is **de Gerlache Rim 2 at 41.5 km** from `site_01` (-89.7°, +17.7°). The median distance to the nearest USGS polygon across the top-20 is **135.1 km**. Every top site sits outside every USGS polygon under global ranking — the geometric separation persists against the authoritative reference. v1.3.0's per-region reframing resolves this by searching inside the polygons rather than globally.

### Per-region results vs the 15 km disks (legacy, v1.0.0 — v1.1.0)

The disk-based table is preserved for continuity with the v1.0.0 / v1.1.0 history; the v1.2 USGS table above is the authoritative metric.

| NASA disk | nearest site | dist to centroid (km) | dist to disk edge (km) | inside disk? |
| --- | --- | ---: | ---: | --- |
| Cabeus B | site_18 | **27.3** | **12.3** | no |
| Haworth | site_01 | 97.6 | 82.6 | no |
| Malapert Massif | site_01 | 115.4 | 100.4 | no |
| Mons Mouton | site_08 | 59.1 | 44.1 | no |
| Mons Mouton Plateau | site_08 | 85.1 | 70.1 | no |
| Nobile Rim 1 | site_01 | 127.6 | 112.6 | no |
| Nobile Rim 2 | site_17 | 130.3 | 115.3 | no |
| de Gerlache Rim 2 | site_01 | 47.8 | 32.8 | no |
| Slater Plain | site_01 | 55.4 | 40.4 | no |

The disk-based "Cabeus B" entry shows a top site within 12.3 km of the disk edge — but the corresponding USGS polygon ("Peak Near Cabeus B", at a different geographic location) is 67.3 km from the same `site_18`, because the legacy disk centroid sits ~150 km north-east of where USGS places the actual region. Two of the disk-table closest distances (Cabeus B and de Gerlache Rim 2) are misleading once measured against the right geometry.

## Quantitative comparison against Wueller et al. 2026

[Wueller, F., et al. (2026)](https://doi.org/10.1029/2025JE009434) published in *Journal of Geophysical Research: Planets* a peer-reviewed analysis identifying **130 candidate Artemis III landing sites** within NASA's candidate regions using essentially the same outer methodology selene-base implements: NASA HLS hard filters (slope < 8°, ≥ 100 m buffer to steeper terrain) followed by within-region selection. v1.4.2 ships the **quantitative comparison** against the authors' Zenodo data deposit ([doi:10.5281/zenodo.17084058](https://doi.org/10.5281/zenodo.17084058), CC-BY 4.0, 130-site shapefile bundled in-repo); v1.5 reruns the pipeline at 20 m and v1.8 activates the eighth criterion on the same catalog.

> **Current headline (v1.8, eight criteria active, 20 m per-region tiled):** 56 / 69 selene sites (81.2 %) match within 5 km of an in-scope Wueller site; 46 / 73 (63.0 %) Wueller-side; **median matched-pair distance 1.69 km**. Tracked through the per-release subsections below: v1.4.2 baseline (240 m, seven criteria) → v1.5 (20 m, seven criteria) → v1.8 (20 m, eight criteria).

### Headline at v1.4.2 (240 m baseline, seven criteria, default weights, threshold = 5 km, in-scope only)

| metric | value |
| --- | ---: |
| selene per-region sites | 70 |
| Wueller 2026 sites (total catalog) | 130 |
| Wueller sites in USGS scope (NASA Oct 2024 nine) | **73** |
| Wueller sites out of USGS scope | 57 |
| **selene sites matched within 5 km of an in-scope Wueller site** | **56 / 70 (80 %)** |
| Wueller in-scope sites matched within 5 km of a selene site | 46 / 73 (63 %) |
| **median matched-pair distance** | **1.88 km** |
| match threshold (regional-granularity scale) | 5 km |

**Fifty-six of seventy selene sites land within 5 km of a peer-reviewed candidate site identified by an independent group using the same outer methodology on a higher-resolution DEM.** Median match distance is 1.88 km — well inside the 1–5 km regional-granularity scale at which candidate-site *selection* operates (NASA HLS landing accuracy is 100 m, but selecting *which* terrain to land on is the kilometre-scale question). At v1.4.2's `n_per_region = 10` default, both directions of the comparison are now near-saturated: 80 % of selene sites match a Wueller site, and 63 % of in-scope Wueller sites have a selene neighbour within 5 km — up from 41 % at v1.4.1's `n_per_region = 3`.

### Per-region agreement

Selene ranks up to 10 sites per region (`n_per_region = 10` is the v1.4.2 default; NMS-at-2km caps the actual count below 10 for small polygons). Wueller publishes 5–11 sites per region. The "matched" column counts selene sites whose nearest **same-region** Wueller site is within 5 km.

| USGS region | selene n | Wueller in-scope n | matched | median match dist (km) |
| --- | ---: | ---: | ---: | ---: |
| Haworth | 10 | 11 | 7 / 10 | 2.09 |
| Mons Mouton | 10 | 10 | **10 / 10** | 2.30 |
| Mons Mouton Plateau | 10 | 11 | 5 / 10 | 2.38 |
| Nobile Rim 1 | 10 | 9 | 8 / 10 | 1.30 |
| Nobile Rim 2 | 10 | 9 | 9 / 10 | 1.14 |
| Peak Near Cabeus B | 9 | 5 | **9 / 9** | 1.94 |
| Slater Plain | 9 | 11 | 6 / 9 | 2.16 |
| de Gerlache Rim 2 | 2 | 7 | 0 / 2 | — |
| **Total** | **70** | **73** | **56 / 70** | **1.88** |

**Two regions agree at 100 % match-within-region** (Mons Mouton, Peak Near Cabeus B); Nobile Rim 1 and Nobile Rim 2 reach 80 %+ at the new default. The structural outliers are unchanged from v1.4.1:

- **de Gerlache Rim 2 (0/2 matched).** selene's two HLS-compliant sites are 8.8 km from the nearest Wueller dGR2 site — outside the 5 km threshold but still within an order of magnitude. This is the only region where selene caps below `n_per_region` regardless of the cap (the polygon contains only 2 cells passing every HLS filter at 240 m resolution at 2 km NMS).
- **Mons Mouton Plateau (5/10 matched), Slater Plain (6/9 matched), Haworth (7/10 matched).** The largest disagreements are still 5.6–8.3 km — *just* outside the threshold. At a 6 km threshold most flip to match.

The longest selene-only nearest distance at the new default is 12.5 km (a Mons Mouton Plateau outlier vs Wueller MMP09). The longest Wueller-only nearest distance is 18.6 km (dGR202 vs selene's two dGR2 sites) — the dGR2 terrain divergence between selene's HLS-eligible cells and Wueller's catalog is genuine, not a sample-size artefact.

[![selene-base v1.4.2 sites (cyan) vs Wueller 2026 sites (yellow) with matched-pair connector lines, on USGS Artemis III polygon outlines](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/headline_v141.png)](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/headline_v141.png)

### In-scope vs out-of-scope split

Wueller's 130 sites span 16 region codes; **57 are outside NASA's October 2024 down-selected nine** (Wueller pre-dates the down-selection and retained the earlier 13-region list, plus Amundsen Rim and Mons Malapert which are not in the USGS list). selene-base only ranks within the down-selected nine, so the apples-to-apples comparison is against the 73 in-scope sites. The eight out-of-scope Wueller regions:

| Wueller code | Wueller region | n |
| --- | --- | ---: |
| AR | Amundsen Rim | 10 |
| CR | Connecting Ridge | 9 |
| CRE | Connecting Ridge Extension | 6 |
| FRA | Faustini Rim A | 8 |
| MMA | Mons Malapert (distinct from USGS Malapert Massif) | 5 |
| PNS | Peak Near Shackleton | 6 |
| dGKM | de Gerlache-Kocher Massif | 7 |
| dGR1 | de Gerlache Rim 1 | 6 |

Run `selene compare-wueller --no-filter-to-usgs-scope` to compare against all 130; under that mode selene's 56 still match within 5 km (the 14 unmatched selene sites stay unmatched — none of them have a closer match outside their own USGS region).

### Outputs

`selene compare-wueller` writes:

- `data/outputs/wueller_comparison.json` — full result dict (per-site, per-region, in-scope counts, headline metrics).
- `data/outputs/wueller_comparison.csv` — flat per-site distance table for inspection.
- Stdout summary: headline counts (with in-scope/out-of-scope split), per-region agreement table, notable disagreements.

Distances are computed in lunar south-polar stereographic metres (`+proj=stere +lat_0=-90 +lat_ts=-90 +R=1737400`), conformal at the pole and sub-percent error vs great-circle for sub-100 km offsets.

### Data provenance

The catalog ships in-repo at [`src/selene_base/validation/data/wueller_2026/`](src/selene_base/validation/data/wueller_2026/) as a six-file shapefile bundle (`LandingSites.shp/.shx/.dbf/.prj/.cpg`, ~60 KB total) plus a [README](src/selene_base/validation/data/wueller_2026/README.md) carrying the full attribution. Source is the authors' Zenodo deposit "Complementary Data for Wueller et al. (2026)" (Wueller, Berger, Christopher, Sugimoto, Thaker, Carton, Jo, Lee, Pedrelli, Sanchez, & Kring, 2025), [doi:10.5281/zenodo.17084058](https://doi.org/10.5281/zenodo.17084058), licensed CC-BY 4.0. The deposit also contains an 884 MB HLS slope raster (`HLS_LandingAreas_8°_-100m_buffer.tif`), exploration-area polygons, and an illumination-modeling layer set; only the 130-site point catalog is bundled in-repo, the rest stays in the upstream archive and may be integrated in a future release.

### What ships in v1.4.2

- **Default `n_per_region` raised from 3 → 10** in [`scoring/ranking.py`](src/selene_base/scoring/ranking.py), [`pipeline/rank_per_region.py`](src/selene_base/pipeline/rank_per_region.py), and the `selene rank-per-region` CLI option. NMS-at-2km caps the small regions automatically (Slater Plain, Peak Near Cabeus B at 9; de Gerlache Rim 2 at 2). Total sites at the new default: 70.
- [`src/selene_base/validation/data/wueller_2026/`](src/selene_base/validation/data/wueller_2026/) — bundled shapefile and README (v1.4.1, unchanged).
- [`src/selene_base/validation/wueller_comparison.py`](src/selene_base/validation/wueller_comparison.py) — `load_wueller_sites` (defaults to the shapefile, falls back to the legacy synthetic CSV with a `DeprecationWarning`), `compare_sites` (with `filter_to_usgs_scope` parameter, default True), `WUELLER_TO_USGS_REGION_MAP`, `WUELLER_CODE_TO_NAME`, `is_synthetic_placeholder`, `render_summary`. Distance computation in lunar polar stereographic metres (conformal at the pole, sub-percent error vs great-circle for sub-100 km offsets — verified by an explicit known-offset unit test).
- [`src/selene_base/pipeline/compare_wueller.py`](src/selene_base/pipeline/compare_wueller.py) and the [`selene compare-wueller`](src/selene_base/cli.py) CLI subcommand, with `--filter-to-usgs-scope/--no-filter-to-usgs-scope` and a backward-compat `--wueller-csv` alias.
- [`tests/test_wueller_comparison.py`](tests/test_wueller_comparison.py) — 20 tests passing, with the real-data assertion bound updated from `<= 23` to `<= 70` for the new default.
- [`notebooks/08_wueller_comparison.py`](notebooks/08_wueller_comparison.py) and [`notebooks/09_headline_v141.py`](notebooks/09_headline_v141.py) — regenerate the headline overlay map ([docs/img/selene_vs_wueller.png](docs/img/selene_vs_wueller.png)), the per-region dot-and-line diagnostic ([docs/img/headline_v141.png](docs/img/headline_v141.png)), the distance histogram ([docs/img/wueller_distance_hist.png](docs/img/wueller_distance_hist.png)), and the per-region match-count bar chart ([docs/img/wueller_per_region_bars.png](docs/img/wueller_per_region_bars.png)) at the new default.

### Resolution analysis (v1.5 — 20 m Wueller-class)

v1.5 reruns the per-region HLS-filtered ranking on the 20 m LOLA DEM (`ldem_80s_20m`, 30400 × 30400 polar-stereo), the same native resolution Wueller et al. 2026 use. The global 240 m grid is not feasible at 20 m (~132 GB for the LOS horizon profile alone), so v1.5 ships a per-region tiled driver: each USGS polygon's bounding box plus a 100 km horizon buffer is windowed and reprojected to a local 20 m grid, the horizon profile and Earth-LOS visibility are derived on the GPU (NVIDIA GB10 / CuPy), and the four HLS hard filters run at 20 m inside the polygon. The 240 m run remains the default; v1.5 is opt-in via `--tiled-per-region --resolution 20`.

| metric | v1.4.2 (240 m) | **v1.5 (20 m)** | Δ |
| --- | ---: | ---: | ---: |
| selene per-region sites | 70 | **69** | -1 |
| selene matched within 5 km of an in-scope Wueller site | 56 / 70 (80.0 %) | **56 / 69 (81.2 %)** | +1.2 pp |
| in-scope Wueller sites matched within 5 km of a selene site | 46 / 73 (63.0 %) | 44 / 73 (60.3 %) | -2.7 pp |
| **median matched-pair distance** | 1.88 km | **1.76 km** | **-6.5 %** |
| max matched-pair distance | 4.99 km | 4.81 km | -3.4 % |

**The methodology converges at peer-reviewed resolution.** The median matched-pair distance tightens by 6.5 % (1.88 km → 1.76 km) and the maximum tightens by 3.4 % at the 5 km threshold; selene → Wueller match rate is essentially flat (80.0 % → 81.2 %, well within the ±1-site noise band). The Wueller → selene rate slips by 2.7 pp (46/73 → 44/73) — selene drops one site overall (70 → 69 because Mons Mouton Plateau, Slater Plain, and Peak Near Cabeus B re-rank slightly under the finer slope/buffer), and the new top-10 in those three regions covers two fewer Wueller targets within the 5 km radius. None of those three flip *out* of agreement; they shift to a different cell within the same polygon. The structural outlier from v1.4.2 — de Gerlache Rim 2 with only 2 HLS-compliant cells, 0 within 5 km of Wueller's dGR2 cluster — persists: at 20 m the polygon yields the same 2 cells, *not* a new HLS-eligible band the 240 m grid was hiding. **The resolution-induced disagreement at de Gerlache Rim 2 is genuine terrain, not a 240 m sampling artefact.**

[![selene-base v1.5 (20 m, per-region tiled): 69 sites, 81 % within 5 km of Wueller 2026 (median 1.76 km), with USGS polygon outlines and matched-pair connectors](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/selene_vs_wueller_20m.png)](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/selene_vs_wueller_20m.png)

Per-region detail (240 m vs 20 m, side-by-side):

[![Per-region selene-matched percentage and median matched-pair distance, v1.4.2 (240 m, slate) vs v1.5 (20 m, teal)](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/resolution_sensitivity_v15.png)](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/resolution_sensitivity_v15.png)

| USGS region | selene n (v1.5) | Wueller in-scope n | matched (v1.5) | median match dist v1.4.2 → v1.5 (km) |
| --- | ---: | ---: | ---: | ---: |
| Haworth | 10 | 11 | 8 / 10 | 2.09 → **1.77** |
| Mons Mouton | 10 | 10 | **10 / 10** | 2.30 → **1.91** |
| Mons Mouton Plateau | 10 | 11 | 4 / 10 | 2.38 → **2.03** |
| Nobile Rim 1 | 10 | 9 | 8 / 10 | 1.30 → 1.32 |
| Nobile Rim 2 | 10 | 9 | 9 / 10 | 1.14 → **1.13** |
| Peak Near Cabeus B | 9 | 5 | **9 / 9** | 1.94 → 2.04 |
| Slater Plain | 8 | 11 | 5 / 8 | 2.16 → **2.03** |
| de Gerlache Rim 2 | 2 | 7 | 0 / 2 | — → — |
| Malapert Massif | 0 | 0 | — | — |
| **Total** | **69** | **73** | **56 / 69** | **1.88 → 1.76** |

Six of eight regions tighten or hold their median matched-pair distance at 20 m. Nobile Rim 1 ticks +0.02 km (within rounding noise of the n=8 sample) and Peak Near Cabeus B drifts +0.10 km — selene's CB sites are picking the same rim cluster but the 5 nearest Wueller CB sites are slightly off-axis, an effect the coarser 240 m grid was averaging away. **Malapert Massif still has zero HLS-compliant cells at 20 m** (vs zero at 240 m); the polygon's terrain is genuinely steep, not a coarse-resolution effect.

The v1.5 driver is opt-in. The 240 m baseline still runs in ~6 min wall-clock end-to-end on a modest workstation; the v1.5 20 m run requires a CUDA-capable GPU (developed and benchmarked on an NVIDIA GB10) and ~140 GB of unified-memory headroom because each per-region tile carries a ~14–22 GB float32 horizon profile through the GPU. End-to-end time on the reference DGX Spark host is ~32 min for `preprocess --tiled-per-region` (all 9 polygons) plus ~22 min for `rank-per-region --tiled-per-region`. CPU-only fallback works on the test cases but is impractically slow for the full 9 regions.

#### What ships in v1.5

- [`src/selene_base/criteria/los_to_earth.py`](src/selene_base/criteria/los_to_earth.py) — `derive_horizon_profile` and `compute_earth_visibility_fraction` gain a `use_gpu` flag; the inner ray-march delegates to `cupyx.scipy.ndimage.map_coordinates` with one block per pixel, falling back transparently to the existing CPU path when CuPy is not importable.
- [`src/selene_base/pipeline/preprocess_tiled.py`](src/selene_base/pipeline/preprocess_tiled.py) — per-USGS-polygon tiled driver. Each tile is the polygon's bbox padded by `DEFAULT_BUFFER_M = 100_000` (matching `los_to_earth.DEFAULT_MAX_HORIZON_KM`, so the ray-march sees the same physical horizon as the 240 m run did). Tiles are processed sequentially with explicit cupy memory-pool drain between regions.
- [`src/selene_base/pipeline/rank_per_region_tiled.py`](src/selene_base/pipeline/rank_per_region_tiled.py) — at-resolution HLS filtering. The 240 m aggregate score raster is upsampled onto each tile via bilinear `reproject_match`; slope and Earth-LOS visibility are *re-derived* at 20 m on the tile grid; HLS distance-to-steep is computed as a 20 m EDT (so the 100 m HLS buffer is now a 5-pixel constraint instead of the 0.42-pixel approximation it was at 240 m); greedy NMS at 2 km within each polygon picks the top-10.
- [`tests/test_los_to_earth_gpu.py`](tests/test_los_to_earth_gpu.py), [`tests/test_preprocess_tiled.py`](tests/test_preprocess_tiled.py), [`tests/test_rank_per_region_tiled.py`](tests/test_rank_per_region_tiled.py) — GPU-path correctness against the CPU reference (mark `cuda` skips on CI), tile-spec geometry, NPZ cache layout, and end-to-end pipeline smoke on a synthetic small region.
- CLI: `selene preprocess --tiled-per-region --resolution 20`, `selene rank-per-region --tiled-per-region --resolution 20`, both with `--region-code` for single-tile reruns and standard `--overwrite` semantics on the per-region NPZ caches.
- [`notebooks/10_resolution_v15.py`](notebooks/10_resolution_v15.py) — regenerates [docs/img/selene_vs_wueller_20m.png](docs/img/selene_vs_wueller_20m.png) and [docs/img/resolution_sensitivity_v15.png](docs/img/resolution_sensitivity_v15.png).
- Run artefacts: `data/outputs/per_region_tiled/sites.{geojson,csv}`, `data/outputs/per_region_tiled/per_region_summary.json`, `data/outputs/v15/wueller_comparison.{json,csv}`.

## Robustness

Anyone reading 0/20 fairly asks: *is that just a function of the default weights?* Run `selene sensitivity --n-samples 200` to find out: it draws 200 weight vectors via Latin hypercube over the active-criteria simplex, runs `aggregate -> top_n_sites -> proximity_analysis` for each, and reports the distribution of "NASA regions matched within 25 km" alongside the default-weight result.

![sensitivity over 200 weight samples](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/sensitivity_distribution.png)

The 7-criterion sensitivity sweep is run against **both** the legacy 15 km disks and the USGS polygons (v1.2.0). Distribution across 200 weight samples:

- **Centroid-distance metric (legacy)**: 157/200 samples (78.5 %) match 0 regions within 25 km of any centroid; 8/200 match 1, 16/200 match 2, 14/200 match 3, 5/200 match 4 — the sensitivity ceiling lifted from 3/9 in v1.0.0 to 4/9 in v1.1.0–v1.2.0.
- **Disk inside/outside (v1.0.0 metric)**: 21/200 samples (10.5 %) put at least one top site *inside* a 15 km disk. Pre-v1.1.0 (six criteria, no LOS) this was 0/200.
- **USGS polygon inside/outside (v1.2.0 metric)**: **194 / 200 samples (97.0 %) put 0 sites inside any USGS polygon**; 5/200 (2.5 %) put 1 inside; 1/200 (0.5 %) puts 2 inside. **Best case across the sweep is 2 / 20 sites inside USGS polygons.** The minimum median-distance-to-nearest-USGS-polygon across the sweep is 41.2 km; the mean is 128 km.

The best weight regime against the legacy centroid metric is `slope = 0.07, illumination = 0.24, coupling = 0.21, thermal = 0.00, ice = 0.16, hazard = 0.31, los_to_earth = 0.08` — every criterion contributes, and 4/9 NASA regions are now within 25 km of a top site. The default-weight global-ranking headline against the **USGS polygons** is **0/20 inside any polygon, 0/9 USGS regions containing a top site, median distance 135.1 km, closest 41.5 km (de Gerlache Rim 2)**.

The structural picture: against the legacy disk metric the geometric gap collapses for ~10 % of the weight simplex; against the authoritative USGS polygons it collapses for only ~3 % of the simplex (max 2/20 inside, never more). The disk approximations were systematically wrong — names, sizes, and one location (Slater Plain) were misplaced — but the *underlying geometric separation* between the model's rim-band optimum and NASA's authoritative regions persists across both validation references and across most of the weight simplex. The diagnostic comparison shows what the criteria are saying about that gap. v1.3.0 reframes the analysis to per-region ranking and resolves this by searching within the polygons.

## Diagnostic comparison

Run `selene compare` to ask: *at NASA's centroids vs at our top-20 (global ranking), which criteria agree and which disagree?* The table is reordered below by signed delta — agreement first, disagreement last — because that's the actual structure of the result.

![Per-criterion score: where we differ from NASA](https://raw.githubusercontent.com/Alex0420W/selene-base/main/docs/img/comparison.png)

| criterion | our top-20 | NASA 9 centroids | delta | \|t\| | reads as |
| --- | --- | --- | ---: | ---: | --- |
| **hazard** | 0.976 ± 0.027 | 0.969 ± 0.023 | **+0.007** | 0.73 | strong agreement |
| **coupling** | 0.025 ± 0.113 | **0.000 ± 0.000** | **+0.025** | 1.00 | both near zero — see below |
| **ice** | 0.995 ± 0.009 | 0.916 ± 0.096 | **+0.079** | 2.45 | strong agreement |
| illumination | 0.747 ± 0.071 | 0.321 ± 0.274 | +0.426 | 4.59 | major disagreement |
| **thermal** | 0.962 ± 0.086 | 0.526 ± 0.269 | +0.437 | 4.77 | responsive (v1.0.0 correction) |
| **los_to_earth** | **1.000 ± 0.000** | **0.525 ± 0.505** | +0.475 | 2.82 | bimodal at NASA — see below |
| slope | 0.948 ± 0.125 | 0.285 ± 0.288 | +0.663 | 6.64 | major disagreement |

**Of the seven criteria, four agree at most-or-very-strongly with NASA's selection.** Hazard agrees almost identically (delta +0.007). Ice agrees within 0.08. Thermal sits in the responsive band on both sides (0.962 vs 0.526; v1.0.0's corrected target moved every cell out of the Gaussian's tail). Coupling agrees in the sense that both sets are near zero — the criterion identifies the rim band, which neither set sits on (see below).

**LOS-to-Earth is the single most-bimodal NASA-side number.** Our top-20 score 1.000 ± 0.000 — every top site has Earth visibility above 50 % of the libration cycle, by construction of how the criterion ramps and how the rank picks. NASA centroids score 0.525 ± 0.505 — essentially uniform-bimodal: about half of NASA's centroids are deep enough inside their candidate craters to score 0 (Earth never above local horizon) and the other half are on rim-adjacent terrain that scores 1 (Earth always above local horizon at favorable libration). NASA's selection knowingly mixes sites with-and-without direct Earth comms because the actual landing footprints within each region target off-centroid rim cells with comms; the centroid point is just the geometric centre of the disk approximation, not necessarily the operational target. **Like the v0.7 coupling number, the LOS comparison is partly diagnostic of the centroid-as-proxy issue.**

**Coupling is still the most informative single number in the project.** *NASA centroids score exactly 0.000 ± 0.000 on coupling.* Our top-20 score 0.025 — small, but non-zero. The criterion is correctly identifying the rim band where PSR meets sunlit ridge at a 5 km coupling distance — and *neither site set sits on that band*. NASA's centroids are *inside* their candidate regions (a centroid of a 15 km disk is, by construction, in the middle of the disk, which means inside a crater for Cabeus B / Haworth, on a massif for Malapert, on a plain for Slater Plain), 5–15 km off the rim that the actual NASA landing footprints would target. **The 0.000 number is a discovery about the validation metric: distance to centroid is the wrong proxy for "match a NASA candidate region," because NASA's preferred landing sites within each region are off-centroid by construction.** v1.0.0 added the polygon-inside metric specifically to test this — and even with that primitive plus the v1.1.0 LOS criterion, the default-weights global-ranking answer is still 0/20, because the rim band is geometrically *outside* the disks for the operations-driven default weights even when the right primitive is used. v1.3.0's per-region reframing resolves this by searching inside the polygons.

**Slope and illumination remain the major disagreement** because the linear-sum aggregator cannot represent the *spatial coupling* between them — the methodology finding from v0.6 still holds. The v0.7 coupling criterion captures the conjunction directly (the product is the AND), but with only 0.18 default weight it cannot dominate the linear sum across the other six criteria. Pushing coupling weight higher would require either changing the default or — more architecturally honest — switching the aggregator itself to TOPSIS, which penalises lop-sided profiles globally rather than at one specific spatial coupling.

The `|t|` column is a Welch two-sample t-statistic, informational only.

## Architecture

```
selene-base/
├── src/selene_base/
│   ├── data/                # download + load + reproject + rasterize
│   ├── criteria/            # eight [0,1] scoring functions (slope, illumination,
│   │                        #                                 hazard, thermal, ice,
│   │                        #                                 coupling, los_to_earth,
│   │                        #                                 seismic — v1.8)
│   ├── scoring/             # normalize, aggregate (renormalising), ranking (NMS, per-region)
│   ├── validation/          # NASA candidate regions + proximity_analysis + Wueller comparison
│   ├── viz/                 # folium webmap + per-site HTML reports
│   ├── pipeline/            # one orchestrator module per CLI subcommand
│   └── cli.py               # typer CLI: download, preprocess, score, rank, rank-per-region,
│                            #            score-wueller-sites, validate, validate-per-region,
│                            #            viz, compare, sensitivity, coupling-sweep,
│                            #            compare-wueller
├── config/                  # region_southpole.yaml, weights_default.yaml
├── data/                    # raw/ processed/ outputs/ (all gitignored)
├── notebooks/               # jupytext .py scripts; one per week
├── tests/                   # synthetic-data unit tests + skipif-guarded data tests
└── .github/workflows/ci.yml
```

The dependency graph is one-way: `data/` is the foundation; `criteria/` reads loaded rasters; `scoring/` aggregates criterion outputs; `validation/` and `viz/` consume scoring outputs; `pipeline/` orchestrates; `cli.py` exposes the orchestrators. Tests follow the same layering.

**381 tests collected** (372 passed, 9 skipped: LEND + Robbins data not present), ~80 % combined branch coverage, all running synthetically in CI on Python 3.11 and 3.12. Real-data tests are guarded with `pytest.mark.skipif(not Path(...).exists())` so the suite stays green without ~900 MB of cached LRO data. CI runs a separate `pipeline-smoke` job on push to `main` that downloads the bundled ~12 MB sample tarball, runs `preprocess -> score -> rank-per-region -> validate-per-region -> compare`, and asserts every output file is on disk and schema-valid.

## Roadmap

> **Note on mission designation (April 2026).** NASA restructured the Artemis program in February 2026: the first crewed lunar landing at the south pole was reassigned from Artemis III to Artemis IV (target early 2028); Artemis III is now an Earth-orbit test mission and Lunar Gateway was canceled in March 2026. The candidate landing regions identified in NASA's October 2024 down-selection are unchanged and remain the validation reference for selene-base. Roadmap entries below were written under the original Artemis III framing; the methodology and data they describe are unchanged. v1.5.1 reframes the project's forward-looking copy and adds Lawrence 2025's NASA Figures-of-Merit framework as a citation; the catalog itself is identical to v1.5.

- **v0.1 — data acquisition.** ✅ `selene download` for Robbins, LOLA, Mazarico illumination; Diviner / LEND / scarps URLs flagged.
- **v0.2 — common grid + slope criterion.** ✅ `reproject_to_grid`, COG cache, slope criterion end-to-end on real data.
- **v0.3 — full scoring + ranking.** ✅ All six criteria (3 on real data, 3 skip cleanly), KDTree crater density, NMS top-N extraction.
- **v0.4 — validation + visualisation.** ✅ NASA Artemis III proximity comparison, interactive folium web map, per-site HTML reports.
- **v0.5 — robustness, diagnostic, sample data.** ✅ Latin-hypercube weight-sensitivity sweep, per-criterion `selene compare` diagnostic, bundled ~12 MB sample tarball, CI pipeline smoke test on the sample.
- **v0.6 — Diviner Polar Resource Product integration.** ✅ PDS4 character-table parser, triangle-mesh-to-grid rasteriser, three new score grids (`temp_avg`, `temp_max`, `ice_depth`) on the common 240 m grid, thermal+ice criteria switched to PRP defaults. Five of six criteria now run on real data. Validation rerun *surfaced the methodology finding*: linear-sum MCDA can't model NASA's spatial-coupling constraint regardless of which criteria are summed.
- **v0.7 — spatial-coupling criterion.** ✅ `criteria/coupling.py` scores cells by joint proximity to a PSR and a sunlit ridge via a multiplicative product of two distance falloffs, encoding the AND that linear-sum aggregation cannot. `selene coupling-sweep` tunes the single `coupling_distance_km` knob across 1–20 km. The fix lifted the sensitivity ceiling from 2/9 -> 3/9 region matches and improved closest-distance from 64.8 km -> 47.8 km, but did not break the 25 km threshold. The diagnostic surfaced the project's sharpest finding: **NASA centroids score 0.000 on coupling** — meaning the validation metric (distance to disk centroid) is itself misaligned with NASA's selection logic. NASA's preferred landing footprints within each region are off-centroid by 5–15 km by construction.
- **v1.0.0 — closing chapter: polygon validation, thermal correction, engineering decisions.** ✅ Three changes, all driven by v0.7's diagnostic. (1) `validation/comparison.py` now computes polygon-inside metrics alongside the legacy centroid-distance metrics — `sites_inside_any_region`, `regions_containing_top_site`, `regions_with_top_site_within_disk_radius`, plus signed `distance_to_edge_km` per site/region. The polygon primitive doesn't move the headline (0/20 sites inside any disk; closest edge 32.8 km) — informative on its own, since it shows the rim band the model identifies is geometrically distinct from NASA's centroid disks regardless of which proximity primitive you choose. (2) Thermal default corrected: `target_temp_k` 230 -> 140, `sigma_k` 50 -> 30. The previous values placed the Gaussian peak *outside* the data support (PRP `temp_avg` peaks at 211 K, median 131 K), so every cell scored in the tail. With the correction, our top-20 thermal mean rises 0.325 -> 0.965 and NASA centroids 0.113 -> 0.526; the criterion contributes discriminative signal again, and the sensitivity ceiling broadens from 1/200 to 11/200 weight regimes hitting 3/9 region matches. (3) New "Engineering decisions" section documenting the seven non-obvious choices.
- **v1.1.0 — Earth line-of-sight criterion.** ✅ `criteria/los_to_earth.py` adds the seventh criterion: a per-pixel Earth-visibility fraction derived from a 36-azimuth, log-spaced 50-distance horizon ray-march on the LOLA elevation grid (with curvature correction for $R = 1737.4\,$km), combined with 24-sample libration-cycle sampling of Earth's sub-Earth point on a $\pm6.5°\times\pm7.9°$ ellipse. Score is a linear ramp from `min_visibility = 0.20` (Apollo crew-safety floor) to `target_visibility = 0.50` (sustained-habitat target). Default weight 0.15 — chosen *before* the validation rerun on physics-and-operations grounds, not validation chasing. Effect: closest-disk-edge dropped from 32.8 km → **12.3 km** (Cabeus B, the first NASA region to come within 1 disk-radius of a top site), sensitivity ceiling lifted 3/9 → **4/9** regions matched, and the polygon-inside count became sample-non-zero for the first time (**21/200 weight regimes** now produce ≥1 inside any disk; 0/200 in v1.0.0). Headline polygon-inside under defaults: still **0/20** — the operations-driven defaults narrow the gap but don't close it. The geometric finding from v0.7–v1.0.0 holds: under physics-driven defaults, our model picks the polar rim band where every criterion (now including LOS) is near-saturated, and NASA's centroids are inside the disks 5–15 km off that band; the gap collapses for ~10 % of the weight simplex.
- **v1.2.0 — USGS authoritative polygon validation.** ✅ Replaced the 15 km disk approximations with USGS's officially-published simplified region envelopes (DOI 10.5066/P1MEQ6UK, McClernan 2024). The polygons ship in-repo at [`src/selene_base/validation/data/nasa_regions_polygons_usgs.geojson`](src/selene_base/validation/data/nasa_regions_polygons_usgs.geojson). `selene validate` now prints three result tables (centroid distance, 15 km disk inside/outside, USGS polygon inside/outside) and `validation.json` carries all three metric families. The disk approximations were systematically wrong: most regions are ~400 km² quadrilaterals (vs the 707 km² disk), Mons Mouton Plateau is **4452 km² — 6× the disk area**, one disk centroid for "Slater Plain" sits ~180° away from the USGS polygon's actual location, and "Cabeus B" was misnamed (USGS publishes "Peak Near Cabeus B", centred on the rim peak, not the crater floor). **Default-weights result against USGS polygons: 0/20 inside, 0/9 USGS regions containing a top site, median distance 135.1 km, closest 41.5 km (de Gerlache Rim 2). Sensitivity sweep: 6/200 weight regimes produce ≥1 site inside any USGS polygon (max 2/20), down from 21/200 against the disks** — the geometric separation is *more* pronounced against the right validation reference, not less, because the disks were inflated outward in places where the actual USGS polygons are not. The "validation metric was the bottleneck" hypothesis from v1.0.0 was *partially* correct (the disk approximations were wrong) but the *underlying* geometric separation between the model's rim-band optimum and NASA's authoritative regions persists.
- **v1.3.0 — per-region ranking with NASA HLS hard filters.** ✅ `selene rank-per-region` searches *within* each USGS polygon and applies NASA's published HLS thresholds (slope ≤ 8°, 100 m buffer, illumination ≥ 33 %, DTE visibility ≥ 50 %) as a precondition before ranking by suitability score. New `top_n_sites_per_region` in [`scoring/ranking.py`](src/selene_base/scoring/ranking.py); new `per_region_compliance_analysis` in [`validation/comparison.py`](src/selene_base/validation/comparison.py); new CLI subcommands `rank-per-region` and `validate-per-region`. **Result: 23 sites across 8/9 USGS regions** (default weights, default HLS thresholds), all guaranteed inside their named polygon and HLS-compliant by construction. **Malapert Massif has zero HLS-compliant cells** — a real terrain-driven finding, not a thresholding artefact. **Mons Mouton Plateau is the highest-scoring region** (best score 0.746, 15.07 % HLS-eligible area). The reframing — global ranking → per-region + HLS — matches the methodology of Wueller et al. 2026 (JGR Planets), which catalogued 130 candidate Artemis-III sites with the same outer framing. The 0/20 inside-polygon count through v1.2.0 reflected the global framing; v1.3.0 produces the NASA-aligned catalog v1.0–v1.2 had been pursuing through the wrong question.
- **v1.4.0 — Wueller 2026 comparison framework (framework only).** ✅ Comparison harness shipped: `selene compare-wueller` plus `wueller_comparison.{load_wueller_sites,compare_sites,render_summary}`, 16 synthetic-only tests, the `notebooks/08_wueller_comparison.py` visualisation set, and the headline plot at [docs/img/selene_vs_wueller.png](docs/img/selene_vs_wueller.png). The Wueller 2026 supplementary catalog was unavailable at the time; the bundled CSV was a 5-row synthetic placeholder explicitly flagged through every output channel, and two "real comparison" tests were skipped pending data. **The framework was the v1.4.0 deliverable; the quantitative agreement number unblocked in v1.4.1.**
- **v1.4.1 — quantitative Wueller comparison.** ✅ Replaces v1.4.0's synthetic placeholder with the real 130-site shapefile from the authors' Zenodo deposit ([doi:10.5281/zenodo.17084058](https://doi.org/10.5281/zenodo.17084058), CC-BY 4.0), bundled in-repo at [`src/selene_base/validation/data/wueller_2026/`](src/selene_base/validation/data/wueller_2026/). `compare_sites` now defaults to a USGS-scope filter (drops the 57 Wueller sites whose region is not in NASA's October 2024 down-selected nine), the two formerly-skipped tests now run against the real bundle, and the legacy CSV path is retained behind a `DeprecationWarning`. **Result at the v1.3 default `n_per_region = 3`: 18 / 23 selene sites match within 5 km of an in-scope Wueller site (78 % agreement); median match distance 1.71 km against 73 in-scope Wueller sites across 8 USGS regions.** Six of eight regions agree at 100 % match-within-region; three regions (de Gerlache Rim 2, Mons Mouton Plateau, Slater Plain) have outliers that miss the threshold by ≤ 6 km. Test suite: 350 collected, 346 passed, 4 skipped (LEND-only).
- **v1.4.2 — `n_per_region` default raised to 10 (matches Wueller's per-region site density).** ✅ Default `n_per_region` changed from 3 to 10. The selene methodology converges at any per-region density: selene-to-Wueller match rate is essentially flat across n (78 % at n=3, 78 % at n=5, 80 % at n=10), but Wueller-to-selene match rate climbs from 41 % to 63 % as more selene sites become available to pair with Wueller's 73-site in-scope catalog. **Result at the new default: 56/70 selene sites (80 %) match within 5 km of an in-scope Wueller site; median match distance 1.88 km; 46/73 (63 %) of in-scope Wueller sites match a selene site within 5 km.** Two regions (Mons Mouton, Peak Near Cabeus B) reach 100 % match-within-region at n=10; Nobile Rim 1 and Nobile Rim 2 reach 80 %+. de Gerlache Rim 2 stays at 0/2 (the polygon's HLS-eligible area caps at 2 sites at 2km NMS — its terrain is genuinely far from Wueller's dGR2 cluster, independent of n). Per-region site counts (NMS-capped at 2km separation): HW 10, MM 10, MMP 10, NR1 10, NR2 10, PCB 9, SP 9, dGR2 2 — total 70.
- **v1.5 — 20 m Wueller-class resolution with GPU acceleration.** ✅ Per-region tiled driver runs the HLS filters at 20 m, the same native LOLA resolution Wueller et al. 2026 use. New CLI flags `--tiled-per-region --resolution 20` on `selene preprocess` and `selene rank-per-region`; new GPU path through `derive_horizon_profile` and `compute_earth_visibility_fraction` via CuPy + `cupyx.scipy.ndimage.map_coordinates`; new `pipeline/preprocess_tiled.py` and `pipeline/rank_per_region_tiled.py` covering each USGS polygon's bbox + 100 km horizon buffer with explicit cupy memory-pool drain between tiles. The 240 m global path is unchanged. **Result against Wueller 2026 at 20 m: 56/69 selene sites (81 %) match within 5 km, median matched-pair distance 1.76 km — tighter than v1.4.2's 1.88 km at 240 m by 6.5 %; selene-to-Wueller match rate flat (80 % → 81 %); Wueller-to-selene rate slips slightly (63 % → 60 %) because the finer slope/buffer re-ranks sites within MMP, SP, and CB to a different cell within the same polygon.** Six of eight regions tighten or hold their median matched-pair distance at 20 m; Malapert Massif still has zero HLS-compliant cells (genuine steep terrain, not a 240 m sampling effect); de Gerlache Rim 2 still produces only 2 sites with 0/2 within 5 km of Wueller's dGR2 cluster (genuine terrain divergence, persists at fine resolution). The methodology converges at peer-reviewed resolution; the residual disagreements at v1.4.2 are not 240 m sampling artefacts. v1.5 is opt-in and requires a CUDA GPU plus ~140 GB unified-memory headroom (developed on NVIDIA GB10 / DGX Spark; full 9-region run ~32 min preprocess + ~22 min rank). See "Resolution analysis (v1.5 — 20 m Wueller-class)" above for the per-region table and side-by-side plots.
- **v1.5.1 — Artemis IV mission designation update + NASA Figures-of-Merit framework citation.** ✅ Docs-only patch release. NASA restructured the Artemis program in February 2026: the first crewed lunar landing at the south pole was reassigned from Artemis III to Artemis IV (target early 2028), Artemis III became an Earth-orbit test mission, and Lunar Gateway was canceled in March 2026. v1.5.1 reframes the project's forward-looking copy from "Artemis III" to "Artemis IV (formerly Artemis III)" across the README, the v1.5 catalog report, and user-facing CLI / module docstrings, while leaving citation titles, the USGS dataset name, NASA's October 2024 announcement URL, the Wueller 2026 paper's own framing, and the historical roadmap entries unchanged. Adds [Lawrence 2025](https://ntrs.nasa.gov/citations/20250008952) (NASA NTRS 20250008952) as a reference and a new "Mapping to NASA's Figures of Merit framework" subsection in the Methodology section that maps selene's seven criteria to NASA's published FOM categories. **The methodology, data, validation, and 69-site catalog are identical to v1.5; the headline numbers (56/69 = 81.2 % within 5 km, median 1.76 km) are unchanged.**
- **v1.7 — TOPSIS aggregator as opt-in alternative to weighted_sum.** ✅ Adds [`scoring/aggregate.topsis`](src/selene_base/scoring/aggregate.py) (Hwang & Yoon 1981) as a sibling of `weighted_sum`, with `selene score --method topsis` and `selene rank-per-region --method topsis` exposing it. Same input contract (per-criterion [0, 1] score grids + weights), same renormalise-on-missing-criteria semantics, same NaN propagation. Default stays `weighted_sum` — v1.5's 81.2 % / 1.76 km headline is unchanged. **TOPSIS result on the v1.5 20 m catalog: same 69 sites (HLS filters identical), same top-5 by score, 55 / 69 (79.7 %) matched within 5 km of an in-scope Wueller site (vs 56 at weighted_sum), median matched-pair distance 1.84 km vs 1.76 km globally** — but **three of seven active regions tighten significantly** (Haworth -0.51 km, Mons Mouton Plateau -0.40 km, Nobile Rim 1 -0.43 km) and **max pair distance improves -0.22 km**. Both aggregators agree on which cells are top-tier; they disagree at the 1.5–4 km margin. TOPSIS outputs land under `data/outputs/topsis/` so the two catalogs sit side-by-side. See "Aggregator: weighted_sum vs TOPSIS" in Methodology for the full table.

- **v1.8 — Activate seismic criterion (eight of eight criteria live).** ✅ Adds distance-to-nearest-scarp as the eighth criterion via a logistic of distance to the nearest mapped lunar lobate scarp (midpoint 25 km, steepness 8 km — tuned to Civilini et al. 2023's documented shallow-moonquake-to-scarp clustering distance). Bundles the **Mishra & Kumar (2022)** south-polar scarp catalog ([GRL doi:10.1029/2022GL098505](https://doi.org/10.1029/2022GL098505), Zenodo [doi:10.5281/zenodo.6624114](https://doi.org/10.5281/zenodo.6624114), CC-BY 4.0) — 704 line segments combining Watters et al. (2015) global mapping with 75 new south-polar scarps Mishra & Kumar mapped from independent LROC NAC analysis — at `src/selene_base/criteria/data/scarps_mishra_kumar_2022/`, with the original LROC SOC `POLAR_SCARP_LOCATIONS` shapefile (citing Watters 2015 verbatim) bundled alongside as an attribution anchor at `scarps_watters_2015_polar/`. Default seismic weight 0.10, with the other seven weights scaled by 0.928 so the simplex still sums to 1.00. Methodology fully realised — eight of eight planned criteria contribute to the aggregate score. **Result: 56 / 69 selene sites (81.2 %) match within 5 km of an in-scope Wueller site — identical to v1.5's headline; 46 / 73 (63.0 %) Wueller-side match (up from 44 / 73 in v1.5, +2.7 pp); median matched-pair distance 1.69 km vs 1.76 km in v1.5 (-0.07 km).** Same site count, same top-3 in every region; 13 sites in HW/MMP/N1/N2 reshuffle within ranks 4-10. The criterion is most active where the Mishra/Kumar catalog has dense south-polar coverage and least active for cells far from any mapped scarp (most of the polar grid). See "Seismic criterion (v1.8)" in Methodology for the logistic scoring details and per-region table.

### Where this goes next

**Wueller 2026 deposit — non-point layers.** The Zenodo deposit also contains an 884 MB HLS slope raster (`HLS_LandingAreas_8°_-100m_buffer.tif`), 2 km exploration-area polygons (`ExplorationArea_2km.*`), and an illumination-modelling layer set. None are bundled through v1.8 (point catalog only); folding any of them in (e.g. comparing selene's HLS-eligible mask against Wueller's 8°-buffer raster) is a near-term candidate.

**v1.8 sensitivity rerun.** The 200-sample Latin-hypercube weight sweep in [`Robustness`](#robustness) was last run under the v1.5 seven-criterion vector. Re-running it on the v1.8 eight-criterion simplex (with seismic at 0.10) is a straightforward follow-up that would also refresh the diagnostic-comparison per-criterion table.

### Smaller follow-ups

- Auto-tune `--min-score` in the rank pipeline based on the aggregate's percentile distribution.
- ML-based criterion inputs (planned as a separate project, `selene-vision`).

## References

- Robbins, S. J. (2019). *A new global database of lunar impact craters >1–2 km: 1. Crater locations and sizes, comparisons with published databases, and global analysis.* Journal of Geophysical Research: Planets, 124, 871–892. [doi:10.1029/2018JE005592](https://doi.org/10.1029/2018JE005592)
- Mazarico, E., Neumann, G. A., Smith, D. E., Zuber, M. T., & Torrence, M. H. (2011). *Illumination conditions of the lunar polar regions using LOLA topography.* Icarus, 211(2), 1066–1081. [doi:10.1016/j.icarus.2010.10.030](https://doi.org/10.1016/j.icarus.2010.10.030)
- Smith, D. E., et al. (2010). *The Lunar Orbiter Laser Altimeter investigation on the Lunar Reconnaissance Orbiter mission.* Space Science Reviews, 150(1–4), 209–241. [doi:10.1007/s11214-009-9512-y](https://doi.org/10.1007/s11214-009-9512-y)
- Barker, M. K., Mazarico, E., Neumann, G. A., Smith, D. E., Zuber, M. T., & Head, J. W. (2021). *Improved LOLA elevation maps for south pole landing sites: Error estimates and their impact on illumination conditions.* Planetary and Space Science, 203, 105119. [doi:10.1016/j.pss.2020.105119](https://doi.org/10.1016/j.pss.2020.105119). (Source of the LOLA south polar 20 m DEM `ldem_80s_20m` used as the elevation input for v1.5's per-region tiled high-resolution analysis.)
- Paige, D. A., et al. (2010). *The Lunar Reconnaissance Orbiter Diviner Lunar Radiometer Experiment.* Space Science Reviews, 150(1–4), 125–160. [doi:10.1007/s11214-009-9529-2](https://doi.org/10.1007/s11214-009-9529-2)
- Williams, J.-P., et al. (2017). *The global surface temperatures of the Moon as measured by the Diviner Lunar Radiometer Experiment.* Icarus, 283, 300–325. (PRP modeled-ice-stability methodology.)
- Diviner Polar Resource Product (PRP), south pole, version 1.0. PDS Geosciences Node Diviner derived bundle: [`dlre_prp_south.tab`](https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_south.tab).
- Mitrofanov, I. G., et al. (2010). *Hydrogen mapping of the lunar south pole using the LRO Neutron Detector Experiment LEND.* Science, 330(6003), 483–486.
- Watters, T. R., Robinson, M. S., Banks, M. E., Tran, T., & Denevi, B. W. (2015). *Global thrust faulting on the Moon and the influence of tidal stresses.* Geology, 43(10), 851–854. [doi:10.1130/G37120.1](https://doi.org/10.1130/G37120.1)
- Civilini, F., Weber, R. C., Jiang, Z., Phillips, D., & Pan, W. (2023). *Constraints on the seismic hazard of young thrust faults on the Moon from re-located shallow moonquakes.* (Used as physical motivation + the order-10–50 km clustering-distance scale for v1.8's seismic-criterion logistic.)
- Mishra, A., & Kumar, P. S. (2022). *Spatial and Temporal Distribution of Lobate Scarps in the Lunar South Polar Region: Evidence for Latitudinal Variation of Scarp Geometry, Kinematics and Formation Ages, Neo-Tectonic Activity and Sources of Potential Seismic Risks at the Artemis Candidate Landing Regions.* Geophysical Research Letters, 49, e2022GL098505. [doi:10.1029/2022GL098505](https://doi.org/10.1029/2022GL098505). Data deposit: [doi:10.5281/zenodo.6624114](https://doi.org/10.5281/zenodo.6624114) (CC-BY 4.0). The 704 south-polar scarp line segments ship in-repo at [`src/selene_base/criteria/data/scarps_mishra_kumar_2022/`](src/selene_base/criteria/data/scarps_mishra_kumar_2022/) and are the primary data source for v1.8's seismic criterion.
- Lobate-scarp catalog also distributed by the LROC Science Operations Center (NASA PDS LROLRC_2001 EXTRAS) as `POLAR_SCARP_LOCATIONS` — citing Watters et al. 2015 verbatim. The 48-named-scarp polar shapefile ships in-repo at [`src/selene_base/criteria/data/scarps_watters_2015_polar/`](src/selene_base/criteria/data/scarps_watters_2015_polar/) as an attribution anchor.
- NASA (October 2024). *Artemis III candidate landing regions.* [https://www.nasa.gov/feature/artemis-iii](https://www.nasa.gov/feature/artemis-iii)
- McClernan, M.T. (2024). *Down Selected Artemis III Candidate Landing Site Navigational Grids.* U.S. Geological Survey data release. [https://doi.org/10.5066/P1MEQ6UK](https://doi.org/10.5066/P1MEQ6UK). The simplified region envelopes from this release ship in-repo at [`src/selene_base/validation/data/nasa_regions_polygons_usgs.geojson`](src/selene_base/validation/data/nasa_regions_polygons_usgs.geojson) and are the authoritative validation reference from v1.2.0 onwards.
- Gracy, S., & Lee, P. (2024). *Update on the Artemis III Reference Mission and Candidate Landing Region Selection.* 55th Lunar and Planetary Science Conference, Abstract #1695. (Source for the four published HLS hard-constraint thresholds used by v1.3.0's `rank-per-region` — slope, slope-buffer, illumination, DTE visibility.)
- Lawrence, S. (2025). *Artemis IV Landing Site Process Overview.* NASA Lunar Site Selection Workshop, NTRS document 20250008952. [https://ntrs.nasa.gov/citations/20250008952](https://ntrs.nasa.gov/citations/20250008952). (Documents NASA's 5-step iterative site-selection process and the multi-directorate Figures of Merit (FOM) framework. Specific FOM weight values are not publicly disclosed; selene-base's "Mapping to NASA's Figures of Merit framework" subsection in Methodology positions the project's seven criteria against this published conceptual framework.)
- Wueller, F., et al. (2026). *Assessing Potential Landing Sites With Favorable Illumination and Accessible, Potentially Volatile-Rich Permanently Shadowed Regions Within Artemis Candidate Landing Regions.* Journal of Geophysical Research: Planets, 131. [doi:10.1029/2025JE009434](https://doi.org/10.1029/2025JE009434). (Peer-reviewed parallel: 130 candidate Artemis III landing sites identified by within-region selection with NASA HLS hard filters. v1.3.0's `rank-per-region` mirrors the outer framing; v1.4.1+ ships the quantitative comparison against this catalog — see README §"Quantitative comparison against Wueller et al. 2026".)
- Wueller, L., Berger, L. M., Christopher, H., Sugimoto, K., Thaker, A., Carton, L., Jo, W., Lee, S., Pedrelli, R., Sanchez, P., & Kring, D. (2025). *Complementary Data for Wueller et al. (2026): Assessing Potential Landing Sites with Favorable Illumination and Accessible, Potentially Volatile-Rich Permanently Shadowed Regions within Artemis Candidate Landing Regions* (Version v1) [Data set]. Zenodo. [doi:10.5281/zenodo.17084058](https://doi.org/10.5281/zenodo.17084058). (CC-BY 4.0 data deposit accompanying the JGR paper. The 130-site point catalog ships in-repo at [`src/selene_base/validation/data/wueller_2026/`](src/selene_base/validation/data/wueller_2026/) and is the comparison reference from v1.4.1 onwards.)
- NASA (2019). *Human Landing System Requirements Document.* (Underlying source for the HLS slope, slope-buffer, illumination, and DTE-visibility thresholds.)

## License

MIT — see [LICENSE](LICENSE).
