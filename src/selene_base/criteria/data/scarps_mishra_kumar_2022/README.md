# Lunar lobate scarp catalog — Mishra & Kumar (2022)

Primary data source for selene-base's seismic criterion (v1.8+).

## Source

Mishra, A., & Kumar, P. S. (2022). *Spatial and Temporal Distribution of
Lobate Scarps in the Lunar South Polar Region: Evidence for Latitudinal
Variation of Scarp Geometry, Kinematics and Formation Ages,
Neo-Tectonic Activity and Sources of Potential Seismic Risks at the
Artemis Candidate Landing Regions.* Geophysical Research Letters, 49,
e2022GL098505. <https://doi.org/10.1029/2022GL098505>

Data deposit: <https://doi.org/10.5281/zenodo.6624114>

License: **Creative Commons Attribution 4.0 International (CC-BY 4.0)**.

## What's bundled

The bundled shapefile family is the *main* "All segments merged"
product from the Zenodo deposit at
`Datasets_GRL_2022GL098505/Data_Shape_files/All Scarp segments/`:

- `main_segments.shp` (~1.3 MB) — 704 LineString / MultiLineString
  features representing south-polar lunar lobate scarp segments. Builds
  on Watters et al. (2015) global mapping plus 75 new scarps Mishra &
  Kumar mapped from independent LROC NAC analysis.
- `main_segments.shx`, `main_segments.dbf`, `main_segments.prj` —
  shapefile sidecars.

Renamed from the original `All segments merged.{shp,shx,dbf,prj}` to
`main_segments.*` for snake_case consistency with the rest of the
repo. No data transformation; bytes are unchanged.

## Coordinate reference system

Polar Stereographic Moon, sphere R = 1737400 m, latitude of origin
−90°, central meridian 0° — matches selene-base's analysis CRS exactly,
no reprojection required at consumption time.

## What's NOT bundled

- `Data_Shape_files/All Ages/All dated Segments.{shp,shx,dbf,prj}`
  (145 dated points) — out of scope for v1.8's distance-only scoring;
  potential v1.9 candidate to weight scarps by formation age.
- `Data_Shape_files/New Lobate scarps/All new scarps merged.{shp,...}`
  (75 new scarps) — already included in the 704-segment merged file
  bundled above.
- `Data_Age_SCC_files/PLT Files/*` (446 small per-scarp plot files)
  and `Data_Figures/*.xlsx` — figure source data, not catalog input.

## Selene-base usage

`selene_base.criteria.seismic` reads `main_segments.shp` via
`geopandas.read_file`, densifies the polylines to ~1 km point spacing,
builds a `cKDTree`, and produces a per-pixel distance-to-nearest-scarp
score in `[0, 1]`. See the seismic-criterion module docstring for the
sigmoid mapping and the Civilini et al. (2023) physical justification.

## Attribution required

Per CC-BY 4.0, any downstream use must credit Mishra & Kumar (2022) +
the Zenodo DOI. The selene-base README's References section carries
this citation; no further action is required for users running the
pipeline locally.
