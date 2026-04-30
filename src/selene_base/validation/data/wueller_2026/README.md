# Wueller et al. (2026) Landing Site Catalog

130 candidate Artemis III landing sites from Wueller et al. (2026)'s
Zenodo data deposit, used by selene-base for v1.4.1 comparison
validation.

## Source

Wueller, L., Berger, L. M., Christopher, H., Sugimoto, K., Thaker, A.,
Carton, L., Jo, W., Lee, S., Pedrelli, R., Sanchez, P., & Kring, D.
(2025). *Complementary Data for Wueller et al. (2026): Assessing
Potential Landing Sites with Favorable Illumination and Accessible,
Potentially Volatile-Rich Permanently Shadowed Regions within Artemis
Candidate Landing Regions* (Version v1) [Data set]. Zenodo.
<https://doi.org/10.5281/zenodo.17084058>

## Original publication

Wueller, F., et al. (2026). *Assessing Potential Landing Sites With
Favorable Illumination and Accessible, Potentially Volatile-Rich
Permanently Shadowed Regions Within Artemis Candidate Landing Regions.*
Journal of Geophysical Research: Planets, 131. doi:10.1029/2025JE009434

## License

Creative Commons Attribution 4.0 International (CC-BY 4.0).
<https://creativecommons.org/licenses/by/4.0/>

Bundled in this repository under the terms of CC-BY 4.0. Cite the
original paper and the Zenodo deposit in any derived work.

## Files

- `LandingSites.shp`, `.shx`, `.dbf`, `.prj`, `.cpg` — shapefile bundle
  of 130 candidate landing sites with attributes from the original
  deposit. Coordinate reference system is `Moon2000_spole` (lunar
  south-polar stereographic, R = 1 737 400 m, false easting/northing
  zero, latitude of origin -90°, standard parallel -70°). Geometry is
  Point in projected metres; `Latitude` and `Longitude` planetocentric
  degrees ship as DBF columns.

## Schema

The DBF carries the upstream column names verbatim. The
selene-base loader (`load_wueller_sites`) renames a small subset to
match the project schema and adds two derived columns:

| upstream column | selene-base column | notes                                 |
| --------------- | ------------------ | ------------------------------------- |
| `Name`          | `wueller_site_id`  | site label, e.g. `MMO01`, `NR101`     |
| `Landing_Re`    | `wueller_region`   | 3-letter region code, e.g. `MMO`      |
| `Latitude`      | `lat`              | planetocentric degrees                |
| `Longitude`     | `lon`              | planetocentric degrees                |
| —               | `region`           | USGS canonical name, or code if none  |
| —               | `in_usgs_scope`    | True if region is one of NASA's 9     |

All other DBF columns (`SunDays25`–`SunDays32`, `SunEarth26`–`SunEarth32`,
`PSR_AREA`, `coldest_te`, `oldest_PSR`, `StereoX`, `StereoY`, `ID`) are
preserved verbatim and may be consumed by future comparison metrics.

## Region distribution

| Code  | Wueller region              | n  | USGS scope                      |
| ----- | --------------------------- | -- | ------------------------------- |
| AR    | Amundsen Rim                | 10 | out (not in USGS 9)             |
| CR    | Connecting Ridge            |  9 | out (removed Oct 2024)          |
| CRE   | Connecting Ridge Extension  |  6 | out (removed Oct 2024)          |
| FRA   | Faustini Rim A              |  8 | out (removed Oct 2024)          |
| HW    | Haworth                     | 11 | **in** — Haworth                |
| MMA   | Mons Malapert               |  5 | out (removed Oct 2024)          |
| MMO   | Mons Mouton                 | 10 | **in** — Mons Mouton            |
| MMP   | Mons Mouton Plateau         | 11 | **in** — Mons Mouton Plateau    |
| NR1   | Nobile Rim 1                |  9 | **in** — Nobile Rim 1           |
| NR2   | Nobile Rim 2                |  9 | **in** — Nobile Rim 2           |
| PCB   | Peak Near Cabeus B          |  5 | **in** — Peak Near Cabeus B     |
| PNS   | Peak Near Shackleton        |  6 | out (removed Oct 2024)          |
| SP    | Slater Plain                | 11 | **in** — Slater Plain           |
| dGKM  | de Gerlache-Kocher Massif   |  7 | out (removed Oct 2024)          |
| dGR1  | de Gerlache Rim 1           |  6 | out (removed Oct 2024)          |
| dGR2  | de Gerlache Rim 2           |  7 | **in** — de Gerlache Rim 2      |
| Total |                             |130 | 73 in scope, 57 out             |

USGS scope = NASA's October 2024 down-selected nine regions, as shipped
in `nasa_regions_polygons_usgs.geojson`. Wueller pre-dates that
down-selection; the 57 out-of-scope sites cluster in regions that were
on NASA's earlier 13-region list but were dropped at down-selection
time. selene-base's `compare-wueller` filters to in-scope regions by
default; pass `--no-filter-to-usgs-scope` (or `filter_to_usgs_scope=False`
to `compare_sites`) to compare against all 130 sites.
