# Lunar polar lobate scarp locations — Watters et al. (2015), via LROC SOC

**Attribution anchor / sanity reference** for selene-base's seismic
criterion. Not the primary data source — see the sibling
`scarps_mishra_kumar_2022/` bundle for the 704-segment line catalog
the criterion actually consumes.

## Source

Watters, T. R., Robinson, M. S., Collins, G. C., Banks, M. E., Daud,
K., Williams, N. R., & Selvans, M. M. (2015). *Global thrust faulting
on the Moon and the influence of tidal stresses.* Geology, 43(10),
851–854. <https://doi.org/10.1130/G37120.1>

Hosted by: LROC Science Operations Center, NASA Planetary Data System
(LROLRC_2001 EXTRAS). Direct URL:
<https://pds.lroc.im-ldi.com/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/EXTRAS/SHAPEFILE/POLAR_SCARP_LOCATIONS/POLAR_SCARP_LOCATIONS_180.ZIP>

License: **NASA public data; attribution required.** No registration
required for download.

## What's bundled

- `polar_scarp_locations.shp` (~1.5 KB) — 48 Point features (named
  scarps with LROC NAC image references). 12 in the south polar
  region; the rest are north polar.
- `polar_scarp_locations.shx`, `polar_scarp_locations.dbf`,
  `polar_scarp_locations.prj` — shapefile sidecars.

Renamed from the original `POLAR_SCARP_LOCATIONS_180.{SHP,SHX,DBF,PRJ}`
to lowercase + snake_case. No data transformation.

## Why this is bundled if it isn't the primary data source

Two reasons:

1. **Verbatim attribution to Watters et al. 2015.** The LROC SOC
   shapefile cites Watters 2015 directly as its data citation. Mishra
   & Kumar 2022 (the primary source consumed by the criterion) builds
   on Watters' work and adds 75 new scarps; bundling the original
   alongside makes the attribution chain auditable rather than
   implicit.
2. **Sanity reference.** With only 48 named points, this catalog is
   too sparse to drive the per-pixel distance-to-nearest-scarp
   criterion (only 3 features south of -80°, vs. 704 line segments in
   the Mishra & Kumar bundle). It is, however, useful as a quick
   cross-check that the seismic criterion's distance metric
   behaves sensibly at the named scarps.

## Coordinate reference system

IAU Moon 2000 geographic (lon/lat), sphere R = 1737400 m. Reproject
to selene-base's polar stereographic before any use; the
seismic-criterion module handles this automatically when the file is
loaded.

## Selene-base usage

Not consumed by `selene score` directly. The seismic criterion module
exposes a helper that loads this file for diagnostic comparison
against the Mishra & Kumar primary catalog.
