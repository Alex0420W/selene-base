"""Validation against NASA's Artemis III candidate landing regions.

Two modules:

- :mod:`selene_base.validation.nasa_regions` — the public-information
  list of NASA's nine announced south-polar candidate regions and a
  helper to materialise them as disk-approximation polygons in any
  target CRS.
- :mod:`selene_base.validation.comparison` — proximity metrics between
  selene-base's ranked sites and the NASA candidates.
"""
