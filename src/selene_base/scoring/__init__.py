"""Scoring subpackage — normalise, aggregate, rank.

Three concerns kept separate so each is independently testable:

- :mod:`selene_base.scoring.normalize` — pure-numpy [0, 1] mappings
  (real implementations; week 1).
- :mod:`selene_base.scoring.aggregate` — weighted combination of per-
  criterion score grids (real implementation; week 1).
- :mod:`selene_base.scoring.ranking` — non-maximum suppression to extract
  geographically-distinct top sites (stub; week 3).
"""
