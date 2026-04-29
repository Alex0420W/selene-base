"""Per-criterion suitability scorers.

Each submodule exports ``compute(grid, **kwargs) -> xr.DataArray`` and
returns a [0, 1] score map aligned with the input grid. Scores follow the
convention "1 is good, 0 is unusable" so that aggregation in
:mod:`selene_base.scoring.aggregate` is a straight weighted sum.

Implementations land week 2 (slope first, alongside reproject) and week 3
(remaining criteria).
"""
