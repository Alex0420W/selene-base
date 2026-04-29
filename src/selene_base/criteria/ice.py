"""Ice criterion — rewards proximity to inferred water-ice deposits.

LEND epithermal-neutron suppression maps trace hydrogen abundance, which
under polar conditions is widely interpreted as buried water ice.
The score is monotonically increasing in the inferred wt% H2O proxy.

Filled in week 3.
"""

from __future__ import annotations

import xarray as xr


def compute(grid: xr.DataArray, **kwargs: object) -> xr.DataArray:
    """Score every cell on inferred water-ice resource potential.

    Args:
        grid: DataArray on the common south-polar grid holding a
            hydrogen-abundance proxy (typically wt% water-equivalent).
        **kwargs: Tuning knobs. Recognised keys:

            * ``hi_wt_pct`` (float, default 0.5) — value mapped to score 1.
            * ``lo_wt_pct`` (float, default 0.0) — value mapped to score 0.

    Returns:
        DataArray of [0, 1] scores aligned with ``grid``; NaN where
        hydrogen data is missing.

    Raises:
        NotImplementedError: Implementation is filled in week 3.
    """
    raise NotImplementedError("filled in week 3")
