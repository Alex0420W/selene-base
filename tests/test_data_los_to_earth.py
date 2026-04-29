"""Smoke test for the LOS-to-Earth visibility-fraction COG (week 9).

Loads the cached visibility raster (after ``selene preprocess``) and
checks that the values are in ``[0, 1]`` and that they vary
non-trivially across the polar grid. Skipped automatically when the
cache file is not present so CI stays green without ~900 MB of LOLA
data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr

PROCESSED = Path("data/processed")
VISIBILITY_COG = PROCESSED / "los_visibility_fraction_southpole_240m.tif"
HORIZON_NPZ = PROCESSED / "lola_horizon_profile_southpole_240m.npz"

pytestmark = pytest.mark.skipif(
    not VISIBILITY_COG.exists(),
    reason=f"LOS visibility cache not present: {VISIBILITY_COG.as_posix()}",
)


def _load_visibility() -> xr.DataArray:
    return rioxarray.open_rasterio(VISIBILITY_COG, masked=True).squeeze("band", drop=True)


def test_visibility_fraction_in_unit_range() -> None:
    da = _load_visibility()
    arr = da.to_numpy()
    finite = arr[np.isfinite(arr)]
    assert finite.size > 0
    assert finite.min() >= 0.0
    assert finite.max() <= 1.0


def test_visibility_distribution_non_degenerate() -> None:
    """The visibility fraction should not be all-zero or all-one — at
    the south pole, libration cycles Earth above/below horizon so we
    expect a non-trivial spread.
    """
    da = _load_visibility()
    arr = da.to_numpy()
    finite = arr[np.isfinite(arr)]
    p10 = float(np.percentile(finite, 10))
    p50 = float(np.percentile(finite, 50))
    p90 = float(np.percentile(finite, 90))
    # The grid spans the polar cap; at least *some* cells must score
    # near zero (deep crater floors blocked by rim) and *some* near one
    # (high points with clear horizons). The medians falls somewhere in
    # the middle.
    assert p10 < p50 < p90
    assert p10 < 0.5  # somebody is in the lower half
    assert p90 > 0.3  # somebody is well above zero


def test_horizon_profile_cache_exists_and_3d() -> None:
    """If preprocess ran the LOS step, the .npz cache should be present too."""
    if not HORIZON_NPZ.exists():
        pytest.skip(f"horizon profile cache not present: {HORIZON_NPZ.as_posix()}")
    with np.load(HORIZON_NPZ) as data:
        assert "horizon_profile_deg" in data.files
        assert "azimuth_deg" in data.files
        arr = data["horizon_profile_deg"]
        assert arr.ndim == 3
        finite = arr[np.isfinite(arr)]
        # Horizon angles in degrees: realistically -10° (curvature dip
        # on flat) to +30° (looking up at a near rim). The -90° sentinel
        # value persists for pixels near the grid edge where every ray
        # at every distance leaves the grid (NaN), so the floor of the
        # distribution is exactly -90; the *median* is the meaningful
        # check that horizon values were actually computed.
        assert finite.min() >= -90.0
        assert finite.max() < 90.0
        median = float(np.median(finite))
        assert -10.0 < median < 30.0
