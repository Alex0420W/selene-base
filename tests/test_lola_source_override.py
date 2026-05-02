"""Tests for the explicit ``--lola-source`` override (v1.8.3).

The auto-detect path of :func:`resolve_lola_source` is hardcoded to
the PDS naming convention (``ldem_80s_{resolution_m}m.lbl``) and would
silently pick a coarser file when a non-PDS source (PGDA mosaic,
NAC stereo DEM, custom DEM) is on disk. The ``override=`` parameter
and the equivalent ``--lola-source`` CLI flag bypass the priority
list with caller-supplied path validation.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr
from typer.testing import CliRunner

from selene_base.cli import app
from selene_base.data.load import LUNAR_SOUTH_POLAR_CRS
from selene_base.pipeline import preprocess_tiled

POLAR_CRS = str(LUNAR_SOUTH_POLAR_CRS)

# CI runners preserve ANSI colour codes and Rich's panel/line-wrapping that
# local terminals strip; normalise both before substring-matching CLI output.
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def _normalise_cli_output(text: str) -> str:
    """Strip ANSI escapes and collapse all whitespace runs to single spaces."""
    return " ".join(_ANSI_RE.sub("", text).split())


def _tiny_synthetic_tif(path: Path, resolution_m: float = 240.0) -> Path:
    """Write a 32×32 synthetic GeoTIFF on the south-polar grid for resolver tests.

    Resolver tests only need a file whose suffix and CRS are valid; the
    pixel content is irrelevant. Kept small (~4 KB) so the test stays
    fast even on cold-cache fixtures.
    """
    size = 32
    half = (size * resolution_m) / 2.0
    xs = np.linspace(-half + resolution_m / 2, half - resolution_m / 2, size)
    ys = np.linspace(half - resolution_m / 2, -half + resolution_m / 2, size)
    z = np.zeros((size, size), dtype=np.float32)
    da = xr.DataArray(z, dims=("y", "x"), coords={"y": ys, "x": xs}, name="elevation_m")
    da = da.rio.write_crs(POLAR_CRS, inplace=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    da.rio.to_raster(path, driver="GTiff", compress="DEFLATE")
    return path


def test_resolve_lola_source_override_uses_explicit_path(tmp_path: Path) -> None:
    src = _tiny_synthetic_tif(tmp_path / "custom" / "pgda_like.tif")

    resolved = preprocess_tiled.resolve_lola_source(raw_dir=tmp_path / "raw", override=src)

    assert resolved == src


def test_resolve_lola_source_override_bypasses_priority_list(tmp_path: Path) -> None:
    """Override wins even when an auto-detectable PDS file is on disk."""
    raw = tmp_path / "raw"
    pds_like = raw / "lola" / "ldem_80s_20m.lbl"
    pds_like.parent.mkdir(parents=True, exist_ok=True)
    pds_like.write_text("")  # would be picked by the priority list

    custom = _tiny_synthetic_tif(tmp_path / "custom" / "pgda_like.tif")

    resolved = preprocess_tiled.resolve_lola_source(raw_dir=raw, override=custom)

    assert resolved == custom


def test_resolve_lola_source_override_missing_path_raises(tmp_path: Path) -> None:
    bogus = tmp_path / "does_not_exist.tif"
    with pytest.raises(FileNotFoundError, match="--lola-source path does not exist"):
        preprocess_tiled.resolve_lola_source(raw_dir=tmp_path / "raw", override=bogus)


def test_resolve_lola_source_override_rejects_directory(tmp_path: Path) -> None:
    a_dir = tmp_path / "lola_dir"
    a_dir.mkdir()
    with pytest.raises(FileNotFoundError, match="not a file"):
        preprocess_tiled.resolve_lola_source(raw_dir=tmp_path / "raw", override=a_dir)


def test_resolve_lola_source_override_rejects_unknown_suffix(tmp_path: Path) -> None:
    bad = tmp_path / "elevation.txt"
    bad.write_text("not a raster")
    with pytest.raises(ValueError, match="suffix '.txt' not recognised"):
        preprocess_tiled.resolve_lola_source(raw_dir=tmp_path / "raw", override=bad)


def test_resolve_lola_source_no_override_preserves_priority(tmp_path: Path) -> None:
    """Default (no override) still uses the auto-detect priority list."""
    raw = tmp_path / "raw"
    pds_20m = raw / "lola" / "ldem_80s_20m.lbl"
    pds_20m.parent.mkdir(parents=True, exist_ok=True)
    pds_20m.write_text("")
    pds_80m = raw / "lola" / "ldem_80s_80m.lbl"
    pds_80m.write_text("")

    resolved = preprocess_tiled.resolve_lola_source(raw_dir=raw)

    assert resolved == pds_20m


def test_resolve_lola_source_no_override_no_files_raises(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    (raw / "lola").mkdir(parents=True, exist_ok=True)
    with pytest.raises(FileNotFoundError, match="No LOLA source"):
        preprocess_tiled.resolve_lola_source(raw_dir=raw)


def test_run_tiled_per_region_routes_source_path_through_resolver(
    tmp_path: Path,
) -> None:
    """Passing source_path= to run_tiled_per_region triggers override validation."""
    bogus = tmp_path / "missing.tif"
    with pytest.raises(FileNotFoundError, match="--lola-source path does not exist"):
        preprocess_tiled.run_tiled_per_region(
            resolution_m=240.0,
            region_codes=["SP"],
            processed_dir=tmp_path / "processed",
            source_path=bogus,
            use_gpu=False,
            echo=lambda _msg: None,
        )


def test_cli_lola_source_flag_outside_tiled_mode_errors(tmp_path: Path) -> None:
    """--lola-source without --tiled-per-region is rejected with a clear message."""
    src = _tiny_synthetic_tif(tmp_path / "custom" / "src.tif")

    result = CliRunner().invoke(app, ["preprocess", "--lola-source", str(src)])

    assert result.exit_code != 0
    # Rich renders the error inside a bordered panel that wraps long lines and
    # injects ANSI styling on CI terminals; check distinctive tokens that
    # cannot straddle a wrap or escape-sequence boundary after normalisation.
    flat = _normalise_cli_output(result.output)
    assert "Invalid value" in flat
    assert "lola-source" in flat
    assert "tiled-per-region" in flat


def test_cli_lola_source_flag_appears_in_help() -> None:
    """The CLI surface advertises the flag so users can discover it."""
    result = CliRunner().invoke(app, ["preprocess", "--help"])
    assert result.exit_code == 0, result.output
    # Rich's help renderer may colourise option names and wrap the help body,
    # so look for the bare flag stem after stripping ANSI + whitespace.
    flat = _normalise_cli_output(result.output)
    assert "lola-source" in flat
