"""End-to-end loader for the Diviner Polar Resource Product.

Parses the PDS4 ``dlre_prp_south.tab`` table once via
:mod:`selene_base.data.pds4_table`, rasterises each derived layer onto
the project's south-polar grid via
:mod:`selene_base.data.triangle_to_grid`, and caches the result as a
GeoTIFF so subsequent calls skip the slow parse step.

The cached layers feed the thermal and ice criteria:

- ``temp_avg`` — annual mean surface temperature at 2 cm depth (K).
- ``temp_max`` — annual maximum surface temperature (K).
- ``ice_depth`` — modeled depth (m) at which water ice is stable
  against sublimation; 0 means stable at the surface, NaN means
  effectively no ice (sentinel ``-999`` in the source file).

Filled in week 6.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import rioxarray  # noqa: F401  (registers .rio accessor)
import typer
import xarray as xr

from selene_base.data.pds4_table import parse_pds4_label, read_pds4_table
from selene_base.data.triangle_to_grid import triangles_to_raster

DEFAULT_RAW_DIR = Path("data/raw/diviner")
DEFAULT_CACHE_DIR = Path("data/processed")

PRP_TAB_FILENAME = "dlre_prp_south.tab"
PRP_XML_FILENAME = "dlre_prp_south.xml"

# Field name in the PDS4 label  ->  rasterisation method.
LAYER_METHODS: dict[str, str] = {
    "temp_avg": "linear",  # smooth field
    "temp_max": "linear",
    # ice_depth is effectively discontinuous; linear smearing across the boundary is unphysical.
    "ice_depth": "nearest",
}


def _cache_path(cache_dir: Path, name: str) -> Path:
    return cache_dir / f"diviner_{name}_southpole_240m.tif"


def load_diviner_prp(
    raw_dir: Path = DEFAULT_RAW_DIR,
    target_grid: xr.DataArray | None = None,
    *,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    overwrite: bool = False,
    echo: Callable[[str], None] = typer.echo,
) -> dict[str, xr.DataArray]:
    """Load the PRP and rasterise its three derived layers.

    Args:
        raw_dir: Directory containing ``dlre_prp_south.{tab,xml}``.
        target_grid: DataArray on the projected analysis grid; required
            unless every cache file is already present (in which case
            we read from cache without ever needing the grid).
        cache_dir: Where the rasterised layers are cached as GeoTIFFs.
        overwrite: Re-parse and re-rasterise even when caches exist.
        echo: Logging sink.

    Returns:
        Mapping from layer name (``temp_avg``, ``temp_max``,
        ``ice_depth``) to an :class:`xr.DataArray` on the target grid.

    Raises:
        FileNotFoundError: If the PRP source files are missing AND any
            of the cache files are missing.
        ValueError: If ``target_grid`` is missing while the cache is
            incomplete.
    """
    raw_dir = Path(raw_dir)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_paths = {name: _cache_path(cache_dir, name) for name in LAYER_METHODS}
    have_all_caches = all(p.exists() for p in cache_paths.values())

    if have_all_caches and not overwrite:
        echo(f"[skip] all PRP caches present in {cache_dir}")
        return {
            name: rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)
            for name, path in cache_paths.items()
        }

    tab_path = raw_dir / PRP_TAB_FILENAME
    xml_path = raw_dir / PRP_XML_FILENAME
    if not tab_path.exists() or not xml_path.exists():
        raise FileNotFoundError(
            f"PRP source missing under {raw_dir}; run `selene download diviner` first"
        )
    if target_grid is None:
        raise ValueError(
            "target_grid is required when the cache is incomplete; "
            "pass an aligned DataArray (e.g. the cached LOLA COG)"
        )

    echo(f"[parse] {tab_path.name}")
    spec = parse_pds4_label(xml_path)
    fields = ["tri_clat", "tri_clon", "temp_avg", "temp_max", "ice_depth"]
    df = read_pds4_table(tab_path, spec, fields=fields)
    echo(f"[parse] {len(df):,} triangle records loaded")

    layers: dict[str, xr.DataArray] = {}
    for layer, method in LAYER_METHODS.items():
        path = cache_paths[layer]
        if path.exists() and not overwrite:
            echo(f"[skip] {layer}: {path.name} already cached")
            layers[layer] = rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)
            continue
        echo(f"[rasterise] {layer} ({method})")
        raster = triangles_to_raster(
            points_lat=df["tri_clat"].to_numpy(),
            points_lon=df["tri_clon"].to_numpy(),
            values=df[layer].to_numpy(),
            target_grid=target_grid,
            method=method,
        )
        from selene_base.data.reproject import cache_processed

        out_path = cache_processed(raster, f"diviner_{layer}", cache_dir, overwrite=overwrite)
        echo(f"[done] {layer} -> {out_path}")
        layers[layer] = rioxarray.open_rasterio(out_path, masked=True).squeeze("band", drop=True)

    return layers
