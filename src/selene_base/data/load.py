"""Disk → typed in-memory loaders for raw products.

Wraps rioxarray / rasterio / geopandas / pandas behind small loader
functions so the rest of the pipeline never sees raw file handles.
Loaders preserve native CRS metadata; coordinate transforms happen in
:mod:`selene_base.data.reproject` (still a stub for week 2).

Lunar geographic CRS (R = 1737400 m, lon/lat) is shared across every
vector product. Raster CRSs are inferred per file from the PDS3 LBL or
GeoTIFF tags.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
import rioxarray  # noqa: F401  (registers the .rio accessor on xarray DataArrays)
import xarray as xr
from pyproj import CRS

LUNAR_GEOGRAPHIC_CRS = CRS.from_proj4("+proj=longlat +R=1737400 +no_defs +type=crs")
LUNAR_SOUTH_POLAR_CRS = CRS.from_proj4(
    "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +R=1737400 +no_defs +type=crs"
)

DEFAULT_RAW_DIR = Path("data/raw")


# ----------------------------------------------------------------------------
# Generic raster loader (LOLA, illumination, Diviner, LEND all rasters)
# ----------------------------------------------------------------------------
def load_raster(path: Path) -> xr.DataArray:
    """Load a single-band raster as a 2-D xarray DataArray.

    Wraps :func:`rioxarray.open_rasterio` and squeezes the trivial band
    dimension. The result has dims ``("y", "x")`` and CRS preserved on
    the ``rio`` accessor. The function does NOT reproject; that's the
    job of :func:`selene_base.data.reproject.reproject_to_grid`.

    Args:
        path: Path to a rasterio-readable file (GeoTIFF, PDS3 IMG with
            sibling LBL, etc.).

    Returns:
        DataArray with dims ``("y", "x")`` in the file's native CRS.
    """
    da = rioxarray.open_rasterio(path, masked=True)
    if isinstance(da, list):  # multi-subdataset case (rare for our products)
        da = da[0]
    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)
    return da


def load_lola_ldem(path: Path = DEFAULT_RAW_DIR / "lola" / "ldem_80s_80m.img") -> xr.DataArray:
    """Load the LOLA LDEM 80°S 80 m DEM as elevation in metres.

    The PDS3 ``.img`` is a 16-bit signed integer raster scaled by 0.5 m.
    rioxarray reads the header from the sibling ``.lbl``; we apply the
    metre scaling so downstream code can assume SI units.

    Args:
        path: Path to ``ldem_80s_80m.img``. The ``.lbl`` must be
            alongside.

    Returns:
        DataArray of elevation in metres on the file's native polar
        stereographic grid.
    """
    da = load_raster(path)
    return (da * 0.5).rio.write_crs(da.rio.crs, inplace=False).rename("elevation_m")


def load_diviner(
    path_max: Path = DEFAULT_RAW_DIR / "diviner" / "diviner_tbol_max_sp.tif",
    path_min: Path = DEFAULT_RAW_DIR / "diviner" / "diviner_tbol_min_sp.tif",
) -> xr.Dataset:
    """Load Diviner annual Tbol max/min as a two-variable Dataset.

    Args:
        path_max: GeoTIFF of annual maximum bolometric temperature (K).
        path_min: GeoTIFF of annual minimum bolometric temperature (K).

    Returns:
        Dataset with variables ``tbol_max`` and ``tbol_min``, each in K.
    """
    tmax = load_raster(path_max).rename("tbol_max")
    tmin = load_raster(path_min).rename("tbol_min")
    return xr.Dataset({"tbol_max": tmax, "tbol_min": tmin})


def load_illumination(
    path: Path = DEFAULT_RAW_DIR / "illumination" / "avgvisib_65s_240m_201608.img",
) -> xr.DataArray:
    """Load the Mazarico average-illumination south-polar product.

    Pixel values are average solar visibility as percent (0–100) over
    the 18.6-year lunar precession cycle. We rescale to a [0, 1]
    fraction here so downstream code is unit-consistent.

    Args:
        path: Path to ``avgvisib_65s_240m_201608.img``. The ``.lbl`` must
            be alongside.

    Returns:
        DataArray of illumination fraction in [0, 1].
    """
    da = load_raster(path) / 100.0
    return da.rename("illumination_fraction")


def load_lend(path: Path = DEFAULT_RAW_DIR / "lend" / "lend_csetn_sp.img") -> xr.DataArray:
    """Load the LEND south-polar epithermal-neutron-flux map.

    Args:
        path: Path to the LEND polar product (URL is TODO-flagged in
            :mod:`selene_base.data.download` until verified).

    Returns:
        DataArray of epithermal neutron count rate or flux in the file's
        native units.
    """
    return load_raster(path).rename("epithermal_neutron_flux")


# ----------------------------------------------------------------------------
# Robbins crater catalog (vector)
# ----------------------------------------------------------------------------
def load_crater_catalog(
    path: Path = DEFAULT_RAW_DIR / "robbins" / "robbins_southpole.csv.gz",
) -> gpd.GeoDataFrame:
    """Load the south-polar Robbins crater slice as a GeoDataFrame.

    Reads the gzipped CSV produced by
    :func:`selene_base.data.download.download_robbins` and constructs
    point geometries in the lunar geographic CRS.

    Args:
        path: Path to ``robbins_southpole.csv.gz``.

    Returns:
        GeoDataFrame with at minimum ``lat``, ``lon``, ``diam_km``, and
        a ``geometry`` column of points in
        :data:`LUNAR_GEOGRAPHIC_CRS`.
    """
    df = pd.read_csv(path)
    geom = gpd.points_from_xy(df["lon"], df["lat"])
    return gpd.GeoDataFrame(df, geometry=geom, crs=LUNAR_GEOGRAPHIC_CRS)
