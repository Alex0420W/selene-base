"""Raw data fetchers, one function per source product.

Each function downloads its dataset into ``data/raw/<dataset>/`` and is
idempotent: a second call short-circuits if the file is already on disk
and meets a minimum-size sanity check. Network access, retry, and
progress reporting all go through :mod:`selene_base.data._http`.

Filled in week 1; sample-data path added week 5.
"""

from __future__ import annotations

import gzip
import io
import shutil
import tarfile
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import typer

from selene_base.data._http import stream_to_file

DEFAULT_RAW_DIR = Path("data/raw")

# ----------------------------------------------------------------------------
# 1. Robbins 2018 lunar crater catalog (verified 2026-04 via USGS Astrogeology
#    CKAN dataset f89f5478-b69a-486c-b9b5-30d7b0c5ad2b).
# ----------------------------------------------------------------------------
ROBBINS_URL = (
    "https://astrogeology.usgs.gov/ckan/dataset/"
    "f89f5478-b69a-486c-b9b5-30d7b0c5ad2b/resource/"
    "c4f25cc2-4f8a-4207-a845-5e176da3ac5a/download/lunar_crater_database_robbins_2018"
)
ROBBINS_RAW_BYTES = 50_000_000  # full archive ~92 MB; conservative floor
ROBBINS_FILTERED_LAT_MAX = -75.0  # south-polar slice we keep in repo cache

# ----------------------------------------------------------------------------
# 2. LOLA LDEM south-polar DEM (verified 2026-04 via PDS Geosciences listing).
#    PDS publishes the 80°S+ polar grid at 80 m and 20 m resolutions in the
#    same archive directory. v1.4 used the 80 m product downsampled to a
#    240 m analysis grid; v1.5 fetches the 20 m product for the tiled
#    high-resolution per-region run.
# ----------------------------------------------------------------------------
_LOLA_BASE = (
    "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/polar/img"
)
LOLA_LDEM_IMG_URL = f"{_LOLA_BASE}/ldem_80s_80m.img"
LOLA_LDEM_LBL_URL = f"{_LOLA_BASE}/ldem_80s_80m.lbl"
LOLA_LDEM_IMG_MIN_BYTES = 100_000_000  # ~115 MB
LOLA_LDEM_LBL_MIN_BYTES = 1_000

LOLA_LDEM_20M_IMG_URL = f"{_LOLA_BASE}/ldem_80s_20m.img"
LOLA_LDEM_20M_LBL_URL = f"{_LOLA_BASE}/ldem_80s_20m.lbl"
LOLA_LDEM_20M_IMG_MIN_BYTES = 1_000_000_000  # ~1.85 GB

LOLA_RESOLUTIONS_M: tuple[int, ...] = (80, 20)

# ----------------------------------------------------------------------------
# 3. Diviner Polar Resource Product — south pole (verified week 6).
#    PDS4 Table_Character: 2.88M triangle-mesh records, ~605 MB TAB +
#    14 KB XML. Provides Tavg, Tmax, modeled ice-stability depth at the
#    centre of every triangular facet — drives both the thermal and ice
#    criteria from a single authoritative source.
# ----------------------------------------------------------------------------
DIVINER_PRP_BASE = (
    "https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp"
)
DIVINER_PRP_TAB_URL = f"{DIVINER_PRP_BASE}/dlre_prp_south.tab"
DIVINER_PRP_XML_URL = f"{DIVINER_PRP_BASE}/dlre_prp_south.xml"
DIVINER_PRP_TAB_MIN_BYTES = 600_000_000  # ~605 MB
DIVINER_PRP_XML_MIN_BYTES = 1_000

# ----------------------------------------------------------------------------
# 4. Mazarico average illumination — south-polar 65S 240m product
#    (verified 2026-04 via PDS Geosciences listing). 65S extent fully
#    covers our 80°S–90°S region of interest; resolution matches the
#    target grid.
# ----------------------------------------------------------------------------
_ILLUM_BASE = (
    "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/"
    "lrolol_1xxx/extras/illumination/release_2016/img"
)
ILLUMINATION_IMG_URL = f"{_ILLUM_BASE}/avgvisib_65s_240m_201608.img"
ILLUMINATION_LBL_URL = f"{_ILLUM_BASE}/avgvisib_65s_240m_201608.lbl"
ILLUMINATION_IMG_MIN_BYTES = 70_000_000  # ~82 MB
ILLUMINATION_LBL_MIN_BYTES = 1_000

# ----------------------------------------------------------------------------
# 5. LEND epithermal neutron flux south-polar map.
# TODO(week1): URL not verified — the LEND archive at
#     https://pds-geosciences.wustl.edu/lro/lro-l-lend-2-edr-v1/lrolen_0xxx/
# ships year-binned EDRs and inventory CSVs but does not publish a single
# pre-made polar neutron-flux map under an obvious filename. Original
# spec referenced "polar map of epithermal neutron flux (or count rate)".
# Likely options to confirm: (a) ODE moon search for LEND CSETN derived
# polar product; (b) Sanin/Mitrofanov supplementary maps. Confirm before
# running selene download lend.
# ----------------------------------------------------------------------------
LEND_URL = (
    "https://pds-geosciences.wustl.edu/lro/lro-l-lend-2-edr-v1/"
    "lrolen_0xxx/data_derived/"  # TODO(week1): replace with real polar product URL
    "lend_csetn_sp.img"
)
LEND_MIN_BYTES = 100_000


# ----------------------------------------------------------------------------
# Per-dataset entry points
# ----------------------------------------------------------------------------
def download_robbins(dest: Path = DEFAULT_RAW_DIR / "robbins") -> Path:
    """Download the Robbins (2018) crater catalog and filter to south polar.

    The published archive ships a TAB-separated table of ~1.3M craters.
    We persist two artefacts:

    * ``lunar_crater_database_robbins_2018`` — the full raw download
      (so the slow-step is cached and rerunnable).
    * ``robbins_southpole.csv.gz`` — the lat ≤ -75° slice used by the
      pipeline, gzipped to keep the repo cache small.

    Args:
        dest: Directory to write into. Created if missing.

    Returns:
        Path to the filtered south-polar CSV.
    """
    dest = Path(dest)
    raw_path = dest / "lunar_crater_database_robbins_2018"
    filtered_path = dest / "robbins_southpole.csv.gz"

    if filtered_path.exists() and filtered_path.stat().st_size > 1_000:
        typer.echo(f"[skip] {filtered_path.name} already present")
        return filtered_path

    stream_to_file(
        ROBBINS_URL,
        raw_path,
        min_bytes=ROBBINS_RAW_BYTES,
        label="robbins (raw)",
    )

    typer.echo(f"[filter] keeping rows with LAT_CIRC_IMG <= {ROBBINS_FILTERED_LAT_MAX}")
    df = _read_robbins_table(raw_path)
    lat_col = _robbins_lat_column(df)
    lon_col = _robbins_lon_column(df)
    diam_col = _robbins_diameter_column(df)

    sp = df[df[lat_col] <= ROBBINS_FILTERED_LAT_MAX].copy()
    sp = sp.rename(columns={lat_col: "lat", lon_col: "lon", diam_col: "diam_km"})
    keep = ["lat", "lon", "diam_km"]
    if "CRATER_ID" in sp.columns:
        sp = sp.rename(columns={"CRATER_ID": "crater_id"})
        keep = ["crater_id", *keep]

    with gzip.open(filtered_path, "wt", encoding="utf-8") as gz:
        sp[keep].to_csv(gz, index=False)
    typer.echo(f"[done] {filtered_path} ({len(sp):,} rows)")
    return filtered_path


def download_lola(
    dest: Path = DEFAULT_RAW_DIR / "lola",
    *,
    resolution_m: int = 80,
) -> Path:
    """Download the LOLA LDEM south-polar DEM (PDS3 IMG + LBL).

    The PDS Geosciences archive publishes the 80°S+ polar grid at both
    80 m and 20 m. v1.4 ran on the 80 m product downsampled to a 240 m
    analysis grid; v1.5 uses the 20 m product for the per-region tiled
    high-resolution run.

    Args:
        dest: Directory to write into.
        resolution_m: Native LOLA grid resolution to fetch. Must be one
            of :data:`LOLA_RESOLUTIONS_M` (80 or 20). 80 m is the v1.4
            default; 20 m is the v1.5 product (~1.85 GB).

    Returns:
        Path to the downloaded directory.

    Raises:
        ValueError: If ``resolution_m`` is not supported.
    """
    if resolution_m not in LOLA_RESOLUTIONS_M:
        raise ValueError(
            f"resolution_m={resolution_m!r} not in supported set {LOLA_RESOLUTIONS_M!r}"
        )
    dest = Path(dest)
    if resolution_m == 80:
        img_url, lbl_url = LOLA_LDEM_IMG_URL, LOLA_LDEM_LBL_URL
        img_min, lbl_min = LOLA_LDEM_IMG_MIN_BYTES, LOLA_LDEM_LBL_MIN_BYTES
    else:  # 20 m
        img_url, lbl_url = LOLA_LDEM_20M_IMG_URL, LOLA_LDEM_20M_LBL_URL
        img_min, lbl_min = LOLA_LDEM_20M_IMG_MIN_BYTES, LOLA_LDEM_LBL_MIN_BYTES
    stream_to_file(
        img_url,
        dest / f"ldem_80s_{resolution_m}m.img",
        min_bytes=img_min,
        label=f"lola img ({resolution_m} m)",
    )
    stream_to_file(
        lbl_url,
        dest / f"ldem_80s_{resolution_m}m.lbl",
        min_bytes=lbl_min,
        label=f"lola lbl ({resolution_m} m)",
    )
    return dest


def download_diviner(dest: Path = DEFAULT_RAW_DIR / "diviner") -> Path:
    """Download the Diviner Polar Resource Product for the south pole.

    Fetches both ``dlre_prp_south.tab`` (~605 MB) and
    ``dlre_prp_south.xml`` (~14 KB) from the PDS Geosciences Diviner
    derived bundle. The PRP is a PDS4 character table indexed by
    triangular-mesh facet centre; downstream parsing happens in
    :mod:`selene_base.data.pds4_table`.

    Args:
        dest: Directory to write into. Created if missing.

    Returns:
        Path to the downloaded directory.
    """
    dest = Path(dest)
    stream_to_file(
        DIVINER_PRP_XML_URL,
        dest / "dlre_prp_south.xml",
        min_bytes=DIVINER_PRP_XML_MIN_BYTES,
        label="diviner prp xml",
    )
    stream_to_file(
        DIVINER_PRP_TAB_URL,
        dest / "dlre_prp_south.tab",
        min_bytes=DIVINER_PRP_TAB_MIN_BYTES,
        label="diviner prp tab",
    )
    return dest


def download_illumination(dest: Path = DEFAULT_RAW_DIR / "illumination") -> Path:
    """Download the Mazarico average-illumination south-polar product.

    Args:
        dest: Directory to write into.

    Returns:
        Path to the downloaded directory.
    """
    dest = Path(dest)
    stream_to_file(
        ILLUMINATION_IMG_URL,
        dest / "avgvisib_65s_240m_201608.img",
        min_bytes=ILLUMINATION_IMG_MIN_BYTES,
        label="illum img",
    )
    stream_to_file(
        ILLUMINATION_LBL_URL,
        dest / "avgvisib_65s_240m_201608.lbl",
        min_bytes=ILLUMINATION_LBL_MIN_BYTES,
        label="illum lbl",
    )
    return dest


def download_scarps(dest: Path = DEFAULT_RAW_DIR / "scarps") -> Path:
    """Download Watters et al. lobate-scarp catalog for the lunar south pole.

    .. note::
       URL is not verified. The Watters catalog is published as
       supplementary material across multiple papers (2010, 2015,
       2019); locate or assemble a consolidated south-polar version,
       or contact the LROC team for the working catalog.

       Place the file at ``data/raw/scarps/scarps_southpole.geojson``
       (or ``.csv``) with at minimum ``geometry`` (line/multiline) plus
       ``length_km`` and ``scarp_id`` columns.

    Args:
        dest: Directory to write into.

    Returns:
        Path to the downloaded catalog.
    """
    raise NotImplementedError(
        "Watters scarp catalog: source URL not verified. See docstring for the "
        "manual placement path the seismic criterion expects."
    )


def download_lend(dest: Path = DEFAULT_RAW_DIR / "lend") -> Path:
    """Download the LEND south-polar epithermal-neutron map.

    .. note::
       URL is TODO-flagged; confirm the polar-map filename via the PDS
       LEND derived archive or the Lunar ODE before first run.

    Args:
        dest: Directory to write into.

    Returns:
        Path to the downloaded directory.
    """
    dest = Path(dest)
    stream_to_file(
        LEND_URL,
        dest / "lend_csetn_sp.img",
        min_bytes=LEND_MIN_BYTES,
        label="lend",
    )
    return dest


# ----------------------------------------------------------------------------
# Aggregate driver
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Sample-data bundle (week 5).
# ----------------------------------------------------------------------------
# TODO(week5): URL valid only after the maintainer creates the
# 'sample-data-v1' GitHub release with sample_data.tar.gz attached.
# Run scripts/build_sample_data.py to regenerate the tarball.
SAMPLE_DATA_URL = (
    "https://github.com/Alex0420W/selene-base/releases/download/sample-data-v1/sample_data.tar.gz"
)
SAMPLE_DATA_MIN_BYTES = 1_000_000  # ~12 MB; floor at 1 MB
_SAMPLE_MARKER = ".sample_extracted"


def download_sample_data(
    raw_dir: Path = DEFAULT_RAW_DIR,
    *,
    url: str = SAMPLE_DATA_URL,
) -> Path:
    """Download and extract the bundled sample dataset.

    Fetches ``sample_data.tar.gz`` from the GitHub release asset,
    extracts its contents alongside the regular per-dataset folders so
    the rest of the pipeline finds them without any extra config:

    - ``data/raw/lola/sample_lola.tif``
    - ``data/raw/illumination/sample_illumination.tif``
    - ``data/raw/robbins/robbins_southpole.csv.gz``

    Idempotent: skips if the marker file
    ``data/raw/.sample_extracted`` exists.

    Args:
        raw_dir: Destination root (defaults to ``data/raw``).
        url: Override URL (used by tests).

    Returns:
        ``raw_dir`` once extraction is complete.

    Raises:
        AssertionError: If the downloaded archive is smaller than
            :data:`SAMPLE_DATA_MIN_BYTES`.
    """
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    marker = raw_dir / _SAMPLE_MARKER
    if marker.exists():
        typer.echo(f"[skip] sample data already extracted ({marker})")
        return raw_dir

    archive_path = raw_dir / "sample_data.tar.gz"
    stream_to_file(
        url,
        archive_path,
        min_bytes=SAMPLE_DATA_MIN_BYTES,
        label="sample_data.tar.gz",
    )

    typer.echo(f"[extract] {archive_path}")
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            # Strip the top-level "sample_data/" prefix so the inner
            # subdirs land directly under raw_dir.
            parts = Path(member.name).parts
            if not parts:
                continue
            relative = Path(*parts[1:]) if parts[0] == "sample_data" else Path(*parts)
            dest = raw_dir / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            with dest.open("wb") as out:
                shutil.copyfileobj(extracted, out)
            typer.echo(f"  -> {dest}")

    # Robbins sample is stored as robbins_sample.csv.gz; canonicalise the
    # name the rest of the pipeline expects.
    rb_sample = raw_dir / "robbins" / "robbins_sample.csv.gz"
    rb_canonical = raw_dir / "robbins" / "robbins_southpole.csv.gz"
    if rb_sample.exists() and not rb_canonical.exists():
        rb_sample.rename(rb_canonical)
        typer.echo(f"  -> {rb_canonical} (canonical robbins filename)")

    marker.write_text("sample-data-v1\n", encoding="utf-8")
    archive_path.unlink(missing_ok=True)
    typer.echo("[done] sample data installed; run `selene preprocess` next.")
    return raw_dir


DATASETS: dict[str, Callable[[Path], Path]] = {
    "robbins": download_robbins,
    "lola": download_lola,
    "diviner": download_diviner,
    "illumination": download_illumination,
    "lend": download_lend,
    "scarps": download_scarps,
}


def download_all(dest: Path = DEFAULT_RAW_DIR) -> dict[str, Path]:
    """Run every per-dataset downloader and return a name → path mapping.

    Args:
        dest: Parent directory under which each dataset gets its own folder.

    Returns:
        Mapping from dataset name to local directory or filtered file.
    """
    dest = Path(dest)
    out: dict[str, Path] = {}
    for name, fn in DATASETS.items():
        out[name] = fn(dest / name)
    return out


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def _read_robbins_table(path: Path) -> pd.DataFrame:
    """Robbins ships either a CSV directly or a PDS4 bundle ZIP.

    Inside the ZIP the actual ~240 MB data table is the largest CSV
    member; smaller CSVs are inventory manifests.
    """
    with path.open("rb") as fh:
        head = fh.read(4)
    if head.startswith(b"PK\x03\x04"):
        import zipfile

        with zipfile.ZipFile(path) as zf:
            csv_members = [
                info
                for info in zf.infolist()
                if info.filename.lower().endswith(".csv") and info.file_size > 0
            ]
            if not csv_members:
                raise RuntimeError(f"no CSV member found inside {path}")
            largest = max(csv_members, key=lambda i: i.file_size)
            typer.echo(f"[zip] reading {largest.filename} ({largest.file_size:,} bytes)")
            with zf.open(largest) as inner:
                return pd.read_csv(io.TextIOWrapper(inner, encoding="utf-8"))
    return pd.read_csv(path)


def _robbins_lat_column(df: pd.DataFrame) -> str:
    return _first_present(df, ["LAT_CIRC_IMG", "LAT_ELLI_IMG", "LAT", "LATITUDE_CIRCLE_IMAGE"])


def _robbins_lon_column(df: pd.DataFrame) -> str:
    return _first_present(df, ["LON_CIRC_IMG", "LON_ELLI_IMG", "LON", "LONGITUDE_CIRCLE_IMAGE"])


def _robbins_diameter_column(df: pd.DataFrame) -> str:
    return _first_present(df, ["DIAM_CIRC_IMG", "DIAM_ELLI_IMG", "DIAMETER", "DIAM_CIRCLE_IMAGE"])


def _first_present(df: pd.DataFrame, candidates: list[str]) -> str:
    for name in candidates:
        if name in df.columns:
            return name
    raise KeyError(f"none of {candidates} present in columns: {list(df.columns)[:20]}")
