"""Build the small bundled sample dataset published as a release asset.

Produces ``data/sample/sample_data.tar.gz`` containing downsampled,
trimmed copies of every weeks-1 product the pipeline needs to run end
to end:

- ``sample_data/lola/sample_lola.tif`` — LOLA elevation, downsampled to
  ~1200 m / pixel (5x coarser than working resolution).
- ``sample_data/illumination/sample_illumination.tif`` — Mazarico
  average illumination, same downsampling.
- ``sample_data/robbins/robbins_sample.csv.gz`` — Robbins crater catalog
  filtered to lat ≤ -85°.

This is a developer utility, not user-facing. Run it once, upload the
resulting tarball as a release asset, and the rest of the project
references it by URL. Run again whenever the source data or sample
schema changes.
"""

from __future__ import annotations

import gzip
import shutil
import sys
import tarfile
from pathlib import Path

import pandas as pd
import rioxarray  # noqa: F401

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "sample"
TARBALL = OUT_DIR / "sample_data.tar.gz"
DOWNSAMPLE_FACTOR = 5
SAMPLE_LAT_MAX = -85.0


def _downsample_raster(src: Path, dest: Path) -> None:
    da = rioxarray.open_rasterio(src, masked=True).squeeze("band", drop=True)
    coarse = da.coarsen(y=DOWNSAMPLE_FACTOR, x=DOWNSAMPLE_FACTOR, boundary="trim").mean()
    coarse = coarse.rio.write_crs(da.rio.crs, inplace=False)
    dest.parent.mkdir(parents=True, exist_ok=True)
    coarse.rio.to_raster(dest, driver="GTiff", compress="DEFLATE")


def _filter_robbins(src: Path, dest: Path) -> None:
    df = pd.read_csv(src)
    sp = df[df["lat"] <= SAMPLE_LAT_MAX].copy()
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wt", encoding="utf-8") as gz:
        sp.to_csv(gz, index=False)


def _human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def main() -> int:
    if not RAW.exists():
        print(f"[error] {RAW} does not exist; run `selene download` first.", file=sys.stderr)
        return 2

    sources = {
        "lola": RAW / "lola" / "ldem_80s_80m.lbl",
        "illumination": RAW / "illumination" / "avgvisib_65s_240m_201608.lbl",
        "robbins": RAW / "robbins" / "robbins_southpole.csv.gz",
    }
    for name, path in sources.items():
        if not path.exists():
            print(
                f"[error] missing {name} source at {path}; run `selene download {name}` first.",
                file=sys.stderr,
            )
            return 2

    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    staging = OUT_DIR / "sample_data"

    print("[downsample] LOLA")
    _downsample_raster(sources["lola"], staging / "lola" / "sample_lola.tif")

    print("[downsample] illumination")
    _downsample_raster(
        sources["illumination"], staging / "illumination" / "sample_illumination.tif"
    )

    print("[filter] robbins (lat <= -85)")
    _filter_robbins(sources["robbins"], staging / "robbins" / "robbins_sample.csv.gz")

    manifest = staging / "MANIFEST.txt"
    manifest.write_text(
        "selene-base sample data bundle\n"
        f"downsample factor: {DOWNSAMPLE_FACTOR}x (vs full-resolution working grid)\n"
        f"robbins lat cutoff: {SAMPLE_LAT_MAX}\n"
        "intended for: smoke-testing the pipeline; NOT for analysis.\n",
        encoding="utf-8",
    )

    print(f"[archive] {TARBALL}")
    with tarfile.open(TARBALL, "w:gz") as tar:
        tar.add(staging, arcname=staging.name)

    size = TARBALL.stat().st_size
    print()
    print(f"sample_data.tar.gz: {_human_bytes(size)}")
    with tarfile.open(TARBALL, "r:gz") as tar:
        for member in tar.getmembers():
            kind = "d" if member.isdir() else "f"
            print(f"  {kind} {member.size:>10}  {member.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
