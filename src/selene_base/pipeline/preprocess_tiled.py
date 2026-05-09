"""Per-region tiled preprocessing for high-resolution analysis (v1.5).

The v1.4 pipeline reprojects every raw raster onto a single 240 m
analysis grid covering the whole 80°S+ region. At 240 m the
``los_to_earth`` horizon profile fits in memory as one ``(36, 2533,
2533)`` float32 array (~0.9 GB) and the 100-km ray-march finishes in a
few seconds on the GPU.

At 20 m the same global grid would be ``(36, 30400, 30400)`` ≈ 132 GB —
larger than even unified memory on a GB10. v1.5 sidesteps this by
processing each USGS region as its own tile: bounding box of the
polygon plus a 100 km buffer (matching :data:`DEFAULT_MAX_HORIZON_KM`,
so the ray-march sees the same physical horizon as the 240 m run did),
windowed and reprojected to a local 20 m polar-stereographic grid, then
fed through :func:`derive_horizon_profile` on the GPU. Tiles are
processed sequentially with explicit GPU memory release between regions.

Polygons in the central polar cluster (Mons Mouton, MMP, Malapert
Massif, Haworth, Nobile Rim 1, Nobile Rim 2) produce overlapping tiles;
recomputing the overlapping pixels is intentional — keeping each tile
self-contained avoids stitching artefacts at tile boundaries and the
total wall-clock cost is dominated by the largest (MMP) tile anyway.

Cache layout::

    data/processed/horizon_profile_southpole_<resolution_m>m_<region_code>.npz

with ``region_code`` lower-cased (``sp``, ``mm``, ``mmp``, ...). The
240 m global path in :mod:`selene_base.pipeline.preprocess` is
unchanged; tiled mode is opt-in via the ``--tiled-per-region`` CLI flag.

Unified-memory + per-azimuth chunking + mmap accumulator (v1.9). At
10 m + 100 km buffer, per-tile dimensions are ~22000×22000 = 484 M
pixels; the full ``(36, height, width)`` float32 horizon-profile
output array is ~70 GB. On unified-memory hardware (NVIDIA GB10 / DGX
Spark, 119 GB shared between CPU and GPU) this exceeds the
unified-memory cap once scratch buffers + CuPy mempool overhead are
factored in, and the OOM killer takes the entire process group
(including the parent tmux session). A pure host-side accumulator
"spills" the GPU output back to RAM but on unified memory there is no
separate large host pool — the host array still consumes the same
unified pool the GPU draws from, so per-azimuth chunking alone is
insufficient on this architecture. v1.9 closes this with a
mmap-backed accumulator: at fine resolutions the per-tile output
array is opened with ``np.lib.format.open_memmap`` against a
temporary on-disk ``.npy`` file, and the chunked horizon-profile
loop in :func:`selene_base.criteria.los_to_earth.derive_horizon_profile`
writes one azimuth-chunk slice at a time. Peak in-memory footprint
is then bounded by the per-chunk GPU buffer (~17 GB at chunk=9, 10 m)
plus the elevation tile (~2 GB) plus mempool overhead, regardless of
the per-tile output size. After the chunked compute completes,
``np.savez_compressed`` streams the mmap to the final NPZ cache file
in 16 MB nditer chunks (numpy's default buffer size); the temporary
``.npy`` is then unlinked. At 240 m / 80 m / 20 m the auto-detected
chunk size is ``n_azimuths`` (no chunking; v1.5 behaviour preserved
verbatim).
"""

from __future__ import annotations

import gc
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rioxarray  # noqa: F401  (registers .rio accessor)
import typer
import xarray as xr

from selene_base.criteria import los_to_earth
from selene_base.data.load import LUNAR_SOUTH_POLAR_CRS, load_lola_ldem, load_raster
from selene_base.data.reproject import reproject_to_grid
from selene_base.validation.nasa_regions import regions_polygons_to_geodataframe

DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_BUFFER_M = 100_000.0  # matches los_to_earth.DEFAULT_MAX_HORIZON_KM


# Order LOLA candidate sources finest → coarsest so a 20 m IMG is
# preferred when both are present. Keys are the on-disk filenames; the
# value is the native pixel resolution (used only for logging).
_LOLA_SOURCE_PRIORITY: tuple[tuple[str, int], ...] = (
    ("ldem_80s_20m.lbl", 20),
    ("ldem_80s_80m.lbl", 80),
    ("sample_lola.tif", 240),
)


@dataclass(frozen=True)
class TileSpec:
    """One per-region tile: polygon bbox + buffer in the polar-stereo CRS."""

    region_name: str
    region_code: str  # lower-case, e.g. 'sp' for Slater Plain
    xmin: float
    ymin: float
    xmax: float
    ymax: float

    @property
    def width_m(self) -> float:
        return self.xmax - self.xmin

    @property
    def height_m(self) -> float:
        return self.ymax - self.ymin

    def shape(self, resolution_m: float) -> tuple[int, int]:
        """Pixel ``(height, width)`` for this tile at the given resolution."""
        return (
            int(round(self.height_m / resolution_m)),
            int(round(self.width_m / resolution_m)),
        )


@dataclass
class TiledHorizonResult:
    """One row of the per-region tiled summary."""

    region_name: str
    region_code: str
    status: str  # 'cached', 'skipped', 'failed'
    output_path: Path | None
    n_pixels: int
    elapsed_s: float


_VALID_LOLA_SUFFIXES: frozenset[str] = frozenset({".tif", ".lbl", ".img"})


def resolve_lola_source(
    raw_dir: Path = DEFAULT_RAW_DIR,
    *,
    prefer_resolution_m: int | None = None,
    override: Path | str | None = None,
) -> Path:
    """Pick the best available LOLA source IMG/LBL/TIF on disk.

    The on-disk priority is finest-resolution first (``ldem_80s_20m.lbl``
    → ``ldem_80s_80m.lbl`` → ``sample_lola.tif``). Pass
    ``prefer_resolution_m`` to require a specific resolution; the call
    raises :class:`FileNotFoundError` if that file is missing rather
    than silently falling back to a coarser source. Pass ``override``
    to bypass the priority list entirely with a caller-supplied path —
    used for non-PDS-named sources like PGDA mosaics
    (``ldem_83s_10mpp_adj.tif``) and arbitrary user-supplied DEMs.

    Args:
        raw_dir: Root of ``data/raw``.
        prefer_resolution_m: Require this native resolution (m). When
            ``None``, any available source is acceptable. Ignored when
            ``override`` is supplied.
        override: Explicit path to a LOLA source file. Bypasses both
            the priority list and ``prefer_resolution_m``. The file
            must exist and have a recognised suffix
            (``.tif``/``.lbl``/``.img``).

    Returns:
        Absolute path to the LOLA source the loaders should open.

    Raises:
        FileNotFoundError: If no acceptable source is on disk, or if
            ``override`` is supplied but the file does not exist.
        ValueError: If ``override`` is supplied with an unrecognised
            file suffix.
    """
    if override is not None:
        path = Path(override)
        if not path.exists():
            raise FileNotFoundError(f"--lola-source path does not exist: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"--lola-source path is not a file: {path}")
        if path.suffix.lower() not in _VALID_LOLA_SUFFIXES:
            raise ValueError(
                f"--lola-source suffix {path.suffix!r} not recognised; "
                f"expected one of {sorted(_VALID_LOLA_SUFFIXES)}"
            )
        return path
    lola_dir = Path(raw_dir) / "lola"
    if prefer_resolution_m is not None:
        path = lola_dir / f"ldem_80s_{prefer_resolution_m}m.lbl"
        if not path.exists():
            raise FileNotFoundError(
                f"LOLA {prefer_resolution_m} m source not found at {path}; "
                f"run `selene download lola --resolution {prefer_resolution_m}m` first."
            )
        return path
    for filename, _res in _LOLA_SOURCE_PRIORITY:
        candidate = lola_dir / filename
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No LOLA source found under {lola_dir}; expected one of "
        f"{[name for name, _ in _LOLA_SOURCE_PRIORITY]}."
    )


def _load_lola_source(path: Path) -> xr.DataArray:
    """Load a LOLA source raster, applying the 0.5 m PDS3 scale when needed."""
    if path.suffix.lower() == ".tif":
        return load_raster(path).rename("elevation_m")
    return load_lola_ldem(path)


def compute_tile_specs(
    target_crs: str = str(LUNAR_SOUTH_POLAR_CRS),
    *,
    buffer_m: float = DEFAULT_BUFFER_M,
    resolution_m: float = 20.0,
    region_codes: Iterable[str] | None = None,
) -> list[TileSpec]:
    """Compute per-USGS-region tile bounding boxes in the analysis CRS.

    The bbox of each USGS polygon is expanded by ``buffer_m`` on every
    side and snapped to a multiple of ``resolution_m`` so the tile aligns
    with the analysis grid. The buffer defaults to 100 km, matching
    :data:`selene_base.criteria.los_to_earth.DEFAULT_MAX_HORIZON_KM`, so
    a pixel anywhere inside the polygon can ray-march the full 100 km
    without leaving the tile.

    Args:
        target_crs: CRS of the output bbox. Defaults to the v1.4
            south-polar stereographic analysis CRS.
        buffer_m: Buffer to add on every side of the polygon bbox.
        resolution_m: Pixel size to snap the tile bounds to. Snapping
            ensures the tile origin lands on an integer multiple of the
            resolution, which keeps adjacent tiles aligned to a
            consistent global grid.
        region_codes: When given, restrict the output to USGS regions
            whose ``RegionCode`` (case-insensitive) is in this set.
            ``None`` returns specs for all 9 NASA candidate regions.

    Returns:
        One :class:`TileSpec` per polygon, in USGS file order.
    """
    if buffer_m < 0:
        raise ValueError(f"buffer_m must be non-negative, got {buffer_m!r}")
    if resolution_m <= 0:
        raise ValueError(f"resolution_m must be positive, got {resolution_m!r}")

    polygons = regions_polygons_to_geodataframe(target_crs=target_crs)
    if region_codes is not None:
        wanted = {c.upper() for c in region_codes}
        polygons = polygons[polygons["RegionCode"].str.upper().isin(wanted)]

    specs: list[TileSpec] = []
    for _, row in polygons.iterrows():
        xmin, ymin, xmax, ymax = row.geometry.bounds
        # Expand by buffer, then snap each edge outward to a multiple of
        # resolution_m so adjacent tiles share a common grid origin.
        xmin = float(np.floor((xmin - buffer_m) / resolution_m) * resolution_m)
        ymin = float(np.floor((ymin - buffer_m) / resolution_m) * resolution_m)
        xmax = float(np.ceil((xmax + buffer_m) / resolution_m) * resolution_m)
        ymax = float(np.ceil((ymax + buffer_m) / resolution_m) * resolution_m)
        specs.append(
            TileSpec(
                region_name=str(row["Region"]),
                region_code=str(row["RegionCode"]).lower(),
                xmin=xmin,
                ymin=ymin,
                xmax=xmax,
                ymax=ymax,
            )
        )
    return specs


def reproject_to_tile(
    src: xr.DataArray,
    tile: TileSpec,
    *,
    target_crs: str = str(LUNAR_SOUTH_POLAR_CRS),
    resolution_m: float = 20.0,
    resampling: str = "bilinear",
) -> xr.DataArray:
    """Window + reproject a source raster onto one tile's local grid.

    Wraps :func:`reproject_to_grid` with the tile's bbox and the target
    resolution. The output grid is rectangular in ``target_crs`` with
    pixel size ``resolution_m`` and the same axis convention as the v1.4
    240 m grid (top-left origin, y descending).
    """
    return reproject_to_grid(
        src,
        target_crs=target_crs,
        bounds_m=(tile.xmin, tile.ymin, tile.xmax, tile.ymax),
        resolution_m=resolution_m,
        resampling=resampling,
    )


def horizon_npz_path(
    processed_dir: Path,
    *,
    resolution_m: float,
    region_code: str,
) -> Path:
    """Per-tile horizon-profile cache filename for one region (legacy NPZ)."""
    return (
        Path(processed_dir)
        / f"horizon_profile_southpole_{int(round(resolution_m))}m_{region_code.lower()}.npz"
    )


def horizon_cache_paths(
    processed_dir: Path,
    *,
    resolution_m: float,
    region_code: str,
) -> tuple[Path, Path, Path]:
    """Return ``(npz, npy, meta_npz)`` cache-path candidates for one region.

    The unchunked v1.5 path writes a single ``<base>.npz`` containing the
    horizon profile and all metadata. The v1.9 chunked + mmap path writes
    ``<base>.npy`` (the bulk float32 array, mmap-loadable) plus
    ``<base>.meta.npz`` (the small metadata fields). Callers detect which
    pair is present on disk via :func:`existing_horizon_cache_path`.
    """
    base = (
        Path(processed_dir)
        / f"horizon_profile_southpole_{int(round(resolution_m))}m_{region_code.lower()}"
    )
    return (
        base.with_suffix(".npz"),
        base.with_suffix(".npy"),
        base.with_name(base.name + ".meta.npz"),
    )


def existing_horizon_cache_path(
    processed_dir: Path,
    *,
    resolution_m: float,
    region_code: str,
) -> Path | None:
    """Return whichever cache file exists on disk: ``.npy`` (mmap path,
    accompanied by ``.meta.npz``), ``.npz`` (legacy single-file path), or
    ``None`` if neither is present.

    The mmap path is preferred when both happen to exist (the rank flow
    would otherwise materialise the full 70 GB compressed array on
    unified-memory hardware). Callers writing a fresh cache should clear
    any stale partner file to avoid format ambiguity.
    """
    npz_path, npy_path, meta_path = horizon_cache_paths(
        processed_dir, resolution_m=resolution_m, region_code=region_code
    )
    if npy_path.exists() and meta_path.exists():
        return npy_path
    if npz_path.exists():
        return npz_path
    return None


def _free_gpu_memory() -> None:
    """Release pooled GPU allocations between tiles. No-op if cupy is absent."""
    try:
        import cupy as cp  # type: ignore[import-not-found]

        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass
    gc.collect()


_DEFAULT_PER_CHUNK_BUDGET_BYTES: int = 18 * 1024**3  # ~18 GiB per-chunk GPU buffer cap.
# Sized so a 22000×22000 tile (10 m + 100 km buffer) yields chunk=9, and a
# 44000×44000 tile (5 m + 100 km buffer) yields chunk=2. v1.5's 20 m tiles
# fit unchunked for 8 of 9 polygons; the largest (Mons Mouton Plateau,
# 12896×13846) is split into two chunks of 27 + 9 azimuths — a small extra
# round-trip through the temp .npy that costs <30 s on NVMe per tile.


def default_azimuth_chunk_size(
    *,
    n_azimuths: int,
    height: int,
    width: int,
    bytes_per_value: int = 4,
    budget_bytes: int = _DEFAULT_PER_CHUNK_BUDGET_BYTES,
) -> int | None:
    """Pick a per-tile per-azimuth chunk size that fits in the GPU memory budget.

    Returns ``None`` (i.e. no chunking) when the full
    ``(n_azimuths, height, width)`` working buffer fits inside
    ``budget_bytes`` — which is the v1.5 case at every resolution we ship
    by default (240 m / 80 m / 20 m). At 10 m + 100 km buffer on a
    22000×22000 tile, the full 36-azimuth buffer is ~70 GB and exceeds
    the 20 GB default cap, so the function returns 9 (4 chunks per
    polygon-tile, ~17 GB working buffer per chunk). At 5 m it would
    return 2 (16 chunks per polygon-tile, ~16 GB working buffer per
    chunk) — supplying the same ``derive_horizon_profile`` machinery for
    finer-resolution work that becomes data-feasible.

    Args:
        n_azimuths: Number of azimuthal directions in the horizon profile
            (36 by default).
        height: Per-tile pixel height.
        width: Per-tile pixel width.
        bytes_per_value: Bytes per output sample (4 for float32).
        budget_bytes: Maximum acceptable per-chunk working-buffer size,
            in bytes. Defaults to ~20 GiB which sits comfortably under
            the GB10's 119 GB unified memory cap once the elevation tile
            and ray-march scratch buffers are factored in.
    """
    full_bytes = n_azimuths * height * width * bytes_per_value
    if full_bytes <= budget_bytes:
        return None
    per_az_bytes = max(1, height * width * bytes_per_value)
    chunk = max(1, budget_bytes // per_az_bytes)
    return min(chunk, n_azimuths)


def run_tiled_per_region(
    *,
    resolution_m: float = 20.0,
    region_codes: Iterable[str] | None = None,
    processed_dir: Path = DEFAULT_PROCESSED_DIR,
    raw_dir: Path = DEFAULT_RAW_DIR,
    overwrite: bool = False,
    use_gpu: bool = True,
    target_crs: str = str(LUNAR_SOUTH_POLAR_CRS),
    buffer_m: float = DEFAULT_BUFFER_M,
    source_path: Path | None = None,
    azimuth_chunk_size: int | None = None,
    echo: Callable[[str], None] = typer.echo,
) -> list[TiledHorizonResult]:
    """Compute the Earth-LOS horizon profile for each USGS region tile.

    Sequentially:

    1. Load the LOLA source DEM at its native resolution (20 m for v1.5,
       80 m for the v1.4-style fallback).
    2. For each USGS polygon: window + reproject to a local
       ``resolution_m`` grid covering the polygon bbox + ``buffer_m``;
       run :func:`derive_horizon_profile` on the GPU; save the result
       as a compressed NPZ; release GPU memory; move to the next tile.

    Args:
        resolution_m: Output pixel size of the per-tile analysis grid.
        region_codes: When given, restrict the run to USGS regions whose
            ``RegionCode`` is in this set (case-insensitive). ``None``
            processes all 9 NASA candidate regions.
        processed_dir: Output directory for the per-tile NPZ files.
        raw_dir: Root of ``data/raw``. Used to resolve the LOLA source
            when ``source_path`` is not given.
        overwrite: When ``False`` (default), skip tiles whose cache
            file already exists.
        use_gpu: When ``True`` (default), run the inner ray-march on
            the GPU via CuPy. Falls back to NumPy if ``False``.
        target_crs: Output CRS for every tile. Default is the v1.4
            south-polar stereographic CRS so per-tile outputs land on a
            grid that is compatible with the global analysis frame.
        buffer_m: Polygon-bbox buffer applied on every side.
        source_path: Override the LOLA source IMG/LBL/TIF path. When
            ``None``, the function picks the finest-resolution source
            present under ``raw_dir/lola``.
        azimuth_chunk_size: Per-azimuth working-buffer size for the
            chunked + mmap-backed horizon-profile path (v1.9). When
            ``None`` (default), :func:`default_azimuth_chunk_size` is
            consulted per tile: at v1.5's resolutions (240 m / 80 m /
            20 m) the function returns ``None`` (no chunking;
            single-pass behaviour preserved verbatim), at 10 m it
            returns 9, at 5 m it returns 2. Pass an explicit positive
            integer < ``n_azimuths`` to override; pass exactly
            ``n_azimuths`` to force the unchunked path even at fine
            resolutions (useful for benchmarking on hardware with
            separate large host RAM, where the mmap path would be
            unnecessary overhead).
        echo: Logging sink (defaults to ``typer.echo``).

    Returns:
        One :class:`TiledHorizonResult` per processed (or skipped) tile.
    """
    import time

    processed_dir = Path(processed_dir)
    processed_dir.mkdir(parents=True, exist_ok=True)
    src_path = resolve_lola_source(raw_dir, override=source_path)
    echo(f"[load] LOLA source from {src_path}")
    elevation = _load_lola_source(src_path)
    echo(f"[load] elevation shape={elevation.shape!r} dtype={elevation.dtype!s}")

    specs = compute_tile_specs(
        target_crs=target_crs,
        buffer_m=buffer_m,
        resolution_m=resolution_m,
        region_codes=region_codes,
    )
    echo(f"[plan] {len(specs)} tile(s) at {resolution_m:g} m + {buffer_m / 1000:g} km buffer")

    results: list[TiledHorizonResult] = []
    for spec in specs:
        npz_path, npy_path, meta_path = horizon_cache_paths(
            processed_dir, resolution_m=resolution_m, region_code=spec.region_code
        )
        existing = existing_horizon_cache_path(
            processed_dir, resolution_m=resolution_m, region_code=spec.region_code
        )
        if existing is not None and not overwrite:
            size = existing.stat().st_size
            echo(f"[skip] {spec.region_name} ({spec.region_code}): {existing.name} cached")
            results.append(
                TiledHorizonResult(
                    region_name=spec.region_name,
                    region_code=spec.region_code,
                    status="cached",
                    output_path=existing,
                    n_pixels=0,
                    elapsed_s=0.0,
                )
            )
            continue

        height_px, width_px = spec.shape(resolution_m)
        n_pixels = height_px * width_px
        echo(
            f"[tile] {spec.region_name} ({spec.region_code}) "
            f"bbox=[{spec.xmin / 1000:.0f}, {spec.ymin / 1000:.0f}, "
            f"{spec.xmax / 1000:.0f}, {spec.ymax / 1000:.0f}] km "
            f"shape={height_px}x{width_px} ({n_pixels:,} px)"
        )
        t0 = time.perf_counter()

        windowed = reproject_to_tile(
            elevation,
            spec,
            target_crs=target_crs,
            resolution_m=resolution_m,
            resampling="bilinear",
        )
        windowed = windowed.rename("elevation_m")

        if azimuth_chunk_size is not None:
            chunk_size = azimuth_chunk_size
        else:
            chunk_size = default_azimuth_chunk_size(
                n_azimuths=los_to_earth.DEFAULT_N_AZIMUTHS,
                height=height_px,
                width=width_px,
            )

        if chunk_size is not None and chunk_size < los_to_earth.DEFAULT_N_AZIMUTHS:
            echo(
                f"[tile] {spec.region_code}: derive_horizon_profile "
                f"(use_gpu={use_gpu}, pixel_size_m={resolution_m:g}, "
                f"azimuth_chunk_size={chunk_size}, mmap accumulator on disk)"
            )
            # Chunked path produces a plain mmap-able .npy plus a small
            # .meta.npz sidecar instead of a single compressed .npz. The
            # rank stage's `_load_horizon_profile_npz` would otherwise
            # decompress the full ~70 GB array into RAM at load time and
            # OOM on unified-memory hardware (the same physical pool the
            # GPU draws from). Saving the bulk array uncompressed lets
            # the rank flow ``np.load(..., mmap_mode='r')`` and slice
            # row-chunks on demand. Bytes-on-disk go up vs the
            # ``np.savez_compressed`` v1.5 path (10–15 % compression on
            # high-entropy horizon angles) but at this scale disk is
            # cheaper than unified memory.
            npy_path.parent.mkdir(parents=True, exist_ok=True)
            for stale in (npz_path, npy_path, meta_path):
                if stale.exists():
                    stale.unlink()
            mmap_arr = np.lib.format.open_memmap(
                npy_path,
                mode="w+",
                dtype=np.float32,
                shape=(los_to_earth.DEFAULT_N_AZIMUTHS, height_px, width_px),
            )
            try:
                horizon = los_to_earth.derive_horizon_profile(
                    windowed,
                    pixel_size_m=resolution_m,
                    use_gpu=use_gpu,
                    azimuth_chunk_size=chunk_size,
                    out=mmap_arr,
                )
                azimuth_deg = horizon["azimuth"].to_numpy().astype(np.float32)
            finally:
                # Drop the mmap reference before any further file-level
                # operations to release the OS handle cleanly.
                del mmap_arr
            np.savez_compressed(
                meta_path,
                azimuth_deg=azimuth_deg,
                x=windowed["x"].to_numpy().astype(np.float64),
                y=windowed["y"].to_numpy().astype(np.float64),
                tile_bounds_m=np.asarray(
                    [spec.xmin, spec.ymin, spec.xmax, spec.ymax], dtype=np.float64
                ),
                resolution_m=np.float64(resolution_m),
                region_name=np.asarray(spec.region_name),
                region_code=np.asarray(spec.region_code),
            )
            out_path = npy_path
        else:
            echo(
                f"[tile] {spec.region_code}: derive_horizon_profile "
                f"(use_gpu={use_gpu}, pixel_size_m={resolution_m:g})"
            )
            horizon = los_to_earth.derive_horizon_profile(
                windowed, pixel_size_m=resolution_m, use_gpu=use_gpu
            )
            for stale in (npy_path, meta_path):
                if stale.exists():
                    stale.unlink()
            np.savez_compressed(
                npz_path,
                horizon_profile_deg=horizon.to_numpy().astype(np.float32),
                azimuth_deg=horizon["azimuth"].to_numpy().astype(np.float32),
                x=windowed["x"].to_numpy().astype(np.float64),
                y=windowed["y"].to_numpy().astype(np.float64),
                tile_bounds_m=np.asarray(
                    [spec.xmin, spec.ymin, spec.xmax, spec.ymax], dtype=np.float64
                ),
                resolution_m=np.float64(resolution_m),
                region_name=np.asarray(spec.region_name),
                region_code=np.asarray(spec.region_code),
            )
            out_path = npz_path
        elapsed = time.perf_counter() - t0
        size = out_path.stat().st_size
        echo(f"[done] {spec.region_code}: {out_path.name} ({size:,} bytes, {elapsed:.1f} s)")
        results.append(
            TiledHorizonResult(
                region_name=spec.region_name,
                region_code=spec.region_code,
                status="cached",
                output_path=out_path,
                n_pixels=n_pixels,
                elapsed_s=elapsed,
            )
        )

        # Drop the windowed / horizon arrays before the next tile so the
        # GPU pool can be flushed. The next iteration will re-window from
        # the source raster which stays resident.
        del horizon, windowed
        _free_gpu_memory()

    return results


def format_summary(results: list[TiledHorizonResult]) -> str:
    """Render a fixed-width summary of tiled-preprocess results."""
    header = f"{'region':<22} {'code':<5} {'status':<8} {'n_pix':>14} {'sec':>8}  path"
    rows = []
    for r in results:
        path_str = str(r.output_path) if r.output_path else "-"
        rows.append(
            f"{r.region_name:<22} {r.region_code:<5} {r.status:<8} "
            f"{r.n_pixels:>14,} {r.elapsed_s:>8.1f}  {path_str}"
        )
    return "\n".join([header, *rows])
