"""Earth line-of-sight criterion — rewards direct radio comms windows.

NASA's actual Artemis IV (formerly Artemis III) selection cares deeply about whether a site
can maintain direct radio contact with Earth without relying on a relay
satellite, for at least part of every operational period. The south
polar candidate regions cluster in the Earth-facing hemisphere
specifically because of this constraint: a site only has direct
line-of-sight to Earth when its local terrain horizon doesn't block
Earth's apparent position during favorable phases of the lunar
libration cycle.

This module derives a ``[0, 1]`` "Earth visibility" score from the
already-cached LOLA elevation grid in three steps:

1. :func:`derive_horizon_profile` ray-marches outward from each pixel
   in ``n_azimuths`` directions up to ``max_horizon_km``, tracking the
   maximum elevation angle of obstructing terrain (with curvature
   correction for the lunar sphere of radius 1737.4 km).
2. :func:`compute_earth_visibility_fraction` samples the libration
   ellipse (Earth's sub-Earth point cycles within roughly ±6.5° in
   latitude and ±7.9° in longitude over ~27 days) and counts the
   fraction of samples for which Earth's apparent elevation exceeds
   the horizon angle in the matching azimuth direction.
3. :func:`compute` maps that fraction to a ``[0, 1]`` score via a
   linear ramp anchored to operational comms thresholds.

The defaults — ``min_visibility = 0.20`` and ``target_visibility = 0.50``
— reflect operational realities: Apollo-era surface ops baselined a
``>20%`` direct-comms duty cycle as a crew safety floor; a sustained
habitat with redundant relay backup targets ``~50%`` for robust
operations. These values are physics-and-operations driven and were
chosen *before* the validation rerun, so they are not validation
chasing.

The libration sampling at ``n_libration_samples = 24`` is a coarse
parametric approximation of the full Lissajous trajectory traced by
the sub-Earth point — physical libration in latitude (period 27.55 d)
and optical libration in longitude (period 27.32 d) trace an ellipse
over a few months. 24 samples on the parametric ellipse is enough for
criterion ranking; documented as such.

All ray-marching and azimuth bookkeeping is done in the projected
**grid frame** (the ``+proj=stere`` polar stereographic CRS used
throughout the project), not the geographic frame. Earth's sub-Earth
direction is computed in the geographic frame and then rotated into
the grid frame at each pixel using the local grid convergence
``γ = atan2(x_p, y_p)`` of south polar stereographic. See the function
docstrings for the conventions.

Filled in week 9.
"""

from __future__ import annotations

import contextlib
import gc
import mmap as _mmap_mod

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from scipy.ndimage import map_coordinates

DEFAULT_MOON_RADIUS_M = 1_737_400.0
DEFAULT_N_AZIMUTHS = 36
DEFAULT_MAX_HORIZON_KM = 100.0
DEFAULT_N_LIBRATION_SAMPLES = 24
EARTH_MAX_LIBRATION_LAT_DEG = 6.5
EARTH_MAX_LIBRATION_LON_DEG = 7.9


def derive_horizon_profile(
    elevation: xr.DataArray,
    *,
    n_azimuths: int = DEFAULT_N_AZIMUTHS,
    max_horizon_km: float = DEFAULT_MAX_HORIZON_KM,
    pixel_size_m: float = 240.0,
    moon_radius_m: float = DEFAULT_MOON_RADIUS_M,
    n_distance_samples: int = 50,
    use_gpu: bool = False,
    azimuth_chunk_size: int | None = None,
    out: np.ndarray | None = None,
) -> xr.DataArray:
    """For each pixel, compute the maximum horizon elevation angle in
    ``n_azimuths`` grid-azimuth directions.

    The output is a 3D DataArray with dims ``(azimuth, y, x)`` in
    degrees. Negative angles mean the horizon dips below local
    horizontal (the lunar curvature drop dominates locally for a flat
    surface; expect ``-1.65°`` at 100 km on a perfectly flat surface).

    **Grid-azimuth convention.** Azimuth ``0`` is the ``+y_grid``
    direction (top of the COG / "screen north"); azimuth ``π/2`` is
    ``+x_grid``. Positive rotation is clockwise (east of grid north).
    For south polar stereographic centred at the south pole, the
    relationship to true geographic azimuth depends on the pixel's
    longitude — see :func:`compute_earth_visibility_fraction` for the
    grid-convergence correction.

    **Algorithm.**
    For each azimuth ``α``, the function ray-marches outward at
    ``n_distance_samples`` log-spaced distances from one pixel to
    ``max_horizon_km``, samples the elevation grid bilinearly via
    :func:`scipy.ndimage.map_coordinates`, and tracks the maximum
    apparent elevation angle observed along the ray:

    .. code-block:: text

        d              = distance along ray (m)
        curvature_drop = d^2 / (2 * R_moon)
        apparent_h     = elevation_at_(ray) - elevation_at_(viewer) - curvature_drop
        angle          = atan2(apparent_h, d)

    Log-spaced sampling balances near-field accuracy (where angular
    resolution per metre is highest) with reasonable runtime far out.

    **Runtime estimate.** On a 2533×2533 grid with the default 36
    azimuths and 50 distance samples, expect ~5–15 minutes on a
    developer laptop, dominated by the bilinear sampling pass.
    Single-threaded; vectorised within each (azimuth, distance) call.

    Args:
        elevation: DataArray of surface elevation in metres on a
            regular ``(y, x)`` grid in the projected polar stereographic
            CRS. NaN cells are treated as no-data and propagate.
        n_azimuths: Number of azimuthal sampling directions. 36 (every
            10°) is the default and is documented in the criterion's
            sanity tests.
        max_horizon_km: Maximum ray distance to consider. 100 km is the
            cutoff because beyond that the lunar curvature drop alone
            exceeds typical polar terrain relief and the rays' apparent
            elevation flattens.
        pixel_size_m: Grid resolution in metres (240 by default).
        moon_radius_m: Lunar radius in metres for the curvature-drop
            term (1737400 m).
        n_distance_samples: Number of log-spaced distances at which to
            sample each ray. 50 is enough that the sampled max is
            within a few hundredths of a degree of the continuous max
            for typical polar relief.
        use_gpu: If ``True``, run the inner ray-march on the GPU via
            CuPy + ``cupyx.scipy.ndimage.map_coordinates``. The output
            DataArray remains numpy-backed for downstream compatibility.
            Numerically equivalent to the CPU path within float32 ULP;
            see ``tests/test_los_to_earth_gpu.py`` for the tolerance.
        azimuth_chunk_size: When set to a positive integer < ``n_azimuths``,
            process azimuths in groups of this size and accumulate results
            into ``out`` slice by slice. Cuts peak GPU/unified-memory
            footprint of the per-chunk working buffer linearly with chunk
            size: at ``n_azimuths=36`` and ``chunk=9``, the working buffer
            is 1/4 the unchunked size. Required for resolutions ≤ 10 m on
            unified-memory GPUs (GB10 / DGX Spark) where the full
            ``(36, height, width)`` output would exceed the unified-memory
            cap. ``None`` (default) preserves the v1.5 behaviour: allocate
            the full output array on device. Numerically identical to the
            unchunked path within float32 precision (verified by
            ``tests/test_per_azimuth_chunking.py``).
        out: Optional pre-allocated host-side ``(n_azimuths, height, width)``
            float32 ndarray to write results into. Useful with
            ``np.lib.format.open_memmap`` to back the accumulator with an
            on-disk file rather than unified memory — required at 10 m on
            GB10 because host memory and GPU memory share the same physical
            pool, so any in-memory accumulator competes with the GPU
            allocation. When ``None`` (default), a regular ``np.full``
            host array is allocated; in chunked mode this still consumes
            unified memory and is unsuitable for large grids on
            unified-memory hardware.

    Returns:
        DataArray with dims ``(azimuth, y, x)`` and values in degrees.
        ``azimuth`` coordinate is the angle in degrees east of grid
        north; e.g. ``0, 10, 20, ..., 350`` for the default 36.

    Raises:
        ValueError: On out-of-range parameters or a non-2D ``elevation``.
    """
    if n_azimuths <= 0:
        raise ValueError(f"n_azimuths must be positive, got {n_azimuths!r}")
    if max_horizon_km <= 0:
        raise ValueError(f"max_horizon_km must be positive, got {max_horizon_km!r}")
    if pixel_size_m <= 0:
        raise ValueError(f"pixel_size_m must be positive, got {pixel_size_m!r}")
    if moon_radius_m <= 0:
        raise ValueError(f"moon_radius_m must be positive, got {moon_radius_m!r}")
    if n_distance_samples <= 0:
        raise ValueError(f"n_distance_samples must be positive, got {n_distance_samples!r}")
    if elevation.ndim != 2:
        raise ValueError(f"elevation must be 2D (y, x), got dims {elevation.dims!r}")
    if azimuth_chunk_size is not None and azimuth_chunk_size <= 0:
        raise ValueError(f"azimuth_chunk_size must be positive or None, got {azimuth_chunk_size!r}")

    if use_gpu:
        try:
            import cupy as xp
            from cupyx.scipy.ndimage import map_coordinates as _map_coords
        except ImportError as exc:
            raise RuntimeError(
                "use_gpu=True but cupy is not installed; "
                "`pip install cupy-cuda13x` (or cupy-cuda12x) on a CUDA host."
            ) from exc
    else:
        xp = np
        _map_coords = map_coordinates

    arr = xp.asarray(elevation.to_numpy(), dtype=xp.float32)
    height, width = arr.shape

    if out is not None:
        if out.shape != (n_azimuths, height, width):
            raise ValueError(
                f"out shape {out.shape!r} does not match expected ({n_azimuths}, {height}, {width})"
            )
        if out.dtype != np.float32:
            raise ValueError(f"out dtype must be float32, got {out.dtype!r}")
        host_out = out
        host_out[...] = -90.0
    else:
        host_out = np.full((n_azimuths, height, width), -90.0, dtype=np.float32)

    chunk = azimuth_chunk_size if azimuth_chunk_size is not None else n_azimuths
    chunk = min(chunk, n_azimuths)

    # distances_m only feeds python scalars into the loop body; keep on host
    # to avoid a device->host sync per (azimuth, distance) iteration.
    distances_m = np.geomspace(
        pixel_size_m, max_horizon_km * 1000.0, n_distance_samples, dtype=np.float32
    )

    row_idx = xp.arange(height, dtype=xp.float32)[:, None]
    col_idx = xp.arange(width, dtype=xp.float32)[None, :]

    for chunk_start in range(0, n_azimuths, chunk):
        chunk_end = min(chunk_start + chunk, n_azimuths)
        chunk_size = chunk_end - chunk_start

        # Per-chunk working buffer: (chunk_size, height, width) float32.
        # On unified memory this is the dominant allocation alongside the
        # scratch arrays inside the inner ray-march loop. At full
        # n_azimuths chunk this matches v1.5's behaviour exactly.
        horizon_chunk = xp.full((chunk_size, height, width), -90.0, dtype=xp.float32)

        for k_local, k in enumerate(range(chunk_start, chunk_end)):
            az_rad = 2.0 * np.pi * k / n_azimuths
            # Step direction in (col, row) per metre. +y_grid points to
            # smaller row indices in a standard COG (positive y resolution
            # ascending), so row_dir flips sign of cos(az).
            col_dir = float(np.sin(az_rad))
            row_dir = float(-np.cos(az_rad))

            horizon_az = horizon_chunk[k_local]
            for d_m in distances_m:
                d_pix = float(d_m) / pixel_size_m
                sample_row = row_idx + row_dir * d_pix
                sample_col = col_idx + col_dir * d_pix
                ele_at_d = _map_coords(
                    arr,
                    xp.stack(
                        [
                            xp.broadcast_to(sample_row, (height, width)),
                            xp.broadcast_to(sample_col, (height, width)),
                        ]
                    ),
                    order=1,
                    mode="constant",
                    # cupyx.scipy.ndimage.map_coordinates JIT only accepts
                    # ``np.nan`` (the module-level singleton) for cval=NaN; passing
                    # a plain ``float('nan')`` or ``np.float32('nan')`` triggers a
                    # codegen bug ("invalid type conversion: out = (Y)nan").
                    # Equivalent on CPU; required on GPU. cupy 14.0.1 / CUDA 13.
                    cval=np.nan,
                )
                curvature_drop = float(d_m) * float(d_m) / (2.0 * moon_radius_m)
                apparent_height = ele_at_d - arr - curvature_drop
                angle_deg = xp.degrees(xp.arctan2(apparent_height, float(d_m)))
                xp.fmax(horizon_az, angle_deg, out=horizon_az)

        # Copy this chunk's azimuths into the host accumulator. When
        # ``host_out`` is mmap-backed (GB10 high-resolution path), the OS
        # pages writes to disk lazily and the in-memory accumulator
        # footprint stays bounded by the per-chunk working buffer.
        if use_gpu:
            host_out[chunk_start:chunk_end] = horizon_chunk.get()
        else:
            host_out[chunk_start:chunk_end] = horizon_chunk

        # Release the per-chunk GPU buffer + free the CuPy mempool blocks
        # so the next chunk starts with a clean slate. Without this, the
        # mempool grows monotonically and the chunked path provides no
        # memory benefit over the unchunked one.
        del horizon_chunk
        if use_gpu:
            try:
                xp.get_default_memory_pool().free_all_blocks()
                xp.get_default_pinned_memory_pool().free_all_blocks()
            except Exception:
                pass

    if isinstance(host_out, np.memmap):
        host_out.flush()
    horizon_np = host_out

    azimuth_deg = np.degrees(2.0 * np.pi * np.arange(n_azimuths) / n_azimuths)
    out = xr.DataArray(
        horizon_np,
        dims=("azimuth", *elevation.dims),
        coords={
            "azimuth": azimuth_deg,
            **{d: elevation.coords[d] for d in elevation.dims if d in elevation.coords},
        },
        name="horizon_profile_deg",
    )
    if elevation.rio.crs is not None:
        out = out.rio.write_crs(elevation.rio.crs, inplace=False)
    return out


def compute_earth_visibility_fraction(
    horizon_profile: xr.DataArray,
    pixel_lat_deg: xr.DataArray,
    pixel_lon_deg: xr.DataArray,
    grid_convergence_rad: xr.DataArray,
    *,
    earth_max_libration_lat_deg: float = EARTH_MAX_LIBRATION_LAT_DEG,
    earth_max_libration_lon_deg: float = EARTH_MAX_LIBRATION_LON_DEG,
    n_libration_samples: int = DEFAULT_N_LIBRATION_SAMPLES,
    use_gpu: bool = False,
    row_chunk_size: int | None = None,
) -> xr.DataArray:
    """For each pixel, compute the fraction of the libration cycle
    during which Earth is above the local horizon.

    For each of ``n_libration_samples`` parametric points on the
    libration ellipse ``(φ_e, λ_e) = (A_lat sin θ, A_lon cos θ)``, the
    function:

    1. Computes Earth's apparent **geographic** elevation and azimuth
       at every pixel using spherical trigonometry on the lunar
       sphere.
    2. Rotates the geographic azimuth into the **grid frame** at each
       pixel via the local grid convergence
       ``γ = atan2(x_p, y_p)`` (south polar stereographic):
       ``α_grid = α_geo - γ``.
    3. Quantises the grid azimuth to the nearest of the
       ``len(horizon_profile.azimuth)`` buckets and compares Earth's
       elevation against the horizon angle for that bucket.

    The pixel's visibility fraction is the count of samples where
    Earth was above the horizon divided by ``n_libration_samples``.

    Args:
        horizon_profile: 3D DataArray ``(azimuth, y, x)`` from
            :func:`derive_horizon_profile` (degrees).
        pixel_lat_deg: 2D DataArray ``(y, x)`` of latitude in degrees,
            aligned with ``horizon_profile``'s ``(y, x)``.
        pixel_lon_deg: 2D DataArray ``(y, x)`` of longitude in degrees.
        grid_convergence_rad: 2D DataArray ``(y, x)`` of grid
            convergence ``γ`` in radians (positive east of grid north).
            For south polar stereographic this is
            ``atan2(x_projected, y_projected)``.
        earth_max_libration_lat_deg: Latitude semi-axis of the
            libration ellipse, default 6.5°.
        earth_max_libration_lon_deg: Longitude semi-axis, default 7.9°.
        n_libration_samples: Number of parametric points on the
            libration ellipse, default 24.
        use_gpu: If ``True``, run the libration sweep on the GPU via
            CuPy. Output is converted back to NumPy for downstream
            compatibility. Numerically equivalent to the CPU path
            within float64 ULP. Falls back to NumPy on hosts without
            CuPy.
        row_chunk_size: When set to a positive integer, iterate over
            row-slabs of that height instead of loading the full horizon
            array on device. Required at high resolutions on
            unified-memory GPUs (GB10): a 22000×22000 tile's full
            ``(36, 22000, 22000)`` float32 horizon profile is ~70 GB;
            row-chunked at ``row_chunk_size=2000`` the per-slab device
            footprint drops to ~6.3 GB. Pair with a ``horizon_profile``
            backed by ``np.lib.format.open_memmap`` (the v1.9 chunked
            preprocess cache) so each slab is paged from disk on demand
            rather than held in unified memory. Numerically identical to
            the unchunked path within float64 ULP.

    Returns:
        2D DataArray ``(y, x)`` of visibility fractions in ``[0, 1]``.
    """
    if row_chunk_size is not None and row_chunk_size <= 0:
        raise ValueError(f"row_chunk_size must be positive or None, got {row_chunk_size!r}")
    if n_libration_samples <= 0:
        raise ValueError(f"n_libration_samples must be positive, got {n_libration_samples!r}")
    if earth_max_libration_lat_deg <= 0 or earth_max_libration_lon_deg <= 0:
        raise ValueError(
            f"libration semi-axes must be positive, got "
            f"({earth_max_libration_lat_deg!r}, {earth_max_libration_lon_deg!r})"
        )

    if use_gpu:
        try:
            import cupy as xp
        except ImportError as exc:
            raise RuntimeError(
                "use_gpu=True but cupy is not installed; "
                "`pip install cupy-cuda13x` (or cupy-cuda12x) on a CUDA host."
            ) from exc
    else:
        xp = np

    horizon_np = horizon_profile.to_numpy()
    if horizon_np.ndim != 3:
        raise ValueError(f"horizon_profile must be 3D (azimuth, y, x), got {horizon_np.shape!r}")
    n_az = horizon_np.shape[0]

    lat_np = pixel_lat_deg.to_numpy()
    lon_np = pixel_lon_deg.to_numpy()
    gamma_np = grid_convergence_rad.to_numpy()
    if not (lat_np.shape == lon_np.shape == gamma_np.shape == horizon_np.shape[1:]):
        raise ValueError(
            f"shape mismatch: lat={lat_np.shape!r} lon={lon_np.shape!r} "
            f"gamma={gamma_np.shape!r} horizon-yx={horizon_np.shape[1:]!r}"
        )

    height, width = lat_np.shape
    chunk = row_chunk_size if row_chunk_size is not None else height
    chunk = min(chunk, height)

    # Libration parametric coordinates kept on host — they are length-24
    # 1-D vectors that feed scalars into the loop body, so a device
    # round-trip per iteration is wasted bandwidth.
    theta = 2.0 * np.pi * np.arange(n_libration_samples) / n_libration_samples
    libration_lat_rad = np.deg2rad(earth_max_libration_lat_deg) * np.sin(theta)
    libration_lon_rad = np.deg2rad(earth_max_libration_lon_deg) * np.cos(theta)

    az_step_rad = 2.0 * np.pi / n_az

    visibility_count = np.zeros((height, width), dtype=np.int32)

    for y0 in range(0, height, chunk):
        y1 = min(y0 + chunk, height)
        # Slice the row-slab from each input. A mmap-backed
        # ``horizon_profile`` pages just these rows from disk; in-memory
        # arrays just view them.
        horizon_slab = xp.asarray(horizon_np[:, y0:y1, :], dtype=xp.float32)
        lat_rad = xp.deg2rad(xp.asarray(lat_np[y0:y1, :], dtype=xp.float64))
        lon_rad = xp.deg2rad(xp.asarray(lon_np[y0:y1, :], dtype=xp.float64))
        gamma_rad = xp.asarray(gamma_np[y0:y1, :], dtype=xp.float64)

        cos_lat = xp.cos(lat_rad)
        sin_lat = xp.sin(lat_rad)
        slab_visibility = xp.zeros((y1 - y0, width), dtype=xp.int32)

        for phi_e, lambda_e in zip(libration_lat_rad, libration_lon_rad, strict=True):
            cos_phi_e = float(np.cos(phi_e))
            sin_phi_e = float(np.sin(phi_e))
            cos_lon_diff = xp.cos(lon_rad - float(lambda_e))
            sin_lon_diff = xp.sin(lon_rad - float(lambda_e))

            # Earth elevation: arcsin(e · u) where u is local up.
            e_dot_u = cos_lat * cos_phi_e * cos_lon_diff + sin_lat * sin_phi_e
            earth_elev_deg = xp.degrees(xp.arcsin(xp.clip(e_dot_u, -1.0, 1.0)))

            # Earth's *geographic* azimuth at the pixel: project Earth's 3D
            # direction onto the local north/east basis and atan2.
            earth_n = -sin_lat * cos_phi_e * cos_lon_diff + cos_lat * sin_phi_e
            earth_e = -cos_phi_e * sin_lon_diff
            earth_geo_az_rad = xp.arctan2(earth_e, earth_n)

            # Convert to grid azimuth via the pixel's grid convergence.
            # α_grid = (α_geo - γ) mod 2π, then quantise.
            earth_grid_az_rad = xp.mod(earth_geo_az_rad - gamma_rad, 2.0 * np.pi)
            k_idx = xp.mod(xp.round(earth_grid_az_rad / az_step_rad).astype(xp.int64), n_az)

            # Fancy-index the horizon at each pixel's azimuth bucket.
            horizon_at_az = xp.take_along_axis(horizon_slab, k_idx[None, :, :], axis=0)[0]
            slab_visibility += (earth_elev_deg > horizon_at_az).astype(xp.int32)

        if use_gpu:
            visibility_count[y0:y1, :] = slab_visibility.get()
        else:
            visibility_count[y0:y1, :] = slab_visibility

        # Free the slab's GPU footprint before paging in the next slab.
        del horizon_slab, lat_rad, lon_rad, gamma_rad, cos_lat, sin_lat, slab_visibility
        if use_gpu:
            try:
                xp.get_default_memory_pool().free_all_blocks()
                xp.get_default_pinned_memory_pool().free_all_blocks()
            except Exception:
                pass
        # Drop OS page-cache pages backing the mmap slab we just consumed.
        # Without this, the kernel keeps the strided 36-azimuth-region
        # reads in unified memory and steady-state pressure climbs slab
        # by slab until OOM. ``MADV_DONTNEED`` is a hint, not a guarantee
        # — but on Linux it's effective for backing-file mmap regions.
        if isinstance(horizon_np, np.memmap):
            with contextlib.suppress(AttributeError, OSError, ValueError):
                horizon_np._mmap.madvise(_mmap_mod.MADV_DONTNEED)
        gc.collect()

    visibility_fraction = visibility_count.astype(np.float64) / float(n_libration_samples)

    out = xr.DataArray(
        visibility_fraction,
        coords={
            d: pixel_lat_deg.coords[d] for d in pixel_lat_deg.dims if d in pixel_lat_deg.coords
        },
        dims=pixel_lat_deg.dims,
        name="los_visibility_fraction",
    )
    if pixel_lat_deg.rio.crs is not None:
        out = out.rio.write_crs(pixel_lat_deg.rio.crs, inplace=False)
    return out


def compute(
    visibility_fraction: xr.DataArray,
    *,
    min_visibility: float = 0.20,
    target_visibility: float = 0.50,
) -> xr.DataArray:
    """Map Earth-visibility fraction to a ``[0, 1]`` line-of-sight score.

    Linear ramp:

    .. code-block:: text

        visibility < min_visibility    -> 0.0
        visibility >= target_visibility -> 1.0
        between                         -> linear interpolation

    The defaults are physics-and-operations driven, not validation-tuned:

    - ``min_visibility = 0.20`` is the operational floor below which a
      site loses direct-Earth comms more than 80% of the libration
      cycle. Apollo-era surface ops baselined ``>20%`` direct comms as
      a crew-safety floor.
    - ``target_visibility = 0.50`` is the operational target at which
      direct-Earth comms is available a majority of the cycle, removing
      relay-satellite single-point-of-failure dependence for most
      activities.

    Args:
        visibility_fraction: DataArray of fractions in ``[0, 1]`` from
            :func:`compute_earth_visibility_fraction`. NaN values
            propagate.
        min_visibility: Lower threshold below which the score is 0.
        target_visibility: Upper threshold at and above which the
            score saturates at 1.

    Returns:
        DataArray of ``[0, 1]`` scores aligned with the input.

    Raises:
        ValueError: If thresholds are out of range or not strictly
            ordered.
    """
    if not (0 <= min_visibility < target_visibility <= 1):
        raise ValueError(
            f"need 0 <= min_visibility < target_visibility <= 1, "
            f"got min={min_visibility!r}, target={target_visibility!r}"
        )

    arr = visibility_fraction.to_numpy().astype(np.float64)
    score = (arr - min_visibility) / (target_visibility - min_visibility)
    score = np.clip(score, 0.0, 1.0)
    score = np.where(np.isnan(arr), np.nan, score)

    out = xr.DataArray(
        score,
        coords=visibility_fraction.coords,
        dims=visibility_fraction.dims,
        name="los_to_earth_score",
    )
    if visibility_fraction.rio.crs is not None:
        out = out.rio.write_crs(visibility_fraction.rio.crs, inplace=False)
    return out
