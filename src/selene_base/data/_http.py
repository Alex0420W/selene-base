"""Internal HTTP download helper shared by every dataset fetcher.

Two contracts:

- ``stream_to_file`` is idempotent: if ``dest`` already exists and meets
  ``min_bytes``, the function returns immediately. Otherwise it streams
  the URL to a ``.part`` file with a tqdm progress bar and atomically
  renames on success.
- All user-facing logging goes through ``typer.echo`` so the CLI and the
  inventory notebook see the same messages.
"""

from __future__ import annotations

from pathlib import Path

import requests
import typer
from tqdm import tqdm

DEFAULT_TIMEOUT = 60
DEFAULT_CHUNK = 1 << 20  # 1 MiB


def stream_to_file(
    url: str,
    dest: Path,
    *,
    min_bytes: int,
    label: str | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    chunk_size: int = DEFAULT_CHUNK,
) -> Path:
    """Download ``url`` to ``dest`` if not already present and large enough.

    Args:
        url: HTTP(S) URL to fetch.
        dest: Destination path. Parent directories are created.
        min_bytes: Minimum acceptable file size; the function asserts the
            resulting file is at least this large.
        label: Optional human-readable label for the progress bar; falls
            back to the destination filename.
        timeout: Per-request timeout in seconds.
        chunk_size: Streaming chunk size in bytes.

    Returns:
        ``dest`` (always — present file or freshly downloaded).

    Raises:
        AssertionError: If the downloaded file is smaller than ``min_bytes``.
        requests.HTTPError: If the server returns an error response.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        typer.echo(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")

    typer.echo(f"[fetch] {url}")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        total = int(response.headers.get("Content-Length", 0)) or None
        bar_label = label or dest.name
        with (
            part.open("wb") as fh,
            tqdm(
                total=total,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=bar_label,
            ) as bar,
        ):
            for chunk in response.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                fh.write(chunk)
                bar.update(len(chunk))

    size = part.stat().st_size
    if size < min_bytes:
        part.unlink(missing_ok=True)
        raise AssertionError(
            f"download too small: got {size:,} bytes, expected at least {min_bytes:,} (url={url})"
        )
    part.replace(dest)
    typer.echo(f"[done] {dest} ({size:,} bytes)")
    return dest
