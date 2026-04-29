"""Always-on smoke test: every dataset loader is importable and callable.

Catches accidental breakage from refactors even when no data has been
downloaded — runs in CI on every push.
"""

from __future__ import annotations

import pytest


@pytest.mark.parametrize(
    "name",
    [
        "load_raster",
        "load_lola_ldem",
        "load_diviner",
        "load_illumination",
        "load_lend",
        "load_crater_catalog",
    ],
)
def test_loader_callable(name: str) -> None:
    from selene_base.data import load

    fn = getattr(load, name)
    assert callable(fn), f"{name} is not callable"


@pytest.mark.parametrize(
    "name",
    [
        "download_robbins",
        "download_lola",
        "download_diviner",
        "download_illumination",
        "download_lend",
        "download_all",
    ],
)
def test_downloader_callable(name: str) -> None:
    from selene_base.data import download

    fn = getattr(download, name)
    assert callable(fn), f"{name} is not callable"


def test_dataset_registry_complete() -> None:
    from selene_base.cli import Dataset
    from selene_base.data.download import DATASETS

    cli_names = {d.value for d in Dataset if d.value != "all"}
    assert cli_names == set(DATASETS), (
        f"CLI Dataset enum and DATASETS registry diverged: "
        f"cli={cli_names!r} vs registry={set(DATASETS)!r}"
    )
