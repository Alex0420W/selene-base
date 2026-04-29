"""Smoke tests: every module imports and the CLI ``--help`` exits 0.

These tests catch the trivial-but-painful breakages: a typo in a stub
docstring breaking import, a missing dependency, a malformed typer
decorator. They do not exercise any pipeline logic.
"""

from __future__ import annotations

import importlib

import pytest
from typer.testing import CliRunner

MODULES = [
    "selene_base",
    "selene_base.cli",
    "selene_base.data",
    "selene_base.data.download",
    "selene_base.data.load",
    "selene_base.data.reproject",
    "selene_base.criteria",
    "selene_base.criteria.slope",
    "selene_base.criteria.illumination",
    "selene_base.criteria.thermal",
    "selene_base.criteria.ice",
    "selene_base.criteria.hazard",
    "selene_base.criteria.seismic",
    "selene_base.scoring",
    "selene_base.scoring.normalize",
    "selene_base.scoring.aggregate",
    "selene_base.scoring.ranking",
    "selene_base.viz",
    "selene_base.viz.webmap",
    "selene_base.viz.site_report",
]


@pytest.mark.parametrize("module_name", MODULES)
def test_module_imports(module_name: str) -> None:
    importlib.import_module(module_name)


def test_cli_root_help_exits_zero() -> None:
    from selene_base.cli import app

    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0, result.output
    assert "Usage" in result.output


@pytest.mark.parametrize("subcommand", ["download", "preprocess", "score", "rank", "viz"])
def test_cli_subcommand_help_exits_zero(subcommand: str) -> None:
    from selene_base.cli import app

    result = CliRunner().invoke(app, [subcommand, "--help"])
    assert result.exit_code == 0, result.output
    assert "Usage" in result.output
