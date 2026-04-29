"""Pipeline orchestration glue.

CLI subcommands (``selene preprocess``, ``selene score``) call the
``run`` function in the matching submodule. Keeping the heavy lifting
out of :mod:`selene_base.cli` lets us test the pipeline end-to-end on
synthetic input without going through typer.
"""
