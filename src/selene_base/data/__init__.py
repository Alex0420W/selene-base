"""Data ingestion subpackage.

Provides three concerns kept deliberately separate:
- :mod:`selene_base.data.download` — fetch raw products from external archives.
- :mod:`selene_base.data.load` — read files from disk into typed in-memory objects.
- :mod:`selene_base.data.reproject` — warp rasters onto the common analysis grid.

Filled in across weeks 1 (download/load) and 2 (reproject).
"""
