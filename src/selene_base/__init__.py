"""selene-base: multi-criteria habitat suitability for the lunar south pole.

A ranked-site analyzer that fuses LRO topography, illumination, thermal,
hydrogen, crater-catalog, and lobate-scarp data with an Apollo-derived
shallow-moonquake context to score Artemis candidate base sites.

The public surface is the :mod:`selene_base.cli` typer application; library
consumers can also import the criteria and scoring modules directly.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["__version__"]
