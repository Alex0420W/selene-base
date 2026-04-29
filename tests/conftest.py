"""Pytest fixtures and global configuration.

Pin matplotlib to the non-interactive ``Agg`` backend before any test
imports a module that pulls in pyplot — otherwise our viz tests can
race a stale Tk root left over from earlier tests and trip a TclError.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
