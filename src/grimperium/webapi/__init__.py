"""Read-only HTTP API behind the local dashboard.

Nothing here starts CREST/MOPAC or writes to a CSV or to the Semi-Imperium
store: every endpoint reads files that other workflows produce. FastAPI is
an optional dependency (extra ``dashboard``), so import :func:`create_app`
from :mod:`grimperium.webapi.app`.
"""

from grimperium.webapi.config import DashboardPaths

__all__ = ["DashboardPaths"]
