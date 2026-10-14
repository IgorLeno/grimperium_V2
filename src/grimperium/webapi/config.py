"""File locations read by the dashboard API."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# src/grimperium/webapi/config.py -> repository root. The API reads the
# working-tree data, so it is meant to run from a source checkout.
_REPO_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class DashboardPaths:
    """Every input file of the API; tests point these at fixtures."""

    pm7_csv: Path
    semi_imperium_store: Path
    exp_csv: Path
    exp_vs_cbs_csv: Path
    cbs_csv: Path | None = None
    """CHON CBS reference table; ``None`` reports it as not configured."""

    @classmethod
    def default(cls, root: Path | None = None) -> DashboardPaths:
        """Paths of the active datasets under ``root/data``.

        ``GRIMPERIUM_ROOT`` points the API at another checkout-shaped folder
        (used to rehearse the routine on a copy of the data).
        """
        base = root or Path(os.environ.get("GRIMPERIUM_ROOT", _REPO_ROOT))
        data = base / "data"
        return cls(
            pm7_csv=data / "pm7_fila_pipeline.csv",
            semi_imperium_store=data / "semi_imperium",
            exp_csv=data / "experimental.csv",
            exp_vs_cbs_csv=data / "experimental_vs_cbs.csv",
            cbs_csv=data / "cbs_referencia.csv",
        )
