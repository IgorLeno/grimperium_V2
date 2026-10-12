"""File locations read by the dashboard API."""

from __future__ import annotations

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
    def default(cls, root: Path = _REPO_ROOT) -> DashboardPaths:
        """Paths of the active datasets under ``root/data``."""
        data = root / "data"
        return cls(
            pm7_csv=data / "thermo_pm7.csv",
            semi_imperium_store=data / "semi_imperium",
            exp_csv=data / "thermo_exp_hf.csv",
            exp_vs_cbs_csv=data / "thermo_exp_hf_vs_cbs.csv",
            cbs_csv=data / "thermo_cbs_chon.csv",
        )
