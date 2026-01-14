"""
Constants for GRIMPERIUM CLI.

Centralized paths and configuration values.
"""

from pathlib import Path

# Project root (assuming structure: src/grimperium/cli/constants.py)
_PROJECT_ROOT = Path(__file__).resolve().parents[3]

# Data paths for Phase A results
PHASE_A_RESULTS_DIR = _PROJECT_ROOT / "data" / "molecules_pm7" / "computed"
PHASE_A_RESULTS_FILE = PHASE_A_RESULTS_DIR / "phase_a_results.json"
