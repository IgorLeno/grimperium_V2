"""
Constants for GRIMPERIUM CLI.

Centralized paths and configuration values.
"""

from pathlib import Path


def _validate_project_root(candidate_path: Path) -> bool:
    """
    Validate that a path is the actual project root.

    Args:
        candidate_path: Path to validate as project root.

    Returns:
        True if path appears to be the project root, False otherwise.
    """
    # Check for expected directories/files that indicate project root
    expected_markers = [
        candidate_path / "data",
        candidate_path / "src" / "grimperium",
        candidate_path / "pyproject.toml",
    ]
    return any(marker.exists() for marker in expected_markers)


def _get_project_root_from_source() -> Path | None:
    """
    Compute project root from file location (development mode).

    Returns:
        Project root path if valid, None otherwise.
    """
    # Assuming structure: src/grimperium/cli/constants.py
    candidate = Path(__file__).resolve().parents[3]
    if _validate_project_root(candidate):
        return candidate
    return None


def _get_project_root_from_package() -> Path | None:
    """
    Locate project root using package resources (installed mode).

    Returns:
        Project root path if discoverable, None otherwise.
    """
    import logging

    logger = logging.getLogger(__name__)

    try:
        # Use importlib.resources (Python 3.9+)
        from importlib.resources import files

        package_path = files("grimperium")
        # If installed as editable, this might resolve to source
        if hasattr(package_path, "__fspath__"):
            candidate = Path(package_path.__fspath__()).resolve()
            # Walk up to find project root
            for parent in [candidate] + list(candidate.parents):
                if _validate_project_root(parent):
                    return parent
    except (ImportError, TypeError, AttributeError) as e:
        logger.debug("Package-based project root resolution failed: %s", e)
    return None


def get_project_root() -> Path:
    """
    Get validated project root path with fallback logic.

    This function attempts to locate the project root using multiple strategies:
    1. Compute from source file location (development mode)
    2. Use package resources API (installed mode)

    Returns:
        Validated project root path.

    Raises:
        RuntimeError: If project root cannot be located or validated.

    Example:
        >>> root = get_project_root()
        >>> data_dir = root / "data"
    """
    # Try source-based resolution first (most common in development)
    root = _get_project_root_from_source()
    if root is not None:
        return root

    # Try package-based resolution (installed package)
    root = _get_project_root_from_package()
    if root is not None:
        return root

    # All strategies failed
    raise RuntimeError(
        "Cannot locate GRIMPERIUM project root. "
        "Expected directory structure with 'data/', 'src/grimperium/', or 'pyproject.toml'. "
        "The CLI must be run from source (development mode) or properly installed as a package. "
        f"Current file location: {Path(__file__).resolve()}"
    )


# Project root with validation and fallback
_PROJECT_ROOT = get_project_root()

# Data paths for Phase A results
PHASE_A_RESULTS_DIR = _PROJECT_ROOT / "data" / "molecules_pm7" / "computed"
PHASE_A_RESULTS_FILE = PHASE_A_RESULTS_DIR / "phase_a_results.json"

# Dataset paths
DATA_DIR = _PROJECT_ROOT / "data"

# Available CBS Reference datasets
AVAILABLE_DATASETS = {
    "CBS Reference (CHON-only)": {
        "file": DATA_DIR / "thermo_cbs_chon.csv",
        "molecules": 29568,
        "last_updated": "2026-01-15",
        "status": "Ready",
        "description": "Filtered for C, H, O, N only (removed B, P, As, Ge)",
        "recommended": True,  # Default for Phase A
        "is_source": True,  # Source of truth for DatasetManager
    },
}

# Default dataset for batch processing
DEFAULT_DATASET = "CBS Reference (CHON-only)"
DEFAULT_DATASET_PATH = AVAILABLE_DATASETS[DEFAULT_DATASET]["file"]
