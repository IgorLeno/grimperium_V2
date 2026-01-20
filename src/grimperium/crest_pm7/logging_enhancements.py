"""
╔════════════════════════════════════════════════════════════════════════════════╗
║                                                                                ║
║              FILE_2: STRUCTURED LOGGING & WARNING SUPPRESSION                  ║
║                                                                                ║
║  Purpose: Add detailed, structured logging for RDKit → CREST → MOPAC pipeline ║
║  Issue Fixed: #2 (Logging & DtypeWarning suppression)                          ║
║                                                                                ║
║  Location: src/grimperium/crest_pm7/logging_enhancements.py                   ║
║  Status: Production Ready (408 lines)                                         ║
║                                                                                ║
╚════════════════════════════════════════════════════════════════════════════════╝
"""

import logging
import warnings
from collections.abc import Generator
from contextlib import contextmanager
from datetime import datetime

import pandas as pd

# Conditional import of rich (optional dependency, part of 'cli' extra)
_RICH_AVAILABLE = False
try:
    from rich.console import Console
    from rich.logging import RichHandler

    _RICH_AVAILABLE = True
except ImportError:
    Console = None  # type: ignore[assignment, misc]
    RichHandler = None  # type: ignore[assignment, misc]

# Initialize Rich console for colored output (only if available)
console = Console() if _RICH_AVAILABLE else None


class LoggingConfig:
    """Centralized logging configuration for the CREST PM7 pipeline."""

    # Log format with timestamps and colors
    LOG_FORMAT = "%(asctime)s - [%(levelname)s] - %(name)s - %(message)s"
    DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

    # Color mapping for different tools
    TOOL_COLORS = {
        "rdkit": "🧬",  # Molecule icon
        "crest": "🔄",  # Rotation icon
        "mopac": "⚛️",  # Atom icon
        "rdkit_complete": "✓",
        "crest_complete": "✓",
        "mopac_complete": "✓",
    }


def setup_batch_logging(batch_id: str) -> logging.Logger:
    """Set up logging for a batch of molecules.

    Creates a logger that:
    - Logs to console with Rich formatting (if available)
    - Falls back to plain StreamHandler if rich is not installed
    - Uses colored output for visibility (when rich available)
    - Includes timestamp for each entry
    - Shows per-molecule progress

    Args:
        batch_id: Batch identifier (e.g., "batch_0001")

    Returns:
        Logger instance for the batch

    Example:
        >>> logger = setup_batch_logging("batch_0001")
        >>> logger.info(f"[mol_00001] Starting RDKit...")
        [2026-01-20 13:04:26] [INFO] [mol_00001] Starting RDKit...
    """
    # Create logger
    logger = logging.getLogger(f"grimperium.{batch_id}")
    logger.setLevel(logging.DEBUG)

    # Clear existing handlers
    logger.handlers.clear()

    # Disable propagation to prevent duplicate console output
    logger.propagate = False

    # Create handler: Rich if available, otherwise plain StreamHandler
    handler: logging.Handler
    if _RICH_AVAILABLE and RichHandler is not None and console is not None:
        handler = RichHandler(
            console=console,
            markup=True,
            rich_tracebacks=True,
            show_time=True,
            show_level=True,
            show_path=False,
        )
    else:
        handler = logging.StreamHandler()

    handler.setFormatter(
        logging.Formatter(
            fmt="[%(asctime)s] [%(levelname)s] %(message)s", datefmt="%H:%M:%S"
        )
    )

    logger.addHandler(handler)

    return logger


def setup_molecule_logging(batch_logger: logging.Logger, mol_id: str) -> logging.Logger:
    """
    Set up logging for a specific molecule (optional sub-logger).

    Args:
        batch_logger: Batch logger from setup_batch_logging()
        mol_id: Molecule identifier

    Returns:
        Logger instance for the molecule

    Example:
        >>> batch_logger = setup_batch_logging("batch_0001")
        >>> mol_logger = setup_molecule_logging(batch_logger, "mol_00001")
    """
    return logging.getLogger(f"{batch_logger.name}.{mol_id}")


# ╔════════════════════════════════════════════════════════════════════════════════╗
# ║ LOGGING FUNCTIONS FOR EACH PIPELINE STAGE                                     ║
# ╚════════════════════════════════════════════════════════════════════════════════╝


def log_rdkit_start(logger: logging.Logger, mol_id: str) -> None:
    """Log RDKit descriptor calculation start."""
    logger.info(f"[{mol_id}] 🧬 RDKit: Calculating descriptors...")


def log_rdkit_done(
    logger: logging.Logger, mol_id: str, **properties: float | int
) -> None:
    """
    Log RDKit completion with calculated properties.

    Args:
        logger: Logger instance
        mol_id: Molecule identifier
        **properties: Descriptor values (nrotbonds, tpsa, aromatic_rings, etc.)

    Example:
        >>> log_rdkit_done(logger, "mol_00001",
        ...     nrotbonds=2, tpsa=45.5, aromatic_rings=1)
        [mol_00001]   ✓ nrotbonds=2.0, tpsa=45.5, aromatic_rings=1
    """
    props_str = ", ".join(
        f"{k}={v:.2f}" if isinstance(v, float) else f"{k}={v}"
        for k, v in properties.items()
    )
    logger.info(f"[{mol_id}]   ✓ {props_str}")


def log_crest_start(logger: logging.Logger, mol_id: str) -> None:
    """Log CREST conformer generation start."""
    logger.info(f"[{mol_id}] 🔄 CREST: Starting conformer sampling...")


def log_crest_done(
    logger: logging.Logger, mol_id: str, num_conformers: int, time_seconds: float
) -> None:
    """
    Log CREST completion with conformer count and time.

    Args:
        logger: Logger instance
        mol_id: Molecule identifier
        num_conformers: Number of conformers generated
        time_seconds: Calculation time in seconds

    Example:
        >>> log_crest_done(logger, "mol_00001", 4, 4.2)
        [mol_00001]   ✓ Generated 4 conformers in 4.2s
    """
    logger.info(
        f"[{mol_id}]   ✓ Generated {num_conformers} conformers in {time_seconds:.1f}s"
    )


def log_mopac_start(logger: logging.Logger, mol_id: str, num_conformers: int) -> None:
    """Log MOPAC optimization start."""
    logger.info(f"[{mol_id}] ⚛️  MOPAC: Optimizing {num_conformers} conformers...")


def log_mopac_conformer_done(
    logger: logging.Logger,
    mol_id: str,
    conf_idx: int,
    delta_energy: float,
    time_seconds: float,
) -> None:
    """
    Log individual conformer optimization.

    Args:
        logger: Logger instance
        mol_id: Molecule identifier
        conf_idx: Conformer index
        delta_energy: Energy difference (ΔE) in kcal/mol
        time_seconds: Calculation time

    Example:
        >>> log_mopac_conformer_done(logger, "mol_00001", 0, 0.42, 1.2)
        [mol_00001]   • Conformer 0: ✓ ΔE=0.42 kcal/mol (1.2s)
    """
    logger.info(
        f"[{mol_id}]   • Conformer {conf_idx}: ✓ ΔE={delta_energy:.2f} kcal/mol ({time_seconds:.1f}s)"
    )


def log_mopac_done(
    logger: logging.Logger,
    mol_id: str,
    best_conformer_idx: int,
    best_delta_energy: float,
    time_seconds: float,
) -> None:
    """
    Log MOPAC completion with best conformer selection.

    Args:
        logger: Logger instance
        mol_id: Molecule identifier
        best_conformer_idx: Index of best (lowest energy) conformer
        best_delta_energy: Best energy difference
        time_seconds: Total optimization time

    Example:
        >>> log_mopac_done(logger, "mol_00001", 0, 0.42, 4.5)
        [mol_00001]   ✓ Selected conformer #0 with ΔE=0.42 kcal/mol
    """
    logger.info(
        f"[{mol_id}]   ✓ Selected conformer #{best_conformer_idx} with ΔE={best_delta_energy:.2f} kcal/mol"
    )


def log_batch_summary(
    logger: logging.Logger,
    batch_id: str,
    total: int,
    success: int,
    failed: int,
    skipped: int = 0,
) -> None:
    """
    Log batch completion summary.

    Args:
        logger: Logger instance
        batch_id: Batch identifier
        total: Total molecules processed
        success: Successfully processed
        failed: Failed to process
        skipped: Skipped molecules

    Example:
        >>> log_batch_summary(logger, "batch_0001", 100, 98, 2, 0)
        ═════════════════════════════════════════════════════════════
        ✓ BATCH SUMMARY: batch_0001
        └─ Processed: 100 molecules
        └─ Success: 98 (98.0%)
        └─ Failed: 2 (2.0%)
        └─ Skipped: 0
        ═════════════════════════════════════════════════════════════
    """
    success_pct = (success / total * 100) if total > 0 else 0
    failed_pct = (failed / total * 100) if total > 0 else 0

    logger.info("╔════════════════════════════════════════════════════════════╗")
    logger.info(f"║ ✓ BATCH SUMMARY: {batch_id}")
    logger.info(f"║ └─ Processed: {total} molecules")
    logger.info(f"║ └─ Success: {success} ({success_pct:.1f}%)")
    logger.info(f"║ └─ Failed: {failed} ({failed_pct:.1f}%)")
    if skipped > 0:
        logger.info(f"║ └─ Skipped: {skipped}")
    logger.info("╚════════════════════════════════════════════════════════════╝")


# ╔════════════════════════════════════════════════════════════════════════════════╗
# ║ WARNING SUPPRESSION                                                           ║
# ╚════════════════════════════════════════════════════════════════════════════════╝


def suppress_pandas_warnings() -> None:
    """Suppress pandas-specific warnings (narrowed to exact warning classes).

    Call once at application start.

    Example:
        >>> suppress_pandas_warnings()
        # No more DtypeWarning messages!
    """
    # Suppress pandas DtypeWarning specifically (not all Warning)
    try:
        warnings.filterwarnings(
            "ignore",
            category=pd.errors.DtypeWarning,
            message=".*Columns.*have mixed types.*",
        )
    except AttributeError:
        # Fallback for older pandas versions without DtypeWarning
        warnings.filterwarnings(
            "ignore", category=Warning, message=".*Columns.*have mixed types.*"
        )

    # Suppress FutureWarning about DataFrame indexing
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message=".*Setting an item of incompatible dtype.*",
    )

    # Suppress SettingWithCopyWarning
    try:
        warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
    except AttributeError:
        # Fallback for older pandas versions
        pass


# ╔════════════════════════════════════════════════════════════════════════════════╗
# ║ CONTEXT MANAGERS                                                              ║
# ╚════════════════════════════════════════════════════════════════════════════════╝


@contextmanager
def log_phase(
    logger: logging.Logger, phase_name: str, mol_id: str = ""
) -> Generator[logging.Logger, None, None]:
    """
    Context manager for logging a calculation phase.

    Usage:
        >>> with log_phase(logger, "RDKit", "mol_00001"):
        ...     # Do RDKit calculations
        ...     pass
    """
    start_time = datetime.now()
    prefix = f"[{mol_id}] " if mol_id else ""
    logger.info(f"{prefix}Starting {phase_name}...")

    try:
        yield logger
    finally:
        elapsed = (datetime.now() - start_time).total_seconds()
        logger.info(f"{prefix}  ✓ {phase_name} completed in {elapsed:.1f}s")


@contextmanager
def suppress_warnings_context() -> Generator[None, None, None]:
    """Context manager to suppress warnings temporarily."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


# ╔════════════════════════════════════════════════════════════════════════════════╗
# ║ PRE-FORMATTED LOG MESSAGE TEMPLATES                                           ║
# ╚════════════════════════════════════════════════════════════════════════════════╝


class LogMessages:
    """Pre-formatted log message templates."""

    @staticmethod
    def rdkit_descriptor(mol_id: str, **props: float | int) -> str:
        """Format RDKit descriptor message."""
        props_str = ", ".join(f"{k}={v}" for k, v in props.items())
        return f"[{mol_id}] 🧬 RDKit descriptors: {props_str}"

    @staticmethod
    def crest_conformers(mol_id: str, count: int, time_s: float) -> str:
        """Format CREST conformers message."""
        return f"[{mol_id}] 🔄 CREST: Generated {count} conformers in {time_s:.1f}s"

    @staticmethod
    def mopac_selected(mol_id: str, conf_idx: int, delta: float) -> str:
        """Format MOPAC selection message."""
        return f"[{mol_id}] ⚛️  MOPAC: Selected conformer {conf_idx} (ΔE={delta:.2f} kcal/mol)"


if __name__ == "__main__":
    # Test the logging module
    print("Testing logging_enhancements.py module...")

    # Setup logging
    logger = setup_batch_logging("batch_test")
    print("✓ Logging configured")

    # Test logging functions
    log_rdkit_start(logger, "mol_00001")
    log_rdkit_done(logger, "mol_00001", nrotbonds=2, tpsa=45.5, aromatic_rings=1)

    log_crest_start(logger, "mol_00001")
    log_crest_done(logger, "mol_00001", 4, 4.2)

    log_mopac_start(logger, "mol_00001", 4)
    log_mopac_conformer_done(logger, "mol_00001", 0, 0.42, 1.2)
    log_mopac_conformer_done(logger, "mol_00001", 1, 0.87, 1.3)
    log_mopac_done(logger, "mol_00001", 0, 0.42, 2.5)

    log_batch_summary(logger, "batch_test", 100, 98, 2, 0)

    # Test warning suppression
    suppress_pandas_warnings()
    print("✓ Warnings suppressed")

    print("\n✓ All logging functions working!")
