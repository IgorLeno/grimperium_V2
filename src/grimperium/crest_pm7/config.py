"""CREST-PM7 Pipeline Configuration.

This module defines all configuration dataclasses and enums for the
CREST + PM7 computational chemistry pipeline.
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


class QualityGrade(str, Enum):
    """Quality grade for molecule processing results."""

    A = "A"  # Reliable - high confidence
    B = "B"  # Acceptable - minor issues
    C = "C"  # Suspect - significant issues
    FAILED = "FAILED"  # Critical error


class HOFConfidence(str, Enum):
    """Confidence level for Heat of Formation extraction."""

    HIGH = "HIGH"  # Patterns 1-2: FINAL/HEAT OF FORMATION
    MEDIUM = "MEDIUM"  # Patterns 3-4: compact/thermal
    LOW = "LOW"  # Pattern 5: fallback


class TimeoutConfidence(str, Enum):
    """Confidence level for timeout prediction."""

    LOW = "LOW"  # < 20 samples (heuristic)
    MEDIUM = "MEDIUM"  # 20-50 samples (model + 50% margin)
    HIGH = "HIGH"  # 50+ samples (model + 30% margin)


class CRESTStatus(str, Enum):
    """Status of CREST conformer generation."""

    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    NOT_ATTEMPTED = "NOT_ATTEMPTED"


class MOPACStatus(str, Enum):
    """Status of MOPAC PM7 optimization."""

    SUCCESS = "SUCCESS"
    TIMEOUT = "TIMEOUT"
    SCF_FAILED = "SCF_FAILED"
    GEOMETRY_ERROR = "GEOMETRY_ERROR"
    NOT_ATTEMPTED = "NOT_ATTEMPTED"


class AlertLevel(str, Enum):
    """Alert level for threshold monitoring."""

    INFO = "INFO"
    WARNING = "WARNING"
    CRITICAL = "CRITICAL"


@dataclass
class PM7Config:
    """Configuration for the CREST-PM7 pipeline.

    Attributes:
        phase: Current processing phase (A, B, C, or production)
        crest_executable: Path to CREST executable
        mopac_executable: Path to MOPAC executable
        temp_dir: Directory for temporary files
        output_dir: Directory for output files
        max_conformers: Maximum conformers to generate/process
        energy_window: Energy window for conformer selection (kcal/mol)
        crest_timeout: Timeout for CREST execution (seconds)
        mopac_timeout_base: Base timeout for MOPAC execution (seconds)
        mopac_timeout_margin: Margin multiplier for MOPAC timeout
        nrotbonds_threshold_rigid_to_medium: Threshold for rigid -> medium
        nrotbonds_threshold_medium_to_flexible: Threshold for medium -> flexible
        timeout_predictor_recalibrate_interval: Molecules between recalibration
        success_rate_warning: Success rate threshold for warning
        success_rate_critical: Success rate threshold for critical alert
        monitor_window_size: Window size for monitoring metrics
    """

    # Phase control
    phase: str = "A"

    # Executables
    crest_executable: str = "crest"
    mopac_executable: str = "mopac"

    # Paths
    temp_dir: Path = field(default_factory=lambda: Path("/tmp/crest_pm7"))
    output_dir: Path = field(default_factory=lambda: Path("data/molecules_pm7/computed"))

    # CREST settings
    max_conformers: int = 10
    energy_window: float = 6.0  # kcal/mol
    crest_timeout: float = 300.0  # seconds

    # MOPAC settings
    mopac_timeout_base: float = 120.0  # seconds
    mopac_timeout_margin: float = 1.3

    # Conformer selection thresholds (calibrated in Phase C)
    nrotbonds_threshold_rigid_to_medium: int = 1
    nrotbonds_threshold_medium_to_flexible: int = 4

    # Timeout predictor
    timeout_predictor_recalibrate_interval: int = 50

    # Monitoring thresholds
    success_rate_warning: float = 0.80
    success_rate_critical: float = 0.70
    monitor_window_size: int = 50

    def __post_init__(self) -> None:
        """Convert string paths to Path objects."""
        if isinstance(self.temp_dir, str):
            self.temp_dir = Path(self.temp_dir)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

    def ensure_directories(self) -> None:
        """Create necessary directories if they don't exist."""
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "logs").mkdir(parents=True, exist_ok=True)
