"""CREST-PM7 Pipeline Configuration.

This module defines all configuration dataclasses and enums for the
CREST + PM7 computational chemistry pipeline.
"""

import tempfile
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


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
    ERROR = "ERROR"  # Generic error state for unexpected failures


class MoleculeStatus(str, Enum):
    """Status of molecule in the processing pipeline."""

    PENDING = "Pending"
    CURRENT_BATCH = "Current Batch"
    RUNNING = "Running"
    OK = "OK"
    RERUN = "Rerun"
    SKIP = "Skip"


class AlertLevel(str, Enum):
    """Alert level for threshold monitoring."""

    INFO = "INFO"
    WARNING = "WARNING"
    CRITICAL = "CRITICAL"


class Phase(str, Enum):
    """Processing phase for the pipeline."""

    A = "A"  # Initial validation phase
    B = "B"  # Extended testing phase
    C = "C"  # Calibration phase
    PRODUCTION = "PRODUCTION"  # Production phase


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
        crest_rmsd_threshold: RMSD threshold for CREST deduplication
        crest_opt_level: CREST optimization level (0-2)
        crest_optlev_label: CREST optimization level label
        crest_threads: Number of CREST threads
        crest_method: CREST method (gfn2, gfnff, gfn2//gfnff)
        crest_quick_mode: CREST quick mode (off/quick/squick/mquick)
        crest_v3: Use CREST v3 algorithm
        crest_nci: Enable CREST NCI mode
        crest_timeout: Timeout for CREST execution (seconds)
        mopac_timeout_base: Base timeout for MOPAC execution (seconds)
        mopac_timeout_margin: Margin multiplier for MOPAC timeout
        mopac_precise_scf: Enable MOPAC PRECISE SCF
        mopac_scf_threshold: MOPAC SCF convergence threshold
        xtb_preopt: Enable xTB pre-optimization
        nrotbonds_threshold_rigid_to_medium: Threshold for rigid -> medium
        nrotbonds_threshold_medium_to_flexible: Threshold for medium -> flexible
        timeout_predictor_recalibrate_interval: Molecules between recalibration
        success_rate_warning: Success rate threshold for warning
        success_rate_critical: Success rate threshold for critical alert
        monitor_window_size: Window size for monitoring metrics
        hof_extraction_min_samples: Minimum samples for HOF extraction rate check
        hof_extraction_threshold: Threshold for HOF extraction rate
        grade_degradation_min_samples: Minimum samples for grade degradation check
        grade_degradation_threshold: Threshold for grade degradation detection
        timeout_pattern_threshold: Threshold for timeout pattern detection
        scf_pattern_threshold: Threshold for SCF failure pattern detection
        consecutive_failures_warning: Consecutive failures for warning alert
        consecutive_failures_critical: Consecutive failures for critical alert
    """

    # Phase control
    phase: Phase = Phase.A

    # Executables
    crest_executable: str = "crest"
    mopac_executable: str = "mopac"

    # Paths - using cross-platform tempdir and absolute output path
    temp_dir: Path = field(
        default_factory=lambda: Path(tempfile.gettempdir()) / "crest_pm7"
    )
    output_dir: Path = field(
        default_factory=lambda: Path.cwd() / "data/molecules_pm7/computed"
    )

    # CREST settings
    max_conformers: int = 10
    energy_window: float = 6.0  # kcal/mol
    crest_rmsd_threshold: float = 0.125  # Angstrom
    crest_opt_level: int = 2  # 0=loose, 1=normal, 2=tight
    crest_optlev_label: str = "tight"
    crest_threads: int = 4
    crest_method: str = "gfn2"
    crest_quick_mode: str = "off"
    crest_v3: bool = True
    crest_nci: bool = False
    crest_timeout: float = 300.0  # seconds

    # MOPAC settings
    mopac_timeout_base: float = 120.0  # seconds
    mopac_timeout_margin: float = 1.3
    mopac_precise_scf: bool = True
    mopac_scf_threshold: float = 1.0e-4

    # xTB pre-optimization
    xtb_preopt: bool = True

    # Conformer selection thresholds (calibrated in Phase C)
    nrotbonds_threshold_rigid_to_medium: int = 1
    nrotbonds_threshold_medium_to_flexible: int = 4

    # Timeout predictor
    timeout_predictor_recalibrate_interval: int = 50

    # Monitoring thresholds
    success_rate_warning: float = 0.80
    success_rate_critical: float = 0.70
    monitor_window_size: int = 50

    # HOF extraction thresholds (for ThresholdMonitor)
    hof_extraction_min_samples: int = 10
    hof_extraction_threshold: float = 0.8

    # Grade degradation thresholds
    grade_degradation_min_samples: int = 10
    grade_degradation_threshold: float = 0.5

    # Timeout/SCF pattern thresholds
    timeout_pattern_threshold: float = 0.2
    scf_pattern_threshold: float = 0.15

    # Consecutive failures threshold
    consecutive_failures_warning: int = 3
    consecutive_failures_critical: int = 5

    def __post_init__(self) -> None:
        """Convert string paths to Path objects and validate phase."""
        if isinstance(self.temp_dir, str):
            self.temp_dir = Path(self.temp_dir)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

        # Validate and convert phase if it is a string
        if isinstance(self.phase, str):
            valid_phases = {p.value for p in Phase}
            if self.phase not in valid_phases:
                raise ValueError(
                    f"Invalid phase '{self.phase}'. "
                    f"Must be one of: {', '.join(sorted(valid_phases))}"
                )
            self.phase = Phase(self.phase)

        # Normalize and validate CREST settings
        if isinstance(self.crest_method, str):
            self.crest_method = self.crest_method.strip().lower()
        if isinstance(self.crest_quick_mode, str):
            self.crest_quick_mode = self.crest_quick_mode.strip().lower()

        valid_crest_methods = {"gfn2", "gfnff", "gfn2//gfnff"}
        if self.crest_method not in valid_crest_methods:
            raise ValueError(
                f"Invalid crest_method '{self.crest_method}'. "
                f"Must be one of: {', '.join(sorted(valid_crest_methods))}"
            )

        valid_quick_modes = {"off", "quick", "squick", "mquick"}
        if self.crest_quick_mode not in valid_quick_modes:
            raise ValueError(
                f"Invalid crest_quick_mode '{self.crest_quick_mode}'. "
                f"Must be one of: {', '.join(sorted(valid_quick_modes))}"
            )

        # Keep crest_optlev_label derived from crest_opt_level (canonical numeric field)
        if isinstance(self.crest_opt_level, str):
            try:
                self.crest_opt_level = int(self.crest_opt_level)
            except ValueError as e:  # pragma: no cover - defensive
                raise ValueError(
                    f"crest_opt_level must be an integer 0..2, got '{self.crest_opt_level}'"
                ) from e

        optlev_label_map = {0: "loose", 1: "normal", 2: "tight"}
        if self.crest_opt_level not in optlev_label_map:
            raise ValueError(
                f"Invalid crest_opt_level '{self.crest_opt_level}'. "
                f"Must be one of: {', '.join(str(k) for k in sorted(optlev_label_map))}"
            )

        expected_label = optlev_label_map[self.crest_opt_level]
        normalized_label = str(self.crest_optlev_label).strip().lower()
        if normalized_label and normalized_label != expected_label:
            raise ValueError(
                "Inconsistent CREST optlev settings: "
                f"crest_opt_level={self.crest_opt_level} implies "
                f"crest_optlev_label='{expected_label}', got '{normalized_label}'"
            )
        self.crest_optlev_label = expected_label

        # Validate threshold fields
        if not (0 <= self.hof_extraction_threshold <= 1):
            raise ValueError(
                f"hof_extraction_threshold must be between 0 and 1, got {self.hof_extraction_threshold}"
            )
        if not (0 <= self.grade_degradation_threshold <= 1):
            raise ValueError(
                f"grade_degradation_threshold must be between 0 and 1, got {self.grade_degradation_threshold}"
            )
        if not (0 <= self.timeout_pattern_threshold <= 1):
            raise ValueError(
                f"timeout_pattern_threshold must be between 0 and 1, got {self.timeout_pattern_threshold}"
            )
        if not (0 <= self.scf_pattern_threshold <= 1):
            raise ValueError(
                f"scf_pattern_threshold must be between 0 and 1, got {self.scf_pattern_threshold}"
            )
        if not (0 <= self.success_rate_warning <= 1):
            raise ValueError(
                f"success_rate_warning must be between 0 and 1, got {self.success_rate_warning}"
            )
        if not (0 <= self.success_rate_critical <= 1):
            raise ValueError(
                f"success_rate_critical must be between 0 and 1, got {self.success_rate_critical}"
            )
        if self.success_rate_warning <= self.success_rate_critical:
            raise ValueError(
                f"success_rate_warning ({self.success_rate_warning}) must be greater than "
                f"success_rate_critical ({self.success_rate_critical})"
            )

        # Validate min_samples fields
        if self.hof_extraction_min_samples <= 0:
            raise ValueError(
                f"hof_extraction_min_samples must be positive, got {self.hof_extraction_min_samples}"
            )
        if self.grade_degradation_min_samples <= 0:
            raise ValueError(
                f"grade_degradation_min_samples must be positive, got {self.grade_degradation_min_samples}"
            )

        # Validate consecutive failures thresholds
        if self.consecutive_failures_warning <= 0:
            raise ValueError(
                f"consecutive_failures_warning must be positive, got {self.consecutive_failures_warning}"
            )
        if self.consecutive_failures_critical <= 0:
            raise ValueError(
                f"consecutive_failures_critical must be positive, got {self.consecutive_failures_critical}"
            )
        if self.consecutive_failures_warning >= self.consecutive_failures_critical:
            raise ValueError(
                f"consecutive_failures_warning ({self.consecutive_failures_warning}) must be less than "
                f"consecutive_failures_critical ({self.consecutive_failures_critical})"
            )

        # Validate CREST/MOPAC settings
        if self.crest_opt_level not in (0, 1, 2):
            raise ValueError(
                f"crest_opt_level must be 0, 1, or 2, got {self.crest_opt_level}"
            )
        if self.crest_rmsd_threshold <= 0:
            raise ValueError(
                f"crest_rmsd_threshold must be > 0, got {self.crest_rmsd_threshold}"
            )
        if self.crest_threads <= 0:
            raise ValueError(f"crest_threads must be > 0, got {self.crest_threads}")
        if self.mopac_scf_threshold <= 0:
            raise ValueError(
                f"mopac_scf_threshold must be > 0, got {self.mopac_scf_threshold}"
            )

    def ensure_directories(self) -> None:
        """Create necessary directories if they don't exist."""
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "logs").mkdir(parents=True, exist_ok=True)
