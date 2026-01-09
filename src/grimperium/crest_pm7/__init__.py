"""CREST-PM7 Pipeline for computational chemistry.

This module provides a robust pipeline for conformer generation (CREST)
and PM7 optimization (MOPAC) of molecular structures.

Main Components:
- CRESTPM7Pipeline: Main orchestrator class
- PM7Config: Configuration dataclass
- PM7Result: Processing result dataclass
- validate_environment: Environment validation function

Example:
    >>> from grimperium.crest_pm7 import CRESTPM7Pipeline, PM7Config
    >>> config = PM7Config(phase="A", max_conformers=5)
    >>> pipeline = CRESTPM7Pipeline(config)
    >>> if pipeline.validate():
    ...     pipeline.setup()
    ...     result = pipeline.process_molecule("benzene", "c1ccccc1")
    ...     print(f"HOF: {result.most_stable_hof}")
"""

from .config import (
    AlertLevel,
    CRESTStatus,
    HOFConfidence,
    MOPACStatus,
    PM7Config,
    QualityGrade,
    TimeoutConfidence,
)
from .molecule_processor import ConformerData, MoleculeProcessor, PM7Result
from .pipeline import CRESTPM7Pipeline
from .result_evaluator import (
    MoleculeEvaluation,
    PhaseAEvaluation,
    ResultEvaluator,
    TOLERANCE_ABSOLUTE,
)
from .threshold_monitor import Alert, MonitoringMetrics, ThresholdMonitor
from .timeout_predictor import TimeoutPredictor
from .validation import validate_environment, ValidationResult

__all__ = [
    # Main classes
    "CRESTPM7Pipeline",
    "PM7Config",
    "PM7Result",
    # Enums
    "QualityGrade",
    "HOFConfidence",
    "TimeoutConfidence",
    "CRESTStatus",
    "MOPACStatus",
    "AlertLevel",
    # Data classes
    "ConformerData",
    "Alert",
    "MonitoringMetrics",
    "ValidationResult",
    "MoleculeEvaluation",
    "PhaseAEvaluation",
    # Classes
    "MoleculeProcessor",
    "ThresholdMonitor",
    "TimeoutPredictor",
    "ResultEvaluator",
    # Functions
    "validate_environment",
    # Constants
    "TOLERANCE_ABSOLUTE",
]

__version__ = "0.3.2"
