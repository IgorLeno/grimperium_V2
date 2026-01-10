"""Main CREST-PM7 Pipeline orchestrator.

Coordinates the full processing pipeline for multiple molecules.
"""

import json
import logging
from collections.abc import Iterator
from enum import Enum
from pathlib import Path
from typing import Callable, Optional

from .config import MOPACStatus, PM7Config
from .logging_utils import log_molecule_complete, log_molecule_start, setup_logging
from .molecule_processor import MoleculeProcessor, PM7Result
from .result_evaluator import PhaseAEvaluation, ResultEvaluator
from .threshold_monitor import Alert, ThresholdMonitor
from .timeout_predictor import TimeoutPredictor
from .validation import validate_environment

LOG = logging.getLogger("grimperium.crest_pm7.pipeline")


class CRESTPM7Pipeline:
    """Main pipeline orchestrator for CREST + PM7 processing.

    Coordinates:
    - Environment validation
    - Molecule processing
    - Quality monitoring
    - Result evaluation
    - Logging
    """

    def __init__(self, config: Optional[PM7Config] = None) -> None:
        """Initialize pipeline.

        Args:
            config: Pipeline configuration (uses defaults if None)
        """
        self.config = config or PM7Config()
        self.timeout_predictor = TimeoutPredictor(
            recalibrate_interval=self.config.timeout_predictor_recalibrate_interval
        )
        self.monitor = ThresholdMonitor(self.config)
        self.processor = MoleculeProcessor(self.config, self.timeout_predictor)
        self.evaluator = ResultEvaluator()
        self.results: list[PM7Result] = []
        self.logger: Optional[logging.Logger] = None
        self._paused = False

    def validate(self) -> bool:
        """Validate the execution environment.

        Returns:
            True if environment is valid
        """
        result = validate_environment(self.config)

        if not result.valid:
            for error in result.errors:
                LOG.error(f"Validation error: {error}")
            return False

        for warning in result.warnings:
            LOG.warning(f"Validation warning: {warning}")

        LOG.info("Environment validation passed")
        return True

    def _get_phase_value(self) -> str:
        """Get the phase value as a string, handling both enum and string types.

        Returns:
            Phase value as string
        """
        return (
            self.config.phase.value
            if isinstance(self.config.phase, Enum)
            else str(self.config.phase)
        )

    def setup(self, session_name: Optional[str] = None) -> None:
        """Set up pipeline for a processing session.

        Args:
            session_name: Optional session identifier
        """
        self.config.ensure_directories()
        self.logger = setup_logging(self.config, session_name)
        self.results = []
        self._paused = False
        LOG.info(f"Pipeline setup complete for phase {self._get_phase_value()}")

    def load_baseline(self, baseline_path: Path) -> bool:
        """Load baseline expectations for evaluation.

        Args:
            baseline_path: Path to baseline JSON file

        Returns:
            True if loaded successfully
        """
        return self.evaluator.load_baseline(baseline_path)

    def register_alert_callback(self, callback: Callable[[Alert], None]) -> None:
        """Register a callback for quality alerts.

        Args:
            callback: Function to call when alert is generated
        """
        self.monitor.register_callback(callback)

    def _handle_alerts(self, alerts: list[Alert]) -> None:
        """Handle alerts generated during molecule processing.

        Args:
            alerts: List of alerts to handle
        """
        if not alerts:
            return

        for alert in alerts:
            LOG.warning(f"Alert [{alert.level.value}] {alert.pattern}: {alert.message}")

    def process_molecule(
        self,
        mol_id: str,
        smiles: str,
        input_xyz: Optional[Path] = None,
    ) -> PM7Result:
        """Process a single molecule.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            input_xyz: Optional input XYZ file

        Returns:
            PM7Result
        """
        if self.logger:
            log_molecule_start(self.logger, mol_id, smiles)

        result = self.processor.process(mol_id, smiles, input_xyz)
        self.results.append(result)

        # Record in monitor (use direct comparison with enum members)
        timeout_occurred = any(
            c.mopac_status == MOPACStatus.TIMEOUT for c in result.conformers
        )
        scf_failed = any(
            c.mopac_status == MOPACStatus.SCF_FAILED for c in result.conformers
        )

        alerts = self.monitor.record_result(
            success=result.success,
            grade=result.quality_grade,
            hof_extracted=result.most_stable_hof is not None,
            timeout=timeout_occurred,
            scf_failed=scf_failed,
        )

        # Handle alerts returned by the monitor
        self._handle_alerts(alerts)

        if self.logger:
            log_molecule_complete(
                self.logger,
                mol_id,
                result.quality_grade.value,
                result.most_stable_hof,
                result.success,
            )

        # Check if should pause
        if self.monitor.should_pause():
            self._paused = True
            LOG.critical("Pipeline paused due to quality issues")

        return result

    def process_batch(
        self,
        molecules: list[tuple[str, str]],
    ) -> Iterator[PM7Result]:
        """Process a batch of molecules.

        Args:
            molecules: List of (mol_id, smiles) tuples

        Yields:
            PM7Result for each molecule
        """
        for mol_id, smiles in molecules:
            if self._paused:
                LOG.warning(f"Pipeline paused, skipping {mol_id}")
                continue

            result = self.process_molecule(mol_id, smiles)
            yield result

    def evaluate_phase_a(self) -> PhaseAEvaluation:
        """Evaluate Phase A results.

        Returns:
            PhaseAEvaluation
        """
        return self.evaluator.evaluate_phase_a(self.results)

    def save_results(self, output_path: Path) -> bool:
        """Save all results to JSON file.

        Args:
            output_path: Output file path

        Returns:
            True if saved successfully
        """
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)

            data = {
                "phase": self._get_phase_value(),
                "n_molecules": len(self.results),
                "monitor_summary": self.monitor.get_summary(),
                "timeout_predictor_stats": self.timeout_predictor.get_stats(),
                "results": [r.to_dict() for r in self.results],
            }

            with open(output_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)

            LOG.info(f"Results saved to {output_path}")
            return True

        except Exception as e:
            LOG.error(f"Failed to save results: {e}")
            return False

    def save_timeout_predictor(self, output_path: Optional[Path] = None) -> bool:
        """Save timeout predictor model.

        Args:
            output_path: Output path (default: config.output_dir/models/timeout_predictor.pkl)

        Returns:
            True if saved successfully
        """
        if output_path is None:
            output_path = (
                self.config.output_dir.parent / "models" / "timeout_predictor.pkl"
            )

        return self.timeout_predictor.save(output_path)

    def load_timeout_predictor(self, input_path: Path) -> bool:
        """Load timeout predictor model.

        Args:
            input_path: Input path

        Returns:
            True if loaded successfully
        """
        if self.timeout_predictor.load(input_path):
            self.processor.timeout_predictor = self.timeout_predictor
            return True
        return False

    def get_summary(self) -> dict:
        """Get pipeline execution summary.

        Uses single-pass aggregation for efficiency.

        Returns:
            Summary dictionary
        """
        # Single-pass aggregation
        n_successful = 0
        n_failed = 0
        grade_counts = {"A": 0, "B": 0, "C": 0, "FAILED": 0}

        for r in self.results:
            if r.success:
                n_successful += 1
            else:
                n_failed += 1

            # Map quality grade to count
            grade_name = r.quality_grade.name if r.quality_grade else "FAILED"
            if grade_name in grade_counts:
                grade_counts[grade_name] += 1
            else:
                # Track unexpected grades individually instead of collapsing to "UNKNOWN"
                grade_counts.setdefault(grade_name, 0)
                grade_counts[grade_name] += 1

        return {
            "phase": self._get_phase_value(),
            "n_molecules_processed": len(self.results),
            "n_successful": n_successful,
            "n_failed": n_failed,
            "success_rate": self.monitor.metrics.success_rate,
            "hof_extraction_rate": self.monitor.metrics.hof_extraction_rate,
            "grade_distribution": grade_counts,
            "alerts_generated": len(self.monitor.alerts),
            "paused": self._paused,
        }

    @property
    def is_paused(self) -> bool:
        """Whether pipeline is paused due to quality issues."""
        return self._paused

    def resume(self) -> None:
        """Resume pipeline after pause."""
        self._paused = False
        LOG.info("Pipeline resumed")
