"""Batch execution manager for CREST PM7 pipeline.

This module provides BatchExecutionManager for:
- Orchestrating batch execution
- Coordinating CSV tracking, detail files, and processing
- Managing failure policies (PARTIAL_OK, ALL_OR_NOTHING)
- Generating BatchResult with statistics
"""

import json
import logging
import time
from collections.abc import Callable
from datetime import datetime, timezone
from inspect import Parameter, signature
from pathlib import Path
from typing import Any, cast

from grimperium.calculation.contracts.adapters import PM7Result as AdapterPM7Result
from grimperium.calculation.contracts.adapters import (
    canonical_pm7_method_id,
    pm7result_to_canonical,
)
from grimperium.calculation.contracts.models import MoleculeCalculationResult
from grimperium.crest_pm7.batch.artifact_manager import ArtifactManager
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.detail_manager import ConformerDetailManager
from grimperium.crest_pm7.batch.enums import BatchFailurePolicy, MoleculeStatus
from grimperium.crest_pm7.batch.models import Batch, BatchResult
from grimperium.crest_pm7.batch.output_contracts import (
    BatchOutputLayout,
    BatchResultWriteReport,
    write_batch_calculation_results,
)
from grimperium.crest_pm7.batch.processor_adapter import FixedTimeoutProcessor
from grimperium.crest_pm7.batch.result_applier import BatchResultApplier
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.csv_enhancements import (
    BatchSettingsCapture,
    CSVManagerExtensions,
)
from grimperium.crest_pm7.logging_enhancements import (
    setup_batch_logging,
    suppress_pandas_warnings,
)
from grimperium.crest_pm7.molecule_processor import ConformerData
from grimperium.crest_pm7.progress import (
    CREST_STATUS_PREOPT,
    CREST_STATUS_SEARCH,
    MOPAC_STATUS_RUNNING,
    BatchProgressStage,
)

LOG = logging.getLogger("grimperium.crest_pm7.batch.execution_manager")

ResultWriter = Callable[..., BatchResultWriteReport]


class BatchExecutionManager:
    """Orchestrates batch execution of CREST PM7 pipeline.

    Coordinates:
    - BatchStateManager: Operational state tracking
    - BatchCSVManager: Scientific legacy CSV writes
    - ConformerDetailManager: Per-molecule JSON files
    - FixedTimeoutProcessor: Molecule processing with fixed timeouts

    Implements two failure policies:
    - PARTIAL_OK: Failed molecules marked RERUN/SKIP individually
    - ALL_OR_NOTHING: If any fail, reset entire batch to PENDING

    Attributes:
        csv_manager: Manager for scientific CSV writes
        state_manager: Manager for operational state tracking
        detail_manager: Manager for JSON detail files
        pm7_config: PM7 configuration
        processor_adapter: Fixed timeout processor
    """

    def __init__(
        self,
        csv_manager: BatchCSVManager,
        state_manager: BatchStateManager,
        detail_manager: ConformerDetailManager,
        pm7_config: PM7Config,
        processor_adapter: FixedTimeoutProcessor | None = None,
        artifact_manager: ArtifactManager | None = None,
        output_layout: BatchOutputLayout | None = None,
        result_writer: ResultWriter = write_batch_calculation_results,
        write_xlsx: bool = True,
    ) -> None:
        """Initialize batch execution manager.

        Args:
            csv_manager: Manager for scientific CSV writes
            state_manager: Manager for operational state tracking
            detail_manager: Manager for JSON detail files
            pm7_config: PM7 configuration
            processor_adapter: Optional processor adapter (created if None)
            artifact_manager: Optional artifact manager for debug/audit files
            output_layout: Optional layout enabling canonical result files
                (``calculation_results.csv``/``.xlsx``). ``None`` disables them.
            result_writer: Injectable writer for canonical batch results.
            write_xlsx: Whether to also emit the XLSX canonical result file.
        """
        self.csv_manager = csv_manager
        self.state_manager = state_manager
        self.detail_manager = detail_manager
        self.pm7_config = pm7_config
        self.processor_adapter = processor_adapter or FixedTimeoutProcessor(pm7_config)
        self.artifact_manager = artifact_manager
        self._output_layout = output_layout
        self._result_writer = result_writer
        self._write_xlsx = write_xlsx
        self._canonical_results: list[MoleculeCalculationResult] = []
        self.result_applier = BatchResultApplier(
            state_manager=state_manager,
            csv_manager=csv_manager,
        )

        LOG.info("BatchExecutionManager initialized")

    def execute_batch(
        self,
        batch: Batch,
        progress_callback: Callable[[str, int, int], None] | None = None,
    ) -> BatchResult:
        """Execute a batch of molecules.

        Processes each molecule sequentially, updates CSV status,
        saves detail files, and handles failures according to policy.

        Args:
            batch: Batch to execute
            progress_callback: Optional callback(mol_id, current, total) for progress

        Returns:
            BatchResult with execution statistics
        """
        if batch.is_empty:
            LOG.warning(f"Batch {batch.batch_id} is empty, nothing to execute")
            return BatchResult(
                batch_id=batch.batch_id,
                total_count=0,
            )

        LOG.info(
            f"Starting batch {batch.batch_id}: {batch.size} molecules, "
            f"policy={batch.failure_policy.value}"
        )

        # Setup structured logging for this batch
        batch_logger = setup_batch_logging(batch.batch_id)
        batch_logger.info(
            f"🚀 Starting batch {batch.batch_id}: {batch.size} molecules, "
            f"policy={batch.failure_policy.value}"
        )

        # Suppress pandas DtypeWarning and FutureWarning
        suppress_pandas_warnings()

        # Capture batch settings for CSV population
        batch_settings = BatchSettingsCapture.capture_batch_settings(self.pm7_config)
        batch_logger.debug(f"Batch settings: {batch_settings}")

        # Update processor timeouts
        self.processor_adapter.update_timeouts(
            crest_timeout_minutes=batch.crest_timeout_minutes,
            mopac_timeout_minutes=batch.mopac_timeout_minutes,
        )

        # Initialize result
        result = BatchResult(
            batch_id=batch.batch_id,
            total_count=batch.size,
            timestamp_start=datetime.now(timezone.utc),
        )

        # Store logger and settings for use in _process_molecule
        self._batch_logger = batch_logger
        self._batch_settings = batch_settings

        # Track HOF values for statistics
        hof_values: list[tuple[str, float]] = []  # (mol_id, hof)
        start_time = time.time()

        # Reset the canonical scientific accumulator for this batch run.
        self._canonical_results = []

        # Process each molecule
        for i, mol in enumerate(batch.molecules, start=1):
            if progress_callback:
                try:
                    progress_callback(mol.mol_id, i - 1, batch.size)
                except Exception as e:
                    LOG.warning(f"Progress callback error (pre): {e}")

            try:
                self._process_molecule(
                    mol_id=mol.mol_id,
                    smiles=mol.smiles,
                    batch_id=batch.batch_id,
                    batch_order=mol.batch_order,
                    charge=mol.charge,
                    multiplicity=mol.multiplicity,
                    crest_timeout=batch.crest_timeout_minutes,
                    mopac_timeout=batch.mopac_timeout_minutes,
                    result=result,
                    hof_values=hof_values,
                    method_id=batch.method_id,
                    method_version=batch.method_version,
                    method_snapshot=batch.method_snapshot,
                    batch_failure_policy=batch.failure_policy.value,
                )
            except Exception as e:
                LOG.error(
                    f"Unexpected error processing {mol.mol_id}: {e}",
                    exc_info=True,
                )
                result.failed_count += 1
                result.failed_mol_ids.append(mol.mol_id)

            # Progress callback
            if progress_callback:
                try:
                    progress_callback(mol.mol_id, i, batch.size)
                except Exception as e:
                    LOG.warning(f"Progress callback error (post): {e}")

        # Persist canonical scientific results (partial successes included) before
        # any later finalization step can fail. batch_state.csv stays purely
        # operational; the canonical files hold only scientific data.
        self._persist_canonical_results()

        # Finalize timing
        result.total_time = time.time() - start_time
        result.timestamp_end = datetime.now(timezone.utc)

        # Calculate HOF statistics
        if hof_values:
            # Find min and max simultaneously to avoid floating-point equality issues
            min_mol_id, min_hof_val = min(hof_values, key=lambda x: x[1])
            result.min_hof = round(min_hof_val, 2)
            result.min_hof_mol_id = min_mol_id

            max_mol_id, max_hof_val = max(hof_values, key=lambda x: x[1])
            result.max_hof = round(max_hof_val, 2)
            result.max_hof_mol_id = max_mol_id

        # Handle ALL_OR_NOTHING policy
        if batch.failure_policy == BatchFailurePolicy.ALL_OR_NOTHING and (
            result.failed_count > 0 or result.rerun_count > 0
        ):
            LOG.warning(
                f"ALL_OR_NOTHING: Resetting batch {batch.batch_id} due to failures"
            )
            batch_mol_ids = self.state_manager.get_molecules_by_batch_id(batch.batch_id)
            reset_count = self.state_manager.reset_all_or_nothing(batch.batch_id)
            for mol_id in batch_mol_ids:
                self.csv_manager.apply_operational_status(
                    mol_id,
                    MoleculeStatus.PENDING.value,
                )
            result.invalidated = bool(reset_count)
            LOG.info(f"Reset {reset_count} molecules from batch {batch.batch_id}")

        LOG.info(
            f"Batch {batch.batch_id} complete: "
            f"{result.success_count}/{result.total_count} OK, "
            f"{result.rerun_count} rerun, {result.skip_count} skip, "
            f"{result.failed_count} failed, "
            f"time={result.total_time:.1f}s"
        )

        return result

    def _process_molecule(
        self,
        mol_id: str,
        smiles: str,
        batch_id: str,
        batch_order: int,
        crest_timeout: float,
        mopac_timeout: float,
        result: BatchResult,
        hof_values: list[tuple[str, float]],
        charge: int = 0,
        multiplicity: int = 1,
        method_id: str = "pm7_delta_learning",
        method_version: str = "1.0.0",
        method_snapshot: dict[str, Any] | None = None,
        batch_failure_policy: str = "",
    ) -> None:
        """Process a single molecule within batch context.

        Updates CSV status, saves detail file, and updates result statistics.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            batch_id: Current batch ID
            batch_order: Position in batch
            crest_timeout: CREST timeout (minutes)
            mopac_timeout: MOPAC timeout (minutes)
            result: BatchResult to update
            hof_values: List to append HOF values
            method_id: Calculation method registry identifier
            method_version: Calculation method version
            method_snapshot: Validated method definition snapshot
            batch_failure_policy: Failure policy value recorded in state
        """
        LOG.info(f"Processing {mol_id} ({batch_order}/{result.total_count})")

        # Safe logger access with fallback chain
        logger = getattr(
            self,
            "_batch_logger",
            getattr(self, "_logger", logging.getLogger(__name__)),
        )
        logger.info(f"[{mol_id}] Processing ({batch_order}/{result.total_count})")

        # TODO: Add RDKit logging if RDKit processing happens here
        # log_rdkit_start(logger, mol_id)
        # log_rdkit_done(logger, mol_id, nrotbonds=X, tpsa=Y, aromatic_rings=Z)

        # TODO: Add CREST logging in processor_adapter or here if accessible
        # log_crest_start(logger, mol_id)
        # log_crest_done(logger, mol_id, num_conformers=X, time_seconds=Y)

        # TODO: Add MOPAC logging in processor_adapter or here if accessible
        # log_mopac_start(logger, mol_id, num_conformers=X)
        # log_mopac_done(logger, mol_id, best_conformer_idx=X, best_delta_energy=Y, time_seconds=Z)

        state_fields = {
            "smiles": smiles,
            "batch_id": batch_id,
            "batch_order": batch_order,
            "charge": charge,
            "multiplicity": multiplicity,
            "batch_failure_policy": batch_failure_policy,
            "assigned_crest_timeout": crest_timeout,
            "assigned_mopac_timeout": mopac_timeout,
            "method_id": method_id,
            "method_version": method_version,
            "method_definition_snapshot": self._serialize_method_snapshot(
                method_snapshot or {}
            ),
        }
        self.state_manager.update_molecule_status(
            mol_id,
            MoleculeStatus.RUNNING.value,
            extra_fields=state_fields,
        )

        def progress_callback(stage: BatchProgressStage) -> None:
            """Update operational progress stage for current molecule."""
            try:
                if stage == BatchProgressStage.XTB_PREOPT:
                    self.state_manager.update_molecule_status(
                        mol_id,
                        MoleculeStatus.RUNNING.value,
                        extra_fields={"crest_status": CREST_STATUS_PREOPT},
                    )
                elif stage == BatchProgressStage.CREST_SEARCH:
                    self.state_manager.update_molecule_status(
                        mol_id,
                        MoleculeStatus.RUNNING.value,
                        extra_fields={"crest_status": CREST_STATUS_SEARCH},
                    )
                elif stage == BatchProgressStage.MOPAC_CALC:
                    self.state_manager.update_molecule_status(
                        mol_id,
                        MoleculeStatus.RUNNING.value,
                        extra_fields={"mopac_status": MOPAC_STATUS_RUNNING},
                    )
            except Exception as exc:
                logger.warning(
                    f"[{mol_id}] Failed to update progress stage {stage.value}: {exc}"
                )

        try:
            # Process molecule
            pm7_result = self._process_with_charge_metadata(
                mol_id=mol_id,
                smiles=smiles,
                progress_callback=progress_callback,
                charge=charge,
                multiplicity=multiplicity,
            )

            # Save detail file
            detail = self.detail_manager.pm7result_to_detail(
                mol_id=mol_id,
                smiles=smiles,
                result=pm7_result,
                batch_id=batch_id,
            )
            self.detail_manager.save_detail(detail)

            # Save artifacts for debug/audit if artifact manager is configured
            if self.artifact_manager is not None:
                try:
                    # Work directories are typically temp_dir/{mol_id}/crest and /mopac
                    mol_work_dir = self.pm7_config.temp_dir / mol_id
                    crest_work_dir = (
                        mol_work_dir / "crest" if mol_work_dir.exists() else None
                    )
                    mopac_work_dir = (
                        mol_work_dir / "mopac" if mol_work_dir.exists() else None
                    )

                    self.artifact_manager.save_artifacts(
                        mol_id=mol_id,
                        batch_id=batch_id,
                        crest_work_dir=crest_work_dir,
                        mopac_work_dir=mopac_work_dir,
                        success=pm7_result.success,
                        extra_metadata={
                            "quality_grade": pm7_result.quality_grade.value,
                            "most_stable_hof": pm7_result.most_stable_hof,
                            "error_message": pm7_result.error_message,
                        },
                    )
                except Exception as e:
                    LOG.warning(f"Failed to save artifacts for {mol_id}: {e}")

            # Update status based on success
            if pm7_result.success:
                csv_update = self.csv_manager.pm7result_to_csv_update(
                    mol_id=mol_id,
                    result=pm7_result,
                    batch_id=batch_id,
                    batch_order=batch_order,
                    crest_timeout_used=crest_timeout,
                    mopac_timeout_used=mopac_timeout,
                )
                self.result_applier.apply_success(mol_id, csv_update)
                result.success_count += 1

                # Accumulate the canonical scientific result when canonical output
                # is enabled. The batch executes the PM7-only pipeline, so only a
                # BASELINE PM7 estimate is produced here — no delta
                # CORRECTION/FINAL is invented.
                if self._output_layout is not None:
                    self._canonical_results.append(
                        pm7result_to_canonical(
                            cast(AdapterPM7Result, pm7_result),
                            method_id=canonical_pm7_method_id(method_id),
                            method_version=method_version,
                        )
                    )

                # Enhance CSV with delta calculations and batch settings
                # h298_cbs comes from CSV original data (reference_hof), NOT csv_update
                h298_cbs = self.csv_manager.get_reference_hof(mol_id)
                h298_pm7 = pm7_result.most_stable_hof  # Property, may be None

                # Get PM7-selected conformer and its CREST rank
                selected_conformer: ConformerData | None = (
                    pm7_result.get_selected_conformer()
                )
                k_selected_pm7: int | None = pm7_result.k_selected_pm7

                # Safe access to batch_settings
                batch_settings = getattr(self, "_batch_settings", {})

                # Update CSV with descriptors and settings
                success = CSVManagerExtensions.update_molecule_with_mopac_results(
                    csv_manager=self.csv_manager,
                    mol_id=mol_id,
                    h298_cbs=h298_cbs,
                    h298_pm7=h298_pm7,
                    selected_conformer=selected_conformer,
                    k_selected_pm7=k_selected_pm7,
                    batch_settings=batch_settings,
                )

                if success:
                    logger.info(
                        f"[{mol_id}] CSV enhanced with descriptors and settings"
                    )
                else:
                    logger.warning(f"[{mol_id}] CSV enhancement failed")

                # Track HOF for statistics
                if pm7_result.most_stable_hof is not None:
                    hof_values.append((mol_id, pm7_result.most_stable_hof))

                LOG.info(
                    f"{mol_id}: OK (HOF={pm7_result.most_stable_hof}, "
                    f"grade={pm7_result.quality_grade.value})"
                )
            else:
                # Record retry state, or terminal skip when attempts are exhausted.
                error_msg = pm7_result.error_message or "Unknown error"
                decision = self.result_applier.apply_failure(mol_id, error_msg)
                new_status = decision.final_status

                if new_status == MoleculeStatus.SKIP.value:
                    result.skip_count += 1
                else:
                    result.rerun_count += 1
                    result.rerun_mol_ids.append(mol_id)

                LOG.warning(f"{mol_id}: {new_status} - {error_msg}")

        except Exception as e:
            # Unexpected error - record retry state.
            error_msg = f"Processing exception: {str(e)}"
            LOG.error(f"{mol_id}: {error_msg}", exc_info=True)

            decision = self.result_applier.apply_failure(mol_id, error_msg)
            new_status = decision.final_status

            if new_status == MoleculeStatus.SKIP.value:
                result.skip_count += 1
            else:
                result.rerun_count += 1
                result.rerun_mol_ids.append(mol_id)

    def status_summary(self) -> dict[str, int]:
        """Return current operational status counts.

        Returns:
            Dict mapping status to count
        """
        return self.state_manager.status_counts()

    def get_pending_count(self) -> int:
        """Get count of molecules still pending.

        Returns:
            Number of PENDING + RERUN molecules
        """
        counts = self.status_summary()
        return counts.get(MoleculeStatus.PENDING.value, 0) + counts.get(
            MoleculeStatus.RERUN.value, 0
        )

    def is_complete(self) -> bool:
        """Check if all molecules have been processed.

        Returns:
            True if no PENDING or RERUN molecules remain
        """
        return self.get_pending_count() == 0

    def _persist_canonical_results(self) -> None:
        """Write accumulated canonical results to calculation_results.csv/.xlsx.

        No-op when no output layout was configured. The CSV is written before the
        XLSX, so an XLSX failure never destroys a valid CSV. Output failures are
        logged clearly rather than raised, so a failed export cannot discard the
        already-persisted scientific/operational state.
        """
        if self._output_layout is None:
            return

        layout = self._output_layout
        layout.output_dir.mkdir(parents=True, exist_ok=True)
        try:
            report = self._result_writer(
                self._canonical_results,
                layout,
                include_xlsx=self._write_xlsx,
            )
            LOG.info(
                "Canonical batch results written: %d result(s) -> %s",
                report.result_count,
                report.csv_path,
            )
        except Exception as exc:
            LOG.error(
                "Failed to write canonical batch results to %s: %s",
                layout.output_dir,
                exc,
                exc_info=True,
            )

    def _process_with_charge_metadata(
        self,
        *,
        mol_id: str,
        smiles: str,
        progress_callback: Callable[[BatchProgressStage], None],
        charge: int,
        multiplicity: int,
    ) -> Any:
        """Call the processor with charge metadata when its adapter supports it."""
        processor_method = self.processor_adapter.process_with_fixed_timeout
        parameters = signature(processor_method).parameters.values()
        accepts_kwargs = any(
            parameter.kind is Parameter.VAR_KEYWORD for parameter in parameters
        )
        names = {parameter.name for parameter in parameters}
        if accepts_kwargs or {"charge", "multiplicity"}.issubset(names):
            return processor_method(
                mol_id=mol_id,
                smiles=smiles,
                progress_callback=progress_callback,
                charge=charge,
                multiplicity=multiplicity,
            )
        return processor_method(
            mol_id=mol_id,
            smiles=smiles,
            progress_callback=progress_callback,
        )

    def _serialize_method_snapshot(self, method_snapshot: dict[str, Any]) -> str:
        """Serialize method metadata for the operational state CSV."""
        return json.dumps(method_snapshot, sort_keys=True)


def create_execution_manager(
    csv_path: Path,
    state_csv_path: Path | None,
    detail_dir: Path,
    pm7_config: PM7Config,
    crest_timeout_minutes: float = 30.0,
    mopac_timeout_minutes: float = 60.0,
    artifact_dir: Path | None = None,
    preserve_artifacts_on_success: bool = True,
    preserve_artifacts_on_failure: bool = True,
    output_dir: Path | None = None,
    write_xlsx: bool = True,
) -> BatchExecutionManager:
    """Factory function to create BatchExecutionManager with defaults.

    Args:
        csv_path: Path to legacy scientific CSV file
        state_csv_path: Path to split operational state CSV
        detail_dir: Directory for JSON detail files
        pm7_config: PM7 configuration
        crest_timeout_minutes: Default CREST timeout
        mopac_timeout_minutes: Default MOPAC timeout
        artifact_dir: Directory for debug/audit artifacts (None to disable)
        preserve_artifacts_on_success: Save artifacts for successful molecules
        preserve_artifacts_on_failure: Save artifacts for failed molecules
        output_dir: Directory for canonical result files (None to disable)
        write_xlsx: Also emit the XLSX canonical result file

    Returns:
        Configured BatchExecutionManager
    """
    csv_manager = BatchCSVManager(csv_path)
    resolved_state_csv_path = state_csv_path or csv_path.parent / "batch_state.csv"
    state_manager = BatchStateManager(resolved_state_csv_path, pm7_config)
    detail_manager = ConformerDetailManager(detail_dir)
    processor = FixedTimeoutProcessor(
        config=pm7_config,
        crest_timeout_minutes=crest_timeout_minutes,
        mopac_timeout_minutes=mopac_timeout_minutes,
    )

    # Create artifact manager if directory specified
    artifact_manager = None
    if artifact_dir is not None:
        artifact_manager = ArtifactManager(
            artifact_dir=artifact_dir,
            preserve_on_success=preserve_artifacts_on_success,
            preserve_on_failure=preserve_artifacts_on_failure,
        )

    output_layout = BatchOutputLayout(output_dir) if output_dir is not None else None

    return BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=detail_manager,
        pm7_config=pm7_config,
        processor_adapter=processor,
        artifact_manager=artifact_manager,
        output_layout=output_layout,
        write_xlsx=write_xlsx,
    )
