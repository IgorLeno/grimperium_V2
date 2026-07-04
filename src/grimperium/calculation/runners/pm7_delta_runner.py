"""Canonical PM7 + Delta Learning runner (Method B).

Produces a :class:`MoleculeCalculationResult` carrying three estimates:

* ``BASELINE``   — the raw PM7 heat of formation,
* ``CORRECTION`` — the delta learning correction (``final - baseline``),
* ``FINAL``      — the model-predicted (delta-corrected) heat of formation.

The heavy scientific work is delegated to existing, tested components
(``CRESTPM7Pipeline``, ``build_pm7_delta_feature_frame``, the model persistence
helpers and the adapter). Dependencies are injected so the runner can be tested
without real CREST/xTB/MOPAC binaries.
"""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from grimperium.calculation.contracts.adapters import (
    PROPERTY_ID,
    pm7result_to_canonical,
)
from grimperium.calculation.contracts.adapters import (
    PM7Result as AdapterPM7Result,
)
from grimperium.calculation.contracts.enums import (
    PropertyRole,
    StageExecutionStatus,
)
from grimperium.calculation.contracts.models import (
    MoleculeCalculationResult,
    PropertyEstimate,
    StageExecutionRecord,
)
from grimperium.calculation.contracts.quantity import Quantity
from grimperium.calculation.methods import (
    CalculationMethodDefinition,
    get_calculation_method,
)
from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.molecule_processor import PM7Result
from grimperium.ml.persistence import load_model, load_model_metadata

METHOD_ID = "pm7_delta_learning"

PM7Processor = Callable[[str, str], PM7Result]
ModelLoader = Callable[[Path], tuple[Any, Any]]
MetadataLoader = Callable[[Path], dict[str, Any]]
FeatureBuilder = Callable[[str, PM7Result], pd.DataFrame]
ModelValidator = Callable[[Path, CalculationMethodDefinition], None]


class PM7DeltaRunnerError(Exception):
    """Raised when the PM7 + Delta Learning runner cannot produce a result."""


def _default_feature_builder(smiles: str, pm7_result: PM7Result) -> pd.DataFrame:
    # Imported lazily: calculation_features lives in the higher cli layer and
    # importing it eagerly would create a circular import at package load.
    from grimperium.cli.calculation_features import build_pm7_delta_feature_frame

    return build_pm7_delta_feature_frame(smiles, pm7_result)


def _default_model_validator(
    model_path: Path, method: CalculationMethodDefinition
) -> None:
    from grimperium.cli.model_compatibility import validate_model_for_method

    validate_model_for_method(model_path, method)


def _portable_model_path(model_path: Path) -> str:
    """Return a portable (relative, non-secret) representation of a model path."""
    try:
        return model_path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return model_path.name


class PM7DeltaLearningRunner:
    """Run Method B: PM7 baseline plus a delta-learning correction."""

    def __init__(
        self,
        *,
        config: PM7Config | None = None,
        pm7_processor: PM7Processor | None = None,
        model_loader: ModelLoader = load_model,
        metadata_loader: MetadataLoader = load_model_metadata,
        feature_builder: FeatureBuilder | None = None,
        model_validator: ModelValidator | None = None,
    ) -> None:
        self.config = config or PM7Config()
        self._pm7_processor = pm7_processor
        self.model_loader = model_loader
        self.metadata_loader = metadata_loader
        self.feature_builder = feature_builder or _default_feature_builder
        self.model_validator = model_validator or _default_model_validator

    def _process_pm7(self, molecule_id: str, smiles: str) -> PM7Result:
        if self._pm7_processor is not None:
            return self._pm7_processor(molecule_id, smiles)
        # Imported lazily so tests never require the heavy pipeline dependencies.
        from grimperium.crest_pm7.pipeline import CRESTPM7Pipeline

        return CRESTPM7Pipeline(self.config).process_molecule(molecule_id, smiles)

    def calculate_single_smiles(
        self,
        smiles: str,
        *,
        molecule_id: str,
        model_path: Path,
        method: CalculationMethodDefinition | None = None,
        name: str | None = None,
        charge: int = 0,
        multiplicity: int = 1,
        progress_callback: Callable[[str], None] | None = None,
    ) -> MoleculeCalculationResult:
        """Run PM7 + Delta Learning and return a canonical calculation result."""
        method = method or get_calculation_method(METHOD_ID)

        # Model must be validated before any inference is attributed to it.
        self.model_validator(model_path, method)

        if progress_callback:
            progress_callback("⏳ Generating geometry and conformers...")
        pm7_result = self._process_pm7(molecule_id, smiles)
        if not pm7_result.success:
            raise PM7DeltaRunnerError(
                f"PM7 pipeline failed: {pm7_result.error_message or 'unknown error'}"
            )

        h298_pm7 = pm7_result.most_stable_hof
        if h298_pm7 is None:
            raise PM7DeltaRunnerError("No stable conformer found — cannot extract HOF")

        if progress_callback:
            progress_callback("⚡ Loading model...")
        try:
            learner, feature_pipeline = self.model_loader(model_path)
            metadata = self.metadata_loader(model_path)
        except Exception as exc:  # noqa: BLE001 - surfaced as a runner error
            raise PM7DeltaRunnerError(f"Failed to load model: {exc}") from exc
        model_version = str(metadata.get("version", "unknown"))

        feature_started = datetime.now(timezone.utc)
        try:
            frame = self.feature_builder(smiles, pm7_result)
            features = feature_pipeline.transform(frame)
        except Exception as exc:  # noqa: BLE001 - surfaced as a runner error
            raise PM7DeltaRunnerError(f"Feature assembly failed: {exc}") from exc
        feature_completed = datetime.now(timezone.utc)

        if progress_callback:
            progress_callback("🧠 Computing delta correction...")
        inference_started = datetime.now(timezone.utc)
        try:
            # DeltaLearner.predict returns the FINAL H298, not the pure delta.
            y_pm7 = np.array([h298_pm7])
            h298_final = float(learner.predict(features, y_pm7)[0])
        except Exception as exc:  # noqa: BLE001 - surfaced as a runner error
            raise PM7DeltaRunnerError(f"Delta inference failed: {exc}") from exc
        inference_completed = datetime.now(timezone.utc)

        delta = h298_final - h298_pm7
        portable_path = _portable_model_path(model_path)
        conformer_source_id = self._selected_conformer_id(pm7_result)

        # The concrete PM7Result structurally satisfies the adapter Protocol at
        # runtime; the cast bridges nested-Protocol variance mypy cannot prove.
        result = pm7result_to_canonical(
            cast(AdapterPM7Result, pm7_result),
            method_id=METHOD_ID,
            method_version=method.version,
            property_role=PropertyRole.BASELINE,
        )
        # Honor the caller-provided molecule identity/charge state.
        if name is not None:
            result.molecule.name = name
        result.molecule.charge = charge
        result.molecule.multiplicity = multiplicity
        result.estimates.append(
            self._estimate(
                estimate_id=f"{molecule_id}_delta_correction",
                role=PropertyRole.CORRECTION,
                method_version=method.version,
                value=delta,
                conformer_source_id=conformer_source_id,
                model_path=portable_path,
            )
        )
        result.estimates.append(
            self._estimate(
                estimate_id=f"{molecule_id}_final",
                role=PropertyRole.FINAL,
                method_version=method.version,
                value=h298_final,
                conformer_source_id=conformer_source_id,
                model_path=portable_path,
            )
        )
        result.stage_executions.append(
            StageExecutionRecord(
                stage_id="feature_assembly",
                program="grimperium.cli.calculation_features",
                role="features",
                status=StageExecutionStatus.SUCCESS,
                requested=True,
                started_at=feature_started,
                completed_at=feature_completed,
                execution_time_s=(feature_completed - feature_started).total_seconds(),
                program_version=None,
                settings={"n_features": int(frame.shape[1])},
                artifact_ids=[],
                error_message=None,
            )
        )
        result.stage_executions.append(
            StageExecutionRecord(
                stage_id="delta_learning_inference",
                program="grimperium.ml",
                role="correction",
                status=StageExecutionStatus.SUCCESS,
                requested=True,
                started_at=inference_started,
                completed_at=inference_completed,
                execution_time_s=(
                    inference_completed - inference_started
                ).total_seconds(),
                program_version=model_version,
                settings={
                    "model_path": portable_path,
                    "model_version": model_version,
                    "feature_schema_id": metadata.get("feature_schema_id"),
                    "feature_schema_hash": metadata.get("feature_schema_hash"),
                    "feature_columns": self._feature_columns(
                        metadata, feature_pipeline
                    ),
                    "baseline_hamiltonian": method.compatibility.baseline_hamiltonian,
                },
                artifact_ids=[],
                error_message=None,
            )
        )
        return result

    @staticmethod
    def _selected_conformer_id(pm7_result: PM7Result) -> int | None:
        if pm7_result.k_selected_pm7 is not None:
            return pm7_result.k_selected_pm7
        selected = pm7_result.get_selected_conformer()
        return selected.crest_rank if selected is not None else None

    @staticmethod
    def _feature_columns(metadata: dict[str, Any], feature_pipeline: Any) -> list[str]:
        columns = metadata.get("feature_columns")
        if isinstance(columns, list):
            return [str(col) for col in columns]
        legacy = getattr(feature_pipeline, "feature_cols", None)
        if isinstance(legacy, list):
            return [str(col) for col in legacy]
        return []

    @staticmethod
    def _estimate(
        *,
        estimate_id: str,
        role: PropertyRole,
        method_version: str,
        value: float,
        conformer_source_id: int | None,
        model_path: str | None,
    ) -> PropertyEstimate:
        return PropertyEstimate(
            estimate_id=estimate_id,
            property_id=PROPERTY_ID,
            role=role,
            method_id=METHOD_ID,
            method_version=method_version,
            hamiltonian=None,
            value=Quantity(value=value, unit="kcal/mol"),
            value_kcal_mol=value,
            value_kj_mol=None,
            conformer_source_id=conformer_source_id,
            uncertainty=None,
            model_path=model_path,
        )
