"""Model serialization and deserialization for Phase E."""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib

from grimperium import DictStrAny
from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.features import FeaturePipeline

logger = logging.getLogger(__name__)

MODEL_VERSION = "1.0.0"
DEFAULT_MODEL_PROPERTY_ID = "standard_enthalpy_of_formation"
DEFAULT_MODEL_BASELINE_HAMILTONIAN = "PM7"

# The saved bundle has this structure:
# {
#     "learner": DeltaLearner (fitted),
#     "pipeline": FeaturePipeline (fitted),
#     "metrics": {"train": DictStrAny, "test": DictStrAny},
#     "version": str,          e.g. "1.0.0"
#     "trained_at": str,       e.g. "2026-03-15T08:00:00+00:00"
# }


def save_model(bundle: dict[str, Any], path: Path | str) -> None:
    """Serialize DeltaLearner + FeaturePipeline + metrics to a joblib file.

    Args:
        bundle: dict with required keys: learner, pipeline, metrics.
        path: destination path (.joblib). Accepts ``Path`` or ``str``;
              parent directories are created automatically.

    Raises:
        ValueError: if bundle is missing required keys.

    Example:
        >>> save_model({"learner": lr, "pipeline": pl, "metrics": m},
        ...            Path("models/model.joblib"))
    """
    required = {"learner", "pipeline", "metrics"}
    missing = required - bundle.keys()
    if missing:
        msg = f"bundle missing required keys: {missing}"
        raise ValueError(msg)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    from grimperium.calculation.methods.feature_schema import (
        FEATURE_SCHEMA_HASH,
        FEATURE_SCHEMA_ID,
        get_feature_schema,
    )

    schema = get_feature_schema(FEATURE_SCHEMA_ID)
    payload: dict[str, Any] = {
        "learner": bundle["learner"],
        "pipeline": bundle["pipeline"],
        "metrics": bundle["metrics"],
        "version": MODEL_VERSION,
        "trained_at": datetime.now(timezone.utc).isoformat(),
        "property_id": bundle.get("property_id", DEFAULT_MODEL_PROPERTY_ID),
        "baseline_hamiltonian": bundle.get(
            "baseline_hamiltonian",
            DEFAULT_MODEL_BASELINE_HAMILTONIAN,
        ),
        "feature_schema_id": bundle.get("feature_schema_id", FEATURE_SCHEMA_ID),
        "feature_schema_hash": bundle.get(
            "feature_schema_hash",
            FEATURE_SCHEMA_HASH,
        ),
        "feature_columns": bundle.get("feature_columns", list(schema.columns)),
    }
    joblib.dump(payload, path)
    logger.info("Model saved to %s (version %s)", path, MODEL_VERSION)


def load_model(path: Path | str) -> tuple[DeltaLearner, FeaturePipeline]:
    """Load DeltaLearner and FeaturePipeline from a joblib file.

    Security: ``joblib.load`` deserializes objects and can execute
    arbitrary code.  Only load model files from trusted sources.

    Args:
        path: path to the ``.joblib`` file.  Accepts ``Path`` or ``str``.

    Returns:
        (DeltaLearner, FeaturePipeline) ready for predict and transform.

    Raises:
        FileNotFoundError: if the file does not exist.
        ValueError: if the file is not a valid Grimperium bundle.

    Example:
        >>> learner, pipeline = load_model("models/model.joblib")
    """
    path = Path(path)
    if not path.exists():
        msg = f"Model file not found: {path}"
        raise FileNotFoundError(msg)

    bundle: dict[str, Any] = joblib.load(path)

    if "learner" not in bundle or "pipeline" not in bundle:
        msg = (
            f"Invalid model bundle at {path}. "
            "Expected keys: learner, pipeline, metrics, version, trained_at"
        )
        raise ValueError(msg)

    learner: DeltaLearner = bundle["learner"]
    pipeline: FeaturePipeline = bundle["pipeline"]
    logger.info(
        "Model loaded from %s (version %s, trained at %s)",
        path,
        bundle.get("version", "unknown"),
        bundle.get("trained_at", "unknown"),
    )
    return learner, pipeline


def load_model_metadata(path: Path | str) -> DictStrAny:
    """Load only bundle metadata for CLI display.

    Note: joblib deserializes the entire bundle (including learner and
    pipeline). The 'metadata-only' separation is logical, not a
    performance optimization.

    Args:
        path: path to the ``.joblib`` file.  Accepts ``Path`` or ``str``.

    Returns:
        dict with version, trained_at, metrics, n_train, n_test, n_total,
        and test_size.

    Raises:
        FileNotFoundError: if the file does not exist.

    Example:
        >>> meta = load_model_metadata(Path("models/model.joblib"))
        >>> meta["version"]
        '1.0.0'
    """
    path = Path(path)
    if not path.exists():
        msg = f"Model file not found: {path}"
        raise FileNotFoundError(msg)

    bundle: dict[str, Any] = joblib.load(path)
    metrics = bundle.get("metrics", {})
    train_m = metrics.get("train", {}) if isinstance(metrics, dict) else {}
    test_m = metrics.get("test", {}) if isinstance(metrics, dict) else {}

    return {
        "version": bundle.get("version", "unknown"),
        "trained_at": bundle.get("trained_at", "unknown"),
        "metrics": metrics,
        "n_train": train_m.get("n_samples"),
        "n_test": test_m.get("n_samples"),
        "n_total": test_m.get("n_total"),
        "test_size": test_m.get("test_size"),
        "property_id": bundle.get("property_id"),
        "baseline_hamiltonian": bundle.get("baseline_hamiltonian"),
        "feature_schema_id": bundle.get("feature_schema_id"),
        "feature_schema_hash": bundle.get("feature_schema_hash"),
        "feature_columns": bundle.get("feature_columns"),
    }
