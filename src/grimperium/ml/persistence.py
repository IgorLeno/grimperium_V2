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

# The saved bundle has this structure:
# {
#     "learner": DeltaLearner (fitted),
#     "pipeline": FeaturePipeline (fitted),
#     "metrics": {"train": DictStrAny, "test": DictStrAny},
#     "version": str,          e.g. "1.0.0"
#     "trained_at": str,       e.g. "2026-03-15T08:00:00+00:00"
# }


def save_model(bundle: dict[str, Any], path: Path) -> None:
    """Serialize DeltaLearner + FeaturePipeline + metrics to a joblib file.

    Args:
        bundle: dict with required keys: learner, pipeline, metrics.
        path: destination path (.joblib). Parent directories are created
              automatically.

    Raises:
        ValueError: if bundle is missing required keys.
    """
    required = {"learner", "pipeline", "metrics"}
    missing = required - bundle.keys()
    if missing:
        msg = f"bundle missing required keys: {missing}"
        raise ValueError(msg)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    payload: dict[str, Any] = {
        "learner": bundle["learner"],
        "pipeline": bundle["pipeline"],
        "metrics": bundle["metrics"],
        "version": MODEL_VERSION,
        "trained_at": datetime.now(timezone.utc).isoformat(),
    }
    joblib.dump(payload, path)
    logger.info("Model saved to %s (version %s)", path, MODEL_VERSION)


def load_model(path: Path) -> tuple[DeltaLearner, FeaturePipeline]:
    """Load DeltaLearner and FeaturePipeline from a joblib file.

    Args:
        path: path to the .joblib file.

    Returns:
        (DeltaLearner, FeaturePipeline) ready for predict and transform.

    Raises:
        FileNotFoundError: if the file does not exist.
        ValueError: if the file is not a valid Grimperium bundle.
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


def load_model_metadata(path: Path) -> DictStrAny:
    """Load only bundle metadata for CLI display.

    Note: joblib deserializes the entire bundle (including learner and
    pipeline). The 'metadata-only' separation is logical, not a
    performance optimization.

    Returns:
        dict with version, trained_at, metrics.

    Raises:
        FileNotFoundError: if the file does not exist.
    """
    path = Path(path)
    if not path.exists():
        msg = f"Model file not found: {path}"
        raise FileNotFoundError(msg)

    bundle: dict[str, Any] = joblib.load(path)
    return {
        "version": bundle.get("version", "unknown"),
        "trained_at": bundle.get("trained_at", "unknown"),
        "metrics": bundle.get("metrics", {}),
    }
