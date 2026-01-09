"""Timeout prediction for MOPAC calculations.

Uses Huber regression to predict optimal timeouts based on molecule complexity.

Note: Saved models should only be loaded from trusted sources.
This module uses joblib for serialization instead of pickle for better sklearn
compatibility and security. A SHA256 checksum is stored alongside the model
to detect tampering.
"""

import hashlib
import logging
import numbers
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import joblib
import numpy as np
from sklearn.linear_model import HuberRegressor

from .config import TimeoutConfidence

LOG = logging.getLogger("grimperium.crest_pm7.timeout_predictor")

# Constants
MIN_SAMPLES_FOR_FIT = 20  # Minimum samples to train model
TIMEOUT_MIN = 60.0  # Minimum timeout (1 minute)
TIMEOUT_MAX = 3600.0  # Maximum timeout (1 hour)


def _clamp(value: float, min_val: float = TIMEOUT_MIN, max_val: float = TIMEOUT_MAX) -> float:
    """Clamp value to [min, max] range."""
    return max(min_val, min(max_val, value))


def _compute_checksum(filepath: Path) -> str:
    """Compute SHA256 checksum of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


@dataclass
class TimeoutPredictor:
    """Predicts optimal MOPAC timeout based on molecule complexity.

    Uses Huber regression trained on observed execution times.
    Operates in three phases:
    - Phase A (<20 samples): Pure heuristic
    - Phase B (20-50 samples): Model + 50% margin
    - Phase C (50+ samples): Model + 30% margin

    Attributes:
        model: HuberRegressor model
        samples_nheavy: List of nheavy values from observations
        samples_time: List of execution times from observations
        recalibrate_interval: Molecules between recalibration
        epsilon: Huber regressor epsilon parameter
        alpha: Huber regressor alpha parameter
        max_iter: Huber regressor max iterations
        tol: Huber regressor tolerance
    """

    model: Optional[HuberRegressor] = None
    samples_nheavy: list[int] = field(default_factory=list)
    samples_time: list[float] = field(default_factory=list)
    recalibrate_interval: int = 50
    _molecules_since_recalibration: int = 0

    # HuberRegressor hyperparameters
    epsilon: float = 1.35
    alpha: float = 0.0001
    max_iter: int = 100
    tol: float = 1e-5

    @property
    def n_samples(self) -> int:
        """Number of training samples.

        Returns the minimum length if sample lists differ (fail-safe).
        """
        len_nheavy = len(self.samples_nheavy)
        len_time = len(self.samples_time)
        if len_nheavy != len_time:
            LOG.warning(
                f"Sample lists are out of sync: "
                f"nheavy has {len_nheavy}, time has {len_time}. "
                f"Returning min({len_nheavy}, {len_time})."
            )
            return min(len_nheavy, len_time)
        return len_nheavy

    @property
    def is_fitted(self) -> bool:
        """Whether the model has been fitted."""
        return self.model is not None

    def add_observation(self, nheavy: int, execution_time: float) -> None:
        """Add a new observation for model training.

        Args:
            nheavy: Number of heavy atoms (must be non-negative int)
            execution_time: Actual execution time in seconds (must be positive float)

        Raises:
            ValueError: If inputs are invalid
        """
        # Validate inputs (accept all integer-like and real-like types)
        if not isinstance(nheavy, numbers.Integral) or nheavy < 0:
            raise ValueError(
                f"nheavy must be a non-negative integer, got {nheavy} (type: {type(nheavy).__name__})"
            )
        if not isinstance(execution_time, numbers.Real) or execution_time <= 0:
            raise ValueError(
                f"execution_time must be a positive number, got {execution_time}"
            )

        self.samples_nheavy.append(nheavy)
        self.samples_time.append(float(execution_time))
        self._molecules_since_recalibration += 1

        LOG.debug(f"Added observation: nheavy={nheavy}, time={execution_time:.1f}s")

        # Check if recalibration needed - always reset counter
        if self._molecules_since_recalibration >= self.recalibrate_interval:
            self._molecules_since_recalibration = 0  # Reset counter before fitting
            if self.n_samples >= MIN_SAMPLES_FOR_FIT:
                self.fit()

    def fit(self) -> bool:
        """Fit the Huber regression model.

        Returns:
            True if fitting succeeded, False otherwise
        """
        if self.n_samples < MIN_SAMPLES_FOR_FIT:
            LOG.debug(f"Not enough samples to fit: {self.n_samples} < {MIN_SAMPLES_FOR_FIT}")
            return False

        try:
            X = np.array(self.samples_nheavy).reshape(-1, 1)
            y = np.array(self.samples_time)

            self.model = HuberRegressor(
                epsilon=self.epsilon,
                alpha=self.alpha,
                max_iter=self.max_iter,
                tol=self.tol,
            )
            self.model.fit(X, y)
            self._molecules_since_recalibration = 0

            LOG.info(f"Model fitted with {self.n_samples} samples")
            return True

        except Exception as e:
            LOG.warning(f"Model fitting failed: {e}")
            return False

    def _predict_with_model(
        self,
        nheavy: int,
        num_conformers: int,
        margin: float,
        confidence: TimeoutConfidence,
    ) -> tuple[float, TimeoutConfidence]:
        """Predict timeout using the model with validation and fallback.

        Args:
            nheavy: Number of heavy atoms
            num_conformers: Number of conformers to process
            margin: Safety margin multiplier
            confidence: Confidence level to return on success

        Returns:
            Tuple of (timeout_seconds, confidence)
        """
        pred = self.model.predict([[nheavy]])[0]

        # Validate prediction
        if not np.isfinite(pred) or pred < 0:
            LOG.warning(f"Invalid prediction {pred} for nheavy={nheavy}, using heuristic fallback")
            timeout = (120 + 5 * nheavy) * (1 + 0.2 * (num_conformers - 1))
            return _clamp(timeout), TimeoutConfidence.LOW

        timeout = pred * num_conformers * margin
        return _clamp(timeout), confidence

    def predict(self, nheavy: int, num_conformers: int = 1) -> tuple[float, TimeoutConfidence]:
        """Predict optimal timeout for a molecule.

        Args:
            nheavy: Number of heavy atoms
            num_conformers: Number of conformers to process

        Returns:
            Tuple of (timeout_seconds, confidence)
        """
        n = self.n_samples

        if n < MIN_SAMPLES_FOR_FIT or not self.is_fitted:
            # Phase A: Pure heuristic
            timeout = (120 + 5 * nheavy) * (1 + 0.2 * (num_conformers - 1))
            return _clamp(timeout), TimeoutConfidence.LOW

        elif n < 50:
            # Phase B: Model + 50% margin
            return self._predict_with_model(nheavy, num_conformers, 1.5, TimeoutConfidence.MEDIUM)

        else:
            # Phase C+: Model + 30% margin
            return self._predict_with_model(nheavy, num_conformers, 1.3, TimeoutConfidence.HIGH)

    def save(self, filepath: Path) -> bool:
        """Save predictor state to file using joblib.

        Also saves a SHA256 checksum file for tamper detection.

        Args:
            filepath: Path to save file

        Returns:
            True if save succeeded
        """
        try:
            filepath.parent.mkdir(parents=True, exist_ok=True)
            state = {
                "model": self.model,
                "samples_nheavy": self.samples_nheavy,
                "samples_time": self.samples_time,
                "recalibrate_interval": self.recalibrate_interval,
            }
            joblib.dump(state, filepath, compress=3)

            # Save checksum for tamper detection
            checksum = _compute_checksum(filepath)
            checksum_file = filepath.with_suffix(".sha256")
            checksum_file.write_text(checksum + "\n")

            LOG.info(f"Predictor saved to {filepath}")
            return True
        except Exception as e:
            LOG.warning(f"Failed to save predictor: {e}")
            return False

    def load(self, filepath: Path) -> bool:
        """Load predictor state from file using joblib.

        Verifies SHA256 checksum before loading.

        Args:
            filepath: Path to saved file

        Returns:
            True if load succeeded
        """
        try:
            # Verify checksum if exists
            checksum_file = filepath.with_suffix(".sha256")
            if checksum_file.exists():
                expected_checksum = checksum_file.read_text().strip()
                actual_checksum = _compute_checksum(filepath)
                if expected_checksum != actual_checksum:
                    LOG.warning(f"Checksum mismatch for {filepath}, file may be corrupted")
                    return False

            state = joblib.load(filepath)

            # Validate loaded state
            loaded_nheavy = state.get("samples_nheavy", [])
            loaded_time = state.get("samples_time", [])

            # Verify lists have equal length
            if len(loaded_nheavy) != len(loaded_time):
                LOG.warning(
                    f"Corrupted state: samples_nheavy ({len(loaded_nheavy)}) "
                    f"!= samples_time ({len(loaded_time)})"
                )
                return False

            # Verify that lists contain valid values
            if not all(isinstance(x, int) and x >= 0 for x in loaded_nheavy):
                LOG.warning("Corrupted state: samples_nheavy contains invalid values")
                return False

            if not all(isinstance(x, (int, float)) and x > 0 for x in loaded_time):
                LOG.warning("Corrupted state: samples_time contains invalid values")
                return False

            # Verify recalibration interval
            loaded_interval = state.get("recalibrate_interval", 50)
            if not isinstance(loaded_interval, int) or loaded_interval <= 0:
                LOG.warning(f"Invalid recalibrate_interval: {loaded_interval}, using default")
                loaded_interval = 50

            # Verify model
            loaded_model = state.get("model")
            if loaded_model is not None and not hasattr(loaded_model, "predict"):
                LOG.warning("Loaded model does not have predict method")
                return False

            # Assign validated values
            self.model = loaded_model
            self.samples_nheavy = loaded_nheavy
            self.samples_time = loaded_time
            self.recalibrate_interval = loaded_interval
            self._molecules_since_recalibration = 0

            LOG.info(f"Predictor loaded from {filepath} ({self.n_samples} samples)")
            return True
        except Exception as e:
            LOG.warning(f"Failed to load predictor: {e}")
            return False

    def get_stats(self) -> dict:
        """Get predictor statistics.

        Returns:
            Dictionary with stats
        """
        stats = {
            "n_samples": self.n_samples,
            "is_fitted": self.is_fitted,
            "molecules_since_recalibration": self._molecules_since_recalibration,
        }

        # Verificar modelo e extrair coeficientes de forma segura
        if self.is_fitted:
            coef = getattr(self.model, "coef_", None)
            intercept = getattr(self.model, "intercept_", None)

            if coef is not None and hasattr(coef, "__len__") and len(coef) > 0:
                try:
                    stats["coef"] = float(coef[0])
                except (IndexError, TypeError):
                    pass

            if intercept is not None:
                try:
                    stats["intercept"] = float(intercept)
                except (TypeError, ValueError):
                    pass

        if self.n_samples > 0:
            stats["nheavy_range"] = (min(self.samples_nheavy), max(self.samples_nheavy))
            stats["time_range"] = (min(self.samples_time), max(self.samples_time))

        return stats
