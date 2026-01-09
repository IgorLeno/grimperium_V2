"""Timeout prediction for MOPAC calculations.

Uses Huber regression to predict optimal timeouts based on molecule complexity.
"""

import logging
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from sklearn.linear_model import HuberRegressor
import numpy as np

from .config import TimeoutConfidence

LOG = logging.getLogger("grimperium.crest_pm7.timeout_predictor")

# Constants
DEFAULT_TIMEOUT = 300.0  # 5 minutes
MIN_SAMPLES_FOR_FIT = 20  # Minimum samples to train model
TIMEOUT_MIN = 60.0  # Minimum timeout (1 minute)
TIMEOUT_MAX = 3600.0  # Maximum timeout (1 hour)


def _clamp(value: float, min_val: float = TIMEOUT_MIN, max_val: float = TIMEOUT_MAX) -> float:
    """Clamp value to [min, max] range."""
    return max(min_val, min(max_val, value))


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
    """

    model: Optional[HuberRegressor] = None
    samples_nheavy: list[int] = field(default_factory=list)
    samples_time: list[float] = field(default_factory=list)
    recalibrate_interval: int = 50
    _molecules_since_recalibration: int = 0

    @property
    def n_samples(self) -> int:
        """Number of training samples."""
        return len(self.samples_nheavy)

    @property
    def is_fitted(self) -> bool:
        """Whether the model has been fitted."""
        return self.model is not None

    def add_observation(self, nheavy: int, execution_time: float) -> None:
        """Add a new observation for model training.

        Args:
            nheavy: Number of heavy atoms
            execution_time: Actual execution time in seconds
        """
        self.samples_nheavy.append(nheavy)
        self.samples_time.append(execution_time)
        self._molecules_since_recalibration += 1

        LOG.debug(f"Added observation: nheavy={nheavy}, time={execution_time:.1f}s")

        # Check if recalibration needed
        if self._molecules_since_recalibration >= self.recalibrate_interval:
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

            self.model = HuberRegressor()
            self.model.fit(X, y)
            self._molecules_since_recalibration = 0

            LOG.info(f"Model fitted with {self.n_samples} samples")
            return True

        except Exception as e:
            LOG.warning(f"Model fitting failed: {e}")
            return False

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
            pred = self.model.predict([[nheavy]])[0]
            timeout = pred * num_conformers * 1.5
            return _clamp(timeout), TimeoutConfidence.MEDIUM

        else:
            # Phase C+: Model + 30% margin
            pred = self.model.predict([[nheavy]])[0]
            timeout = pred * num_conformers * 1.3
            return _clamp(timeout), TimeoutConfidence.HIGH

    def save(self, filepath: Path) -> bool:
        """Save predictor state to file.

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
            with open(filepath, "wb") as f:
                pickle.dump(state, f)
            LOG.info(f"Predictor saved to {filepath}")
            return True
        except Exception as e:
            LOG.warning(f"Failed to save predictor: {e}")
            return False

    def load(self, filepath: Path) -> bool:
        """Load predictor state from file.

        Args:
            filepath: Path to saved file

        Returns:
            True if load succeeded
        """
        try:
            with open(filepath, "rb") as f:
                state = pickle.load(f)
            self.model = state.get("model")
            self.samples_nheavy = state.get("samples_nheavy", [])
            self.samples_time = state.get("samples_time", [])
            self.recalibrate_interval = state.get("recalibrate_interval", 50)
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

        if self.is_fitted and hasattr(self.model, "coef_"):
            stats["coef"] = float(self.model.coef_[0])
            stats["intercept"] = float(self.model.intercept_)

        if self.n_samples > 0:
            stats["nheavy_range"] = (min(self.samples_nheavy), max(self.samples_nheavy))
            stats["time_range"] = (min(self.samples_time), max(self.samples_time))

        return stats
