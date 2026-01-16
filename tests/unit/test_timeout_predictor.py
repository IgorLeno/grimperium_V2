"""Unit tests for TimeoutPredictor.

Tests timeout prediction using Huber regression.
Requires rdkit to be installed due to crest_pm7 package import chain.
"""

import tempfile
from pathlib import Path

import pytest

# Skip all tests if rdkit is not available (needed due to crest_pm7 import chain)
pytest.importorskip("rdkit", reason="rdkit not available")

from grimperium.crest_pm7.config import TimeoutConfidence
from grimperium.crest_pm7.timeout_predictor import (
    MIN_SAMPLES_FOR_FIT,
    TIMEOUT_MAX,
    TIMEOUT_MIN,
    TimeoutPredictor,
    _clamp,
)


class TestClampFunction:
    """Tests for _clamp helper function."""

    def test_clamp_below_min(self) -> None:
        """Test clamping below minimum."""
        assert _clamp(10.0, 60.0, 3600.0) == 60.0

    def test_clamp_above_max(self) -> None:
        """Test clamping above maximum."""
        assert _clamp(5000.0, 60.0, 3600.0) == 3600.0

    def test_clamp_in_range(self) -> None:
        """Test value in range stays unchanged."""
        assert _clamp(300.0, 60.0, 3600.0) == 300.0

    def test_clamp_defaults(self) -> None:
        """Test with default min/max."""
        assert _clamp(10.0) == TIMEOUT_MIN
        assert _clamp(10000.0) == TIMEOUT_MAX
        assert _clamp(500.0) == 500.0


class TestTimeoutPredictor:
    """Tests for TimeoutPredictor dataclass."""

    def test_init_defaults(self) -> None:
        """Test default initialization."""
        predictor = TimeoutPredictor()
        assert predictor.model is None
        assert predictor.samples_nheavy == []
        assert predictor.samples_time == []
        assert predictor.n_samples == 0
        assert not predictor.is_fitted

    def test_add_observation_valid(self) -> None:
        """Test adding valid observation."""
        predictor = TimeoutPredictor()
        predictor.add_observation(10, 120.5)

        assert predictor.n_samples == 1
        assert predictor.samples_nheavy == [10]
        assert predictor.samples_time == [120.5]

    def test_add_observation_invalid_nheavy(self) -> None:
        """Test adding observation with invalid nheavy."""
        predictor = TimeoutPredictor()

        with pytest.raises(ValueError, match="nheavy must be a non-negative integer"):
            predictor.add_observation(-1, 100.0)

        with pytest.raises(ValueError, match="nheavy must be a non-negative integer"):
            predictor.add_observation(10.5, 100.0)  # type: ignore

    def test_add_observation_invalid_time(self) -> None:
        """Test adding observation with invalid execution_time."""
        predictor = TimeoutPredictor()

        with pytest.raises(
            ValueError, match="execution_time must be a positive number"
        ):
            predictor.add_observation(10, 0)

        with pytest.raises(
            ValueError, match="execution_time must be a positive number"
        ):
            predictor.add_observation(10, -50)

    def test_add_multiple_observations(self) -> None:
        """Test adding multiple observations."""
        predictor = TimeoutPredictor()

        for nheavy, time in [(5, 30.0), (10, 120.0), (15, 300.0)]:
            predictor.add_observation(nheavy, time)

        assert predictor.n_samples == 3
        assert predictor.samples_nheavy == [5, 10, 15]
        assert predictor.samples_time == [30.0, 120.0, 300.0]

    def test_fit_insufficient_samples(self) -> None:
        """Test fit with insufficient samples does not train model."""
        predictor = TimeoutPredictor()

        # Add fewer than MIN_SAMPLES_FOR_FIT
        for i in range(MIN_SAMPLES_FOR_FIT - 5):
            predictor.add_observation(i + 5, 100.0 + i * 10)

        predictor.fit()
        assert not predictor.is_fitted

    def test_fit_sufficient_samples(self) -> None:
        """Test fit with sufficient samples trains model."""
        predictor = TimeoutPredictor()

        # Add enough samples
        for i in range(MIN_SAMPLES_FOR_FIT + 5):
            predictor.add_observation(i + 5, 50.0 + i * 20)

        predictor.fit()
        assert predictor.is_fitted
        assert predictor.model is not None

    def test_predict_no_model(self) -> None:
        """Test predict without fitted model uses heuristic."""
        predictor = TimeoutPredictor()

        timeout, confidence = predictor.predict(10)

        # Should return heuristic-based timeout
        assert timeout >= TIMEOUT_MIN
        assert timeout <= TIMEOUT_MAX
        # Confidence should be LOW for heuristic
        assert confidence == TimeoutConfidence.LOW

    def test_predict_with_model(self) -> None:
        """Test predict with fitted model."""
        predictor = TimeoutPredictor()

        # Add synthetic training data
        for i in range(MIN_SAMPLES_FOR_FIT + 10):
            # Simulate linear relationship: time = 10 * nheavy
            predictor.add_observation(i + 5, 10.0 * (i + 5))

        predictor.fit()
        timeout, confidence = predictor.predict(20)

        # Prediction should be close to 200 (10 * 20) + margin
        assert timeout >= TIMEOUT_MIN
        assert timeout <= TIMEOUT_MAX
        # Should be MEDIUM or HIGH confidence with enough data
        assert confidence != TimeoutConfidence.LOW

    def test_predict_clamping(self) -> None:
        """Test that predictions are clamped to valid range."""
        predictor = TimeoutPredictor()

        # Test with very small nheavy
        timeout_small, _ = predictor.predict(1)
        assert timeout_small >= TIMEOUT_MIN

        # Test with very large nheavy
        timeout_large, _ = predictor.predict(200)
        assert timeout_large <= TIMEOUT_MAX

    def test_save_and_load(self) -> None:
        """Test saving and loading predictor state."""
        predictor = TimeoutPredictor()

        # Add data and fit
        for i in range(MIN_SAMPLES_FOR_FIT + 10):
            predictor.add_observation(i + 5, 10.0 * (i + 5))
        predictor.fit()

        # Save to temp file
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "predictor.joblib"
            predictor.save(filepath)

            # Load into new predictor (now instance method, not classmethod)
            loaded = TimeoutPredictor()
            success = loaded.load(filepath)

            assert success is True
            assert loaded.is_fitted
            assert loaded.n_samples == predictor.n_samples

            # Predictions should be the same
            original_pred = predictor.predict(20)
            loaded_pred = loaded.predict(20)
            assert original_pred[0] == pytest.approx(loaded_pred[0])

    def test_get_stats(self) -> None:
        """Test getting predictor statistics."""
        predictor = TimeoutPredictor()

        # Empty predictor
        stats = predictor.get_stats()
        assert stats["n_samples"] == 0
        assert stats["is_fitted"] is False

        # With data
        for i in range(10):
            predictor.add_observation(i + 5, 100.0 + i * 10)

        stats = predictor.get_stats()
        assert stats["n_samples"] == 10
        # API changed: now uses nheavy_range and time_range instead of mean/std
        assert "nheavy_range" in stats
        assert "time_range" in stats
        assert stats["nheavy_range"] == (5, 14)  # min/max nheavy
        assert stats["time_range"] == (100.0, 190.0)  # min/max time

    def test_predictor_phases(self) -> None:
        """Test predictor behavior across phases A, B, C."""
        predictor = TimeoutPredictor()
        # Phase A: < 20 samples - should use heuristic
        for i in range(15):
            predictor.add_observation(i + 5, 50.0 + i * 10)
        predictor.fit()

        _, confidence_a = predictor.predict(10)
        assert confidence_a == TimeoutConfidence.LOW

        # Phase B: 20-50 samples - should use model with medium confidence
        for i in range(15, 35):
            predictor.add_observation(i + 5, 50.0 + i * 10)
        predictor.fit()

        _, confidence_b = predictor.predict(10)
        assert confidence_b == TimeoutConfidence.MEDIUM

        # Phase C: 50+ samples - should use model with high confidence
        for i in range(35, 60):
            predictor.add_observation(i + 5, 50.0 + i * 10)
        predictor.fit()

        _, confidence_c = predictor.predict(10)
        assert confidence_c == TimeoutConfidence.HIGH
