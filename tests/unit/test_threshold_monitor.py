"""Unit tests for ThresholdMonitor.

Tests the monitoring and alerting functionality for batch processing.
"""

import pytest

# Skip all tests in this module if rdkit is not available
pytest.importorskip("rdkit", reason="rdkit not available")

from grimperium.crest_pm7.config import PM7Config, QualityGrade
from grimperium.crest_pm7.threshold_monitor import (
    Alert,
    AlertLevel,
    MonitoringMetrics,
    ThresholdMonitor,
)


class TestMonitoringMetrics:
    """Tests for MonitoringMetrics dataclass."""

    def test_init_defaults(self) -> None:
        """Test default initialization."""
        metrics = MonitoringMetrics()
        assert metrics.total_processed == 0
        assert metrics.success_count == 0  # Changed: total_successful -> success_count
        assert metrics.hof_extraction_success == 0  # Changed: hof_extractions -> hof_extraction_success
        assert metrics.grade_a_count == 0  # Changed: grades_a -> grade_a_count
        assert metrics.grade_b_count == 0  # Changed: grades_b -> grade_b_count
        assert metrics.grade_c_count == 0  # Changed: grades_c -> grade_c_count

    def test_success_rate_zero_total(self) -> None:
        """Test success_rate when no molecules processed."""
        metrics = MonitoringMetrics()
        assert metrics.success_rate is None  # Changed: 0.0 -> None

    def test_success_rate_calculation(self) -> None:
        """Test success_rate calculation."""
        metrics = MonitoringMetrics()
        metrics.total_processed = 10
        metrics.success_count = 8  # Changed: total_successful -> success_count
        # Returns ratio (0-1), not percentage (0-100)
        assert metrics.success_rate == pytest.approx(0.8)  # Changed: 80.0 -> 0.8

    def test_hof_extraction_rate_calculation(self) -> None:
        """Test hof_extraction_rate calculation."""
        metrics = MonitoringMetrics()
        # Now tracks both success and failure counts
        metrics.hof_extraction_success = 9
        metrics.hof_extraction_failure = 1
        # Returns ratio (0-1), not percentage (0-100)
        assert metrics.hof_extraction_rate == pytest.approx(0.9)  # Changed: 90.0 -> 0.9

    def test_grade_ab_rate_calculation(self) -> None:
        """Test grade_ab_rate calculation."""
        metrics = MonitoringMetrics()
        metrics.total_processed = 10  # Now uses total_processed as denominator
        metrics.grade_a_count = 5  # Changed: grades_a -> grade_a_count
        metrics.grade_b_count = 3  # Changed: grades_b -> grade_b_count
        metrics.grade_c_count = 2  # Changed: grades_c -> grade_c_count
        # Returns ratio (0-1), not percentage (0-100)
        assert metrics.grade_ab_rate == pytest.approx(0.8)  # Changed: 80.0 -> 0.8


class TestAlert:
    """Tests for Alert dataclass."""

    def test_init(self) -> None:
        """Test Alert initialization."""
        alert = Alert(
            level=AlertLevel.WARNING,
            pattern="TEST001",  # Changed: code -> pattern
            message="Test alert",
        )
        assert alert.level == AlertLevel.WARNING
        assert alert.pattern == "TEST001"  # Changed: code -> pattern
        assert alert.message == "Test alert"
        assert alert.timestamp is not None


class TestThresholdMonitor:
    """Tests for ThresholdMonitor class."""

    @pytest.fixture
    def config(self) -> PM7Config:
        """Create a test configuration."""
        # Use existing PM7Config fields
        return PM7Config(
            monitor_window_size=10,
            success_rate_warning=0.8,  # Changed: success_rate_threshold -> success_rate_warning
            hof_extraction_threshold=0.9,
        )

    @pytest.fixture
    def monitor(self, config: PM7Config) -> ThresholdMonitor:
        """Create a test monitor."""
        return ThresholdMonitor(config)

    def test_init(self, monitor: ThresholdMonitor) -> None:
        """Test monitor initialization."""
        assert monitor.metrics.total_processed == 0
        assert len(monitor.alerts) == 0
        assert monitor.window_size == 10

    def test_record_result_success(self, monitor: ThresholdMonitor) -> None:
        """Test recording a successful result."""
        monitor.record_result(
            success=True,
            hof_extracted=True,
            grade=QualityGrade.A,
        )

        assert monitor.metrics.total_processed == 1
        assert monitor.metrics.success_count == 1  # Changed: total_successful -> success_count
        assert monitor.metrics.hof_extraction_success == 1  # Changed: hof_extractions -> hof_extraction_success
        assert monitor.metrics.grade_a_count == 1  # Changed: grades_a -> grade_a_count

    def test_record_result_failure(self, monitor: ThresholdMonitor) -> None:
        """Test recording a failed result."""
        monitor.record_result(
            success=False,
            hof_extracted=False,
            grade=QualityGrade.C,
        )

        assert monitor.metrics.total_processed == 1
        assert monitor.metrics.success_count == 0  # Changed: total_successful -> success_count
        assert monitor.metrics.hof_extraction_failure == 1  # Changed: hof_extractions == 0 -> hof_extraction_failure == 1
        assert monitor.metrics.grade_c_count == 1  # Changed: grades_c -> grade_c_count

    def test_record_multiple_results(self, monitor: ThresholdMonitor) -> None:
        """Test recording multiple results."""
        # 5 successes, 5 failures
        for _ in range(5):
            monitor.record_result(
                success=True,
                hof_extracted=True,
                grade=QualityGrade.A,
            )
        for _ in range(5):
            monitor.record_result(
                success=False,
                hof_extracted=False,
                grade=QualityGrade.C,
            )

        assert monitor.metrics.total_processed == 10
        assert monitor.metrics.success_count == 5  # Changed: total_successful -> success_count
        # Returns ratio (0-1), not percentage (0-100)
        assert monitor.metrics.success_rate == pytest.approx(0.5)  # Changed: 50.0 -> 0.5

    def test_alert_callback(self, monitor: ThresholdMonitor) -> None:
        """Test alert callback registration and triggering."""
        alerts_received: list[Alert] = []

        def callback(alert: Alert) -> None:
            alerts_received.append(alert)

        monitor.register_callback(callback)

        # Force an alert by recording consecutive failures
        # (default consecutive_failure_threshold is typically 5)
        for _ in range(6):
            monitor.record_result(
                success=False,
                hof_extracted=False,
                grade=QualityGrade.C,
            )

        # Should have received at least one alert after exceeding failure threshold
        assert len(alerts_received) >= 1

    def test_get_summary(self, monitor: ThresholdMonitor) -> None:
        """Test getting summary."""
        monitor.record_result(
            success=True,
            hof_extracted=True,
            grade=QualityGrade.A,
        )

        summary = monitor.get_summary()
        assert "total_processed" in summary
        assert "success_rate" in summary
        assert "alert_count" in summary  # Changed: alerts_count -> alert_count
        assert summary["total_processed"] == 1
        # Returns ratio (0-1), not percentage (0-100)
        assert summary["success_rate"] == pytest.approx(1.0)  # Changed: 100.0 -> 1.0

    def test_should_pause_with_no_issues(self, monitor: ThresholdMonitor) -> None:
        """Test should_pause returns False when no issues."""
        # Record some successful results
        for _ in range(5):
            monitor.record_result(
                success=True,
                hof_extracted=True,
                grade=QualityGrade.A,
            )

        assert not monitor.should_pause()

    def test_grade_distribution(self, monitor: ThresholdMonitor) -> None:
        """Test grade distribution tracking."""
        monitor.record_result(success=True, hof_extracted=True, grade=QualityGrade.A)
        monitor.record_result(success=True, hof_extracted=True, grade=QualityGrade.A)
        monitor.record_result(success=True, hof_extracted=True, grade=QualityGrade.B)
        monitor.record_result(success=True, hof_extracted=True, grade=QualityGrade.C)

        assert monitor.metrics.grade_a_count == 2  # Changed: grades_a -> grade_a_count
        assert monitor.metrics.grade_b_count == 1  # Changed: grades_b -> grade_b_count
        assert monitor.metrics.grade_c_count == 1  # Changed: grades_c -> grade_c_count
        # A+B rate should be 75% -> 0.75 (ratio)
        assert monitor.metrics.grade_ab_rate == pytest.approx(0.75)  # Changed: 75.0 -> 0.75
