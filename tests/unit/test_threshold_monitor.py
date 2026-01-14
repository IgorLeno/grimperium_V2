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
        assert metrics.total_successful == 0
        assert metrics.hof_extractions == 0
        assert metrics.grades_a == 0
        assert metrics.grades_b == 0
        assert metrics.grades_c == 0

    def test_success_rate_zero_total(self) -> None:
        """Test success_rate when no molecules processed."""
        metrics = MonitoringMetrics()
        assert metrics.success_rate == 0.0

    def test_success_rate_calculation(self) -> None:
        """Test success_rate calculation."""
        metrics = MonitoringMetrics()
        metrics.total_processed = 10
        metrics.total_successful = 8
        assert metrics.success_rate == pytest.approx(80.0)

    def test_hof_extraction_rate_calculation(self) -> None:
        """Test hof_extraction_rate calculation."""
        metrics = MonitoringMetrics()
        metrics.total_successful = 10
        metrics.hof_extractions = 9
        assert metrics.hof_extraction_rate == pytest.approx(90.0)

    def test_grade_ab_rate_calculation(self) -> None:
        """Test grade_ab_rate calculation."""
        metrics = MonitoringMetrics()
        metrics.total_successful = 10
        metrics.grades_a = 5
        metrics.grades_b = 3
        metrics.grades_c = 2
        assert metrics.grade_ab_rate == pytest.approx(80.0)


class TestAlert:
    """Tests for Alert dataclass."""

    def test_init(self) -> None:
        """Test Alert initialization."""
        alert = Alert(
            level=AlertLevel.WARNING,
            code="TEST001",
            message="Test alert",
        )
        assert alert.level == AlertLevel.WARNING
        assert alert.code == "TEST001"
        assert alert.message == "Test alert"
        assert alert.timestamp is not None


class TestThresholdMonitor:
    """Tests for ThresholdMonitor class."""

    @pytest.fixture
    def config(self) -> PM7Config:
        """Create a test configuration."""
        return PM7Config(
            monitor_window_size=10,
            success_rate_threshold=0.8,
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
        assert monitor.metrics.total_successful == 1
        assert monitor.metrics.hof_extractions == 1
        assert monitor.metrics.grades_a == 1

    def test_record_result_failure(self, monitor: ThresholdMonitor) -> None:
        """Test recording a failed result."""
        monitor.record_result(
            success=False,
            hof_extracted=False,
            grade=QualityGrade.C,
        )

        assert monitor.metrics.total_processed == 1
        assert monitor.metrics.total_successful == 0
        assert monitor.metrics.hof_extractions == 0
        assert monitor.metrics.grades_c == 1

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
        assert monitor.metrics.total_successful == 5
        assert monitor.metrics.success_rate == pytest.approx(50.0)

    def test_alert_callback(self, monitor: ThresholdMonitor) -> None:
        """Test alert callback registration and triggering."""
        alerts_received: list[Alert] = []

        def callback(alert: Alert) -> None:
            alerts_received.append(alert)

        monitor.register_callback(callback)

        # Force an alert by recording many consecutive failures
        for _ in range(6):
            monitor.record_result(
                success=False,
                hof_extracted=False,
                grade=QualityGrade.C,
            )

        # Should have received at least one alert
        assert len(alerts_received) > 0

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
        assert "alerts_count" in summary
        assert summary["total_processed"] == 1
        assert summary["success_rate"] == pytest.approx(100.0)

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

        assert monitor.metrics.grades_a == 2
        assert monitor.metrics.grades_b == 1
        assert monitor.metrics.grades_c == 1
        # A+B rate should be 75%
        assert monitor.metrics.grade_ab_rate == pytest.approx(75.0)
