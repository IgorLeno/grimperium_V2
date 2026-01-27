"""
Tests for ProgressTracker with CSV-driven state machine and Queue pattern.

Coverage requirements: >= 85%

Test categories:
    - ProcessingStage enum values and ordering
    - StageTransition and EVENTS tuple (5 elements, immutable)
    - MoleculeProgress state tracking with _detect_stage()
    - ProgressTracker register, update, render (30 chars exactly)
    - CSVMonitor daemon thread with Queue pattern
    - Thread safety and concurrent operations
    - Batch header statistics
"""

from __future__ import annotations

import threading
import time
from pathlib import Path
from queue import Queue
from typing import TYPE_CHECKING
from unittest.mock import MagicMock

import pandas as pd
import pytest
from rich.console import Console

from grimperium.cli.progress_tracker import (
    EVENTS,
    CSVMonitor,
    MoleculeProgress,
    ProcessingStage,
    ProgressEvent,
    ProgressTracker,
)

if TYPE_CHECKING:
    from queue import Queue as QueueType


# ============== FIXTURES ==============


@pytest.fixture
def mock_console() -> MagicMock:
    """Mock Rich Console with standard width."""
    console = MagicMock(spec=Console)
    console.width = 80
    return console


@pytest.fixture
def progress_tracker(mock_console: MagicMock) -> ProgressTracker:
    """ProgressTracker with mocked console and batch_size=6."""
    return ProgressTracker(console=mock_console, batch_size=6)


@pytest.fixture
def event_queue() -> QueueType[ProgressEvent]:
    """Thread-safe event queue for CSVMonitor."""
    return Queue()


@pytest.fixture
def sample_csv_all_states(tmp_path: Path) -> Path:
    """CSV with molecules in ALL 5 processing states.

    States covered:
        - mol_001: pending (stage 0)
        - mol_002: processing, NOT_ATTEMPTED (stage 1 - RDKit)
        - mol_003: processing, xtb_opt (stage 2 - xTB)
        - mol_004: processing, conformer_search (stage 3 - CREST)
        - mol_005: processing, geometric_opt (stage 4 - MOPAC)
        - mol_006: OK, SUCCESS, success (stage 5 - Complete)
    """
    csv_content = (
        "mol_id,smiles,status,crest_status,mopac_status\n"
        "mol_001,C,pending,NOT_ATTEMPTED,none\n"
        "mol_002,CC,processing,NOT_ATTEMPTED,none\n"
        "mol_003,CCC,processing,xtb_opt,none\n"
        "mol_004,CCCC,processing,conformer_search,none\n"
        "mol_005,CCCCC,processing,conformer_search,geometric_opt\n"
        "mol_006,CCCCCC,OK,SUCCESS,success\n"
    )
    csv_path = tmp_path / "test_batch.csv"
    csv_path.write_text(csv_content)
    return csv_path


@pytest.fixture
def csv_monitor(
    sample_csv_all_states: Path,
    event_queue: QueueType[ProgressEvent],
) -> CSVMonitor:
    """CSVMonitor with test CSV and faster poll interval for tests."""
    monitor = CSVMonitor(
        csv_path=sample_csv_all_states,
        event_queue=event_queue,
        poll_interval_ms=100,  # Faster for tests
    )
    for i in range(1, 7):
        monitor.register_molecule(f"mol_00{i}")
    return monitor


# ============== TEST CLASSES ==============


class TestProcessingStage:
    """Test ProcessingStage enum values and ordering."""

    def test_stage_values(self) -> None:
        """Verify each stage has correct integer value."""
        assert ProcessingStage.NOT_STARTED == 0
        assert ProcessingStage.RDKIT_PARAMS == 1
        assert ProcessingStage.XTB_PREOPT == 2
        assert ProcessingStage.CREST_SEARCH == 3
        assert ProcessingStage.MOPAC_CALC == 4
        assert ProcessingStage.FINAL_CALC == 5

    def test_stage_count(self) -> None:
        """Verify total number of stages (including NOT_STARTED)."""
        assert len(ProcessingStage) == 6

    def test_stage_ordering(self) -> None:
        """Verify stages are ordered correctly for comparison."""
        assert ProcessingStage.NOT_STARTED < ProcessingStage.RDKIT_PARAMS
        assert ProcessingStage.RDKIT_PARAMS < ProcessingStage.XTB_PREOPT
        assert ProcessingStage.XTB_PREOPT < ProcessingStage.CREST_SEARCH
        assert ProcessingStage.CREST_SEARCH < ProcessingStage.MOPAC_CALC
        assert ProcessingStage.MOPAC_CALC < ProcessingStage.FINAL_CALC


class TestStageTransition:
    """Test EVENTS tuple definitions."""

    def test_events_count(self) -> None:
        """EVENTS tuple must have exactly 5 transitions."""
        assert len(EVENTS) == 5

    def test_events_are_immutable(self) -> None:
        """EVENTS tuple should be immutable (frozen)."""
        with pytest.raises(TypeError):
            EVENTS[0] = None  # type: ignore[index]

    def test_event_dataclass_is_frozen(self) -> None:
        """StageTransition instances should be immutable."""
        event = EVENTS[0]
        with pytest.raises(AttributeError):
            event.column = "other"  # type: ignore[misc]

    def test_first_event_is_status_pending_to_processing(self) -> None:
        """First event: status pending -> processing (RDKit)."""
        event = EVENTS[0]
        assert event.column == "status"
        assert event.from_value == "pending"
        assert event.to_value == "processing"
        assert event.stage == ProcessingStage.RDKIT_PARAMS
        assert event.label == "RDKit parameters"

    def test_second_event_is_crest_status_to_xtb_opt(self) -> None:
        """Second event: crest_status NOT_ATTEMPTED -> xtb_opt (xTB)."""
        event = EVENTS[1]
        assert event.column == "crest_status"
        assert event.from_value == "NOT_ATTEMPTED"
        assert event.to_value == "xtb_opt"
        assert event.stage == ProcessingStage.XTB_PREOPT
        assert event.label == "xTB pre-optimization"

    def test_third_event_is_crest_status_to_conformer_search(self) -> None:
        """Third event: crest_status xtb_opt -> conformer_search (CREST)."""
        event = EVENTS[2]
        assert event.column == "crest_status"
        assert event.from_value == "xtb_opt"
        assert event.to_value == "conformer_search"
        assert event.stage == ProcessingStage.CREST_SEARCH
        assert event.label == "CREST conformer search"

    def test_fourth_event_is_mopac_status_to_geometric_opt(self) -> None:
        """Fourth event: mopac_status none -> geometric_opt (MOPAC)."""
        event = EVENTS[3]
        assert event.column == "mopac_status"
        assert event.from_value == "none"
        assert event.to_value == "geometric_opt"
        assert event.stage == ProcessingStage.MOPAC_CALC
        assert event.label == "MOPAC PM7 calculation"

    def test_last_event_is_status_processing_to_ok(self) -> None:
        """Fifth event: status processing -> OK (Final)."""
        event = EVENTS[4]
        assert event.column == "status"
        assert event.from_value == "processing"
        assert event.to_value == "OK"
        assert event.stage == ProcessingStage.FINAL_CALC
        assert event.label == "Final calculations"


class TestMoleculeProgress:
    """Test MoleculeProgress state tracking with _detect_stage()."""

    def test_initial_state(self) -> None:
        """New MoleculeProgress starts at NOT_STARTED with correct defaults."""
        progress = MoleculeProgress(mol_id="test")
        assert progress.current_stage == ProcessingStage.NOT_STARTED
        assert progress.last_csv_state == {
            "status": "pending",
            "crest_status": "NOT_ATTEMPTED",
            "mopac_status": "none",
        }
        assert progress.completed is False
        assert progress.error is None

    def test_detect_stage_pending_to_processing(self) -> None:
        """Detect transition from pending to processing (stage 1)."""
        progress = MoleculeProgress(mol_id="test")
        new_stage = progress._detect_stage(
            {
                "status": "processing",
                "crest_status": "NOT_ATTEMPTED",
                "mopac_status": "none",
            }
        )
        assert new_stage == ProcessingStage.RDKIT_PARAMS

    def test_detect_stage_xtb_opt(self) -> None:
        """Detect transition to xtb_opt (stage 2)."""
        progress = MoleculeProgress(mol_id="test")
        progress.last_csv_state = {
            "status": "processing",
            "crest_status": "NOT_ATTEMPTED",
            "mopac_status": "none",
        }
        new_stage = progress._detect_stage(
            {
                "status": "processing",
                "crest_status": "xtb_opt",
                "mopac_status": "none",
            }
        )
        assert new_stage == ProcessingStage.XTB_PREOPT

    def test_detect_stage_conformer_search(self) -> None:
        """Detect transition to conformer_search (stage 3)."""
        progress = MoleculeProgress(mol_id="test")
        progress.last_csv_state = {
            "status": "processing",
            "crest_status": "xtb_opt",
            "mopac_status": "none",
        }
        new_stage = progress._detect_stage(
            {
                "status": "processing",
                "crest_status": "conformer_search",
                "mopac_status": "none",
            }
        )
        assert new_stage == ProcessingStage.CREST_SEARCH

    def test_detect_stage_geometric_opt(self) -> None:
        """Detect transition to geometric_opt (stage 4)."""
        progress = MoleculeProgress(mol_id="test")
        progress.last_csv_state = {
            "status": "processing",
            "crest_status": "conformer_search",
            "mopac_status": "none",
        }
        new_stage = progress._detect_stage(
            {
                "status": "processing",
                "crest_status": "conformer_search",
                "mopac_status": "geometric_opt",
            }
        )
        assert new_stage == ProcessingStage.MOPAC_CALC

    def test_detect_stage_ok(self) -> None:
        """Detect transition to OK (stage 5)."""
        progress = MoleculeProgress(mol_id="test")
        progress.last_csv_state = {
            "status": "processing",
            "crest_status": "SUCCESS",
            "mopac_status": "success",
        }
        new_stage = progress._detect_stage(
            {
                "status": "OK",
                "crest_status": "SUCCESS",
                "mopac_status": "success",
            }
        )
        assert new_stage == ProcessingStage.FINAL_CALC

    def test_detect_stage_no_change(self) -> None:
        """No stage change when CSV state is unchanged."""
        progress = MoleculeProgress(mol_id="test")
        new_stage = progress._detect_stage(
            {
                "status": "pending",
                "crest_status": "NOT_ATTEMPTED",
                "mopac_status": "none",
            }
        )
        assert new_stage is None

    def test_update_from_csv_row_updates_last_state(self) -> None:
        """update_from_csv_row() must update last_csv_state."""
        progress = MoleculeProgress(mol_id="test")
        row = pd.Series(
            {
                "status": "processing",
                "crest_status": "xtb_opt",
                "mopac_status": "none",
            }
        )
        progress.update_from_csv_row(row)
        assert progress.last_csv_state["status"] == "processing"
        assert progress.last_csv_state["crest_status"] == "xtb_opt"
        assert progress.last_csv_state["mopac_status"] == "none"

    def test_update_from_csv_row_returns_new_stage(self) -> None:
        """update_from_csv_row() returns new stage when transition detected."""
        progress = MoleculeProgress(mol_id="test")
        row = pd.Series(
            {
                "status": "processing",
                "crest_status": "NOT_ATTEMPTED",
                "mopac_status": "none",
            }
        )
        new_stage = progress.update_from_csv_row(row)
        assert new_stage == ProcessingStage.RDKIT_PARAMS

    def test_update_from_csv_row_updates_current_stage(self) -> None:
        """update_from_csv_row() updates current_stage on transition."""
        progress = MoleculeProgress(mol_id="test")
        row = pd.Series(
            {
                "status": "processing",
                "crest_status": "NOT_ATTEMPTED",
                "mopac_status": "none",
            }
        )
        progress.update_from_csv_row(row)
        assert progress.current_stage == ProcessingStage.RDKIT_PARAMS

    def test_completed_on_ok_status(self) -> None:
        """Molecule is marked completed when status is OK."""
        progress = MoleculeProgress(mol_id="test")
        row = pd.Series(
            {"status": "OK", "crest_status": "SUCCESS", "mopac_status": "success"}
        )
        progress.update_from_csv_row(row)
        assert progress.completed is True

    def test_completed_on_skip_status(self) -> None:
        """Molecule is marked completed when status is Skip."""
        progress = MoleculeProgress(mol_id="test")
        row = pd.Series(
            {"status": "Skip", "crest_status": "ERROR", "mopac_status": "none"}
        )
        progress.update_from_csv_row(row)
        assert progress.completed is True

    def test_completed_on_rerun_status(self) -> None:
        """Molecule is marked completed when status is Rerun."""
        progress = MoleculeProgress(mol_id="test")
        row = pd.Series(
            {"status": "Rerun", "crest_status": "timeout", "mopac_status": "none"}
        )
        progress.update_from_csv_row(row)
        assert progress.completed is True


class TestProgressTracker:
    """Test ProgressTracker class methods."""

    def test_tracker_initialization(self, mock_console: MagicMock) -> None:
        """ProgressTracker initializes with correct defaults."""
        tracker = ProgressTracker(console=mock_console, batch_size=10)
        assert tracker.batch_size == 10
        assert tracker.total_processed == 0
        assert tracker.successful == 0
        assert tracker.failed == 0
        assert tracker.skipped == 0
        assert len(tracker._molecules) == 0

    def test_register_molecule(self, progress_tracker: ProgressTracker) -> None:
        """register_molecule() adds molecule to tracking dict."""
        progress_tracker.register_molecule("mol_001")
        assert "mol_001" in progress_tracker._molecules
        assert isinstance(progress_tracker._molecules["mol_001"], MoleculeProgress)

    def test_progress_bar_width_exactly_30_chars(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """Progress bar must be exactly 30 characters."""
        progress_tracker.register_molecule("mol_001")
        bar = progress_tracker.render_progress_bar("mol_001")
        assert len(bar) == 30

    def test_progress_bar_fills_6_chars_per_stage(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """Each stage fills exactly 6 characters in the progress bar."""
        progress_tracker.register_molecule("mol_001")

        # Stage 0: all empty
        bar = progress_tracker.render_progress_bar("mol_001")
        assert bar.count(ProgressTracker.FILLED_CHAR) == 0
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 30

        # Apply stage 1 event
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.RDKIT_PARAMS)
        )
        bar = progress_tracker.render_progress_bar("mol_001")
        assert bar.count(ProgressTracker.FILLED_CHAR) == 6
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 24

        # Apply stage 2 event
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.XTB_PREOPT)
        )
        bar = progress_tracker.render_progress_bar("mol_001")
        assert bar.count(ProgressTracker.FILLED_CHAR) == 12
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 18

        # Apply stage 3 event
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.CREST_SEARCH)
        )
        bar = progress_tracker.render_progress_bar("mol_001")
        assert bar.count(ProgressTracker.FILLED_CHAR) == 18
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 12

        # Apply stage 4 event
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.MOPAC_CALC)
        )
        bar = progress_tracker.render_progress_bar("mol_001")
        assert bar.count(ProgressTracker.FILLED_CHAR) == 24
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 6

        # Apply stage 5 event
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.FINAL_CALC)
        )
        bar = progress_tracker.render_progress_bar("mol_001")
        assert bar.count(ProgressTracker.FILLED_CHAR) == 30
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 0

    def test_apply_event_updates_stage(self, progress_tracker: ProgressTracker) -> None:
        """apply_event() updates molecule's current_stage."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.MOPAC_CALC)
        )
        assert (
            progress_tracker._molecules["mol_001"].current_stage
            == ProcessingStage.MOPAC_CALC
        )

    def test_unknown_event_ignored(self, progress_tracker: ProgressTracker) -> None:
        """apply_event() for unknown molecule does not raise."""
        # Should not raise
        progress_tracker.apply_event(
            ProgressEvent("unknown_mol", ProcessingStage.RDKIT_PARAMS)
        )
        assert "unknown_mol" not in progress_tracker._molecules

    def test_get_stage_label_pending(self, progress_tracker: ProgressTracker) -> None:
        """get_stage_label() returns 'Pending' for NOT_STARTED."""
        progress_tracker.register_molecule("mol_001")
        label = progress_tracker.get_stage_label("mol_001")
        assert label == "Pending"

    def test_get_stage_label_rdkit(self, progress_tracker: ProgressTracker) -> None:
        """get_stage_label() returns correct label for stage 1."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.RDKIT_PARAMS)
        )
        label = progress_tracker.get_stage_label("mol_001")
        assert label == "RDKit parameters"

    def test_get_stage_label_crest(self, progress_tracker: ProgressTracker) -> None:
        """get_stage_label() returns correct label for stage 3."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.CREST_SEARCH)
        )
        label = progress_tracker.get_stage_label("mol_001")
        assert label == "CREST conformer search"

    def test_render_progress_bar_unknown_molecule(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """render_progress_bar() for unknown molecule returns empty bar."""
        bar = progress_tracker.render_progress_bar("unknown")
        assert len(bar) == 30
        assert bar == ProgressTracker.EMPTY_CHAR * 30


class TestProgressTrackerStats:
    """Test batch header statistics tracking."""

    def test_initial_stats(self, progress_tracker: ProgressTracker) -> None:
        """Initial statistics are all zero."""
        assert progress_tracker.total_processed == 0
        assert progress_tracker.successful == 0
        assert progress_tracker.failed == 0
        assert progress_tracker.skipped == 0

    def test_mark_completed_success_updates_stats(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """mark_completed() with success=True updates successful count."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.mark_completed("mol_001", success=True)

        assert progress_tracker.total_processed == 1
        assert progress_tracker.successful == 1
        assert progress_tracker.failed == 0

    def test_mark_completed_failure_updates_stats(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """mark_completed() with success=False updates failed count."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.mark_completed("mol_001", success=False)

        assert progress_tracker.total_processed == 1
        assert progress_tracker.successful == 0
        assert progress_tracker.failed == 1

    def test_mark_completed_removes_molecule(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """mark_completed() removes molecule from active tracking."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.mark_completed("mol_001", success=True)

        assert "mol_001" not in progress_tracker._molecules

    def test_multiple_completions(self, progress_tracker: ProgressTracker) -> None:
        """Track multiple molecule completions correctly."""
        for i in range(1, 6):
            progress_tracker.register_molecule(f"mol_00{i}")

        progress_tracker.mark_completed("mol_001", success=True)
        progress_tracker.mark_completed("mol_002", success=True)
        progress_tracker.mark_completed("mol_003", success=False)
        progress_tracker.mark_completed("mol_004", success=True)
        progress_tracker.mark_completed("mol_005", success=False)

        assert progress_tracker.total_processed == 5
        assert progress_tracker.successful == 3
        assert progress_tracker.failed == 2


class TestCSVMonitor:
    """Test CSV polling daemon thread."""

    def test_daemon_thread(self, csv_monitor: CSVMonitor) -> None:
        """CSVMonitor thread must be a daemon thread."""
        csv_monitor.start()
        assert csv_monitor._thread is not None
        assert csv_monitor._thread.daemon is True
        csv_monitor.stop()

    def test_poll_interval_default(self) -> None:
        """Default poll interval is 500ms."""
        assert CSVMonitor.DEFAULT_POLL_INTERVAL_MS == 500

    def test_stop_gracefully(self, csv_monitor: CSVMonitor) -> None:
        """stop() terminates thread gracefully."""
        csv_monitor.start()
        csv_monitor.stop(timeout=1.0)
        assert csv_monitor._thread is None or not csv_monitor._thread.is_alive()

    def test_start_twice_is_safe(self, csv_monitor: CSVMonitor) -> None:
        """Calling start() twice should not create duplicate threads."""
        csv_monitor.start()
        first_thread = csv_monitor._thread
        csv_monitor.start()
        assert csv_monitor._thread is first_thread
        csv_monitor.stop()

    def test_register_molecule(
        self,
        sample_csv_all_states: Path,
        event_queue: QueueType[ProgressEvent],
    ) -> None:
        """register_molecule() adds molecule to monitoring."""
        monitor = CSVMonitor(
            csv_path=sample_csv_all_states,
            event_queue=event_queue,
        )
        monitor.register_molecule("mol_001")
        assert "mol_001" in monitor._molecule_states


class TestCSVMonitorDetection:
    """Test CSV state change detection and event emission."""

    def test_detects_state_change_and_enqueues(
        self,
        tmp_path: Path,
        event_queue: QueueType[ProgressEvent],
    ) -> None:
        """CSVMonitor detects CSV change and enqueues event."""
        # Create CSV with pending molecule
        csv_path = tmp_path / "test.csv"
        csv_path.write_text(
            "mol_id,status,crest_status,mopac_status\n"
            "mol_001,pending,NOT_ATTEMPTED,none"
        )

        monitor = CSVMonitor(
            csv_path=csv_path,
            event_queue=event_queue,
            poll_interval_ms=50,
        )
        monitor.register_molecule("mol_001")
        monitor.start()

        # Update CSV to trigger transition
        csv_path.write_text(
            "mol_id,status,crest_status,mopac_status\n"
            "mol_001,processing,NOT_ATTEMPTED,none"
        )

        # Wait for event with timeout (deterministic synchronization)
        try:
            event = event_queue.get(timeout=1.0)
        finally:
            monitor.stop()

        # Verify event
        assert event.mol_id == "mol_001"
        assert event.new_stage == ProcessingStage.RDKIT_PARAMS

    def test_no_event_when_no_change(
        self,
        tmp_path: Path,
        event_queue: QueueType[ProgressEvent],
    ) -> None:
        """No event emitted when CSV state is unchanged."""
        csv_path = tmp_path / "test.csv"
        csv_path.write_text(
            "mol_id,status,crest_status,mopac_status\n"
            "mol_001,pending,NOT_ATTEMPTED,none"
        )

        monitor = CSVMonitor(
            csv_path=csv_path,
            event_queue=event_queue,
            poll_interval_ms=50,
        )
        monitor.register_molecule("mol_001")
        monitor.start()

        time.sleep(0.2)  # Let it poll multiple times

        monitor.stop()

        # Queue should be empty (no state change)
        assert event_queue.empty()

    def test_handles_missing_csv(
        self,
        tmp_path: Path,
        event_queue: QueueType[ProgressEvent],
    ) -> None:
        """CSVMonitor handles missing CSV file gracefully."""
        csv_path = tmp_path / "nonexistent.csv"

        monitor = CSVMonitor(
            csv_path=csv_path,
            event_queue=event_queue,
            poll_interval_ms=50,
        )
        monitor.register_molecule("mol_001")
        monitor.start()

        time.sleep(0.15)  # Let it attempt polls

        monitor.stop()

        # Should not crash, queue should be empty
        assert event_queue.empty()

    def test_handles_empty_csv(
        self,
        tmp_path: Path,
        event_queue: QueueType[ProgressEvent],
    ) -> None:
        """CSVMonitor handles empty CSV file gracefully."""
        csv_path = tmp_path / "empty.csv"
        csv_path.write_text("")

        monitor = CSVMonitor(
            csv_path=csv_path,
            event_queue=event_queue,
            poll_interval_ms=50,
        )
        monitor.register_molecule("mol_001")
        monitor.start()

        time.sleep(0.15)

        monitor.stop()

        # Should not crash
        assert event_queue.empty()

    def test_handles_malformed_csv(
        self,
        tmp_path: Path,
        event_queue: QueueType[ProgressEvent],
    ) -> None:
        """CSVMonitor handles malformed CSV gracefully."""
        csv_path = tmp_path / "malformed.csv"
        csv_path.write_text('not,a,valid\ncsv"file')

        monitor = CSVMonitor(
            csv_path=csv_path,
            event_queue=event_queue,
            poll_interval_ms=50,
        )
        monitor.register_molecule("mol_001")
        monitor.start()

        time.sleep(0.15)

        monitor.stop()

        # Should not crash
        assert event_queue.empty()


class TestAllFiveEventsSequence:
    """Test full molecule lifecycle through all 5 stages."""

    def test_all_five_events_in_sequence(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """Molecule progresses through all 5 stages with correct bar fill."""
        progress_tracker.register_molecule("mol_001")

        expected_stages = [
            ProcessingStage.RDKIT_PARAMS,
            ProcessingStage.XTB_PREOPT,
            ProcessingStage.CREST_SEARCH,
            ProcessingStage.MOPAC_CALC,
            ProcessingStage.FINAL_CALC,
        ]

        for i, stage in enumerate(expected_stages):
            progress_tracker.apply_event(ProgressEvent("mol_001", stage))

            bar = progress_tracker.render_progress_bar("mol_001")
            expected_filled = (i + 1) * 6
            assert bar.count(ProgressTracker.FILLED_CHAR) == expected_filled
            assert bar.count(ProgressTracker.EMPTY_CHAR) == 30 - expected_filled

    def test_events_with_correct_labels(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """Each stage has correct label from EVENTS."""
        progress_tracker.register_molecule("mol_001")

        stages_and_labels = [
            (ProcessingStage.RDKIT_PARAMS, "RDKit parameters"),
            (ProcessingStage.XTB_PREOPT, "xTB pre-optimization"),
            (ProcessingStage.CREST_SEARCH, "CREST conformer search"),
            (ProcessingStage.MOPAC_CALC, "MOPAC PM7 calculation"),
            (ProcessingStage.FINAL_CALC, "Final calculations"),
        ]

        for stage, expected_label in stages_and_labels:
            progress_tracker.apply_event(ProgressEvent("mol_001", stage))
            label = progress_tracker.get_stage_label("mol_001")
            assert label == expected_label


class TestMoleculeTransitionResetsProgress:
    """Test that finishing one molecule allows next to start fresh."""

    def test_next_molecule_starts_at_zero(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """New molecule starts with empty progress bar."""
        # First molecule completes
        progress_tracker.register_molecule("mol_001")
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.FINAL_CALC)
        )
        progress_tracker.mark_completed("mol_001", success=True)

        # Second molecule starts fresh
        progress_tracker.register_molecule("mol_002")
        bar = progress_tracker.render_progress_bar("mol_002")

        assert bar.count(ProgressTracker.FILLED_CHAR) == 0
        assert bar.count(ProgressTracker.EMPTY_CHAR) == 30

    def test_completed_molecule_not_in_tracking(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """Completed molecule is removed from active tracking."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.mark_completed("mol_001", success=True)

        assert "mol_001" not in progress_tracker._molecules

    def test_multiple_molecules_independent_progress(
        self, progress_tracker: ProgressTracker
    ) -> None:
        """Multiple molecules track progress independently."""
        progress_tracker.register_molecule("mol_001")
        progress_tracker.register_molecule("mol_002")

        # Advance mol_001 to stage 3
        progress_tracker.apply_event(
            ProgressEvent("mol_001", ProcessingStage.CREST_SEARCH)
        )

        # Advance mol_002 to stage 1
        progress_tracker.apply_event(
            ProgressEvent("mol_002", ProcessingStage.RDKIT_PARAMS)
        )

        bar_001 = progress_tracker.render_progress_bar("mol_001")
        bar_002 = progress_tracker.render_progress_bar("mol_002")

        assert bar_001.count(ProgressTracker.FILLED_CHAR) == 18  # 3 stages
        assert bar_002.count(ProgressTracker.FILLED_CHAR) == 6  # 1 stage


class TestThreadSafety:
    """Test thread safety of Queue-based communication."""

    def test_queue_is_thread_safe(self, event_queue: QueueType[ProgressEvent]) -> None:
        """Queue operations are thread-safe."""
        num_events = 100
        results: list[ProgressEvent] = []

        def producer() -> None:
            for i in range(num_events):
                event_queue.put(ProgressEvent(f"mol_{i}", ProcessingStage.RDKIT_PARAMS))

        def consumer() -> None:
            consumed = 0
            while consumed < num_events:
                try:
                    event = event_queue.get(timeout=1.0)
                    results.append(event)
                    consumed += 1
                except Empty:
                    break

        producer_thread = threading.Thread(target=producer)
        consumer_thread = threading.Thread(target=consumer)

        producer_thread.start()
        consumer_thread.start()

        producer_thread.join()
        consumer_thread.join()

        assert len(results) == num_events

    def test_monitor_and_consumer_concurrent(
        self,
        tmp_path: Path,
        event_queue: QueueType[ProgressEvent],
        mock_console: MagicMock,
    ) -> None:
        """CSVMonitor and main thread can operate concurrently."""
        csv_path = tmp_path / "test.csv"
        csv_path.write_text(
            "mol_id,status,crest_status,mopac_status\n"
            "mol_001,pending,NOT_ATTEMPTED,none"
        )

        tracker = ProgressTracker(console=mock_console, batch_size=1)
        tracker.register_molecule("mol_001")

        monitor = CSVMonitor(
            csv_path=csv_path,
            event_queue=event_queue,
            poll_interval_ms=50,
        )
        monitor.register_molecule("mol_001")
        monitor.start()

        # Simulate main thread consuming events while monitor runs
        time.sleep(0.1)

        # Update CSV
        csv_path.write_text(
            "mol_id,status,crest_status,mopac_status\n"
            "mol_001,processing,NOT_ATTEMPTED,none"
        )

        time.sleep(0.15)

        # Consume events in main thread
        events_consumed = 0
        while not event_queue.empty():
            event = event_queue.get_nowait()
            tracker.apply_event(event)
            events_consumed += 1

        monitor.stop()

        # Verify event was processed
        assert events_consumed >= 1
        assert (
            tracker._molecules["mol_001"].current_stage == ProcessingStage.RDKIT_PARAMS
        )
