"""Progress tracking for batch molecule processing with thread-safe Queue pattern.

This module provides granular progress tracking for the 5-stage molecular
processing pipeline. It uses a CSV-driven state machine with daemon thread
polling and Queue-based event emission for thread-safe Rich.Live updates.

State Machine Events (CSV-driven):
    1. status: Pending/Selected -> Running    -> "RDKit parameters"      (1/5)
    2. crest_status: NOT_ATTEMPTED -> XTB_PREOPT -> "Pre-optimization"   (2/5)
    3. crest_status: XTB_PREOPT -> CREST_SEARCH -> "CREST conformer search" (3/5)
    4. mopac_status: NOT_ATTEMPTED -> RUNNING -> "MOPAC PM7 calculation" (4/5)
    5. status: Running -> OK                 -> "Final calculations"    (5/5)

Thread Safety:
    - CSVMonitor runs as daemon thread, only writes to Queue
    - ProgressTracker operates in main thread only
    - Queue[ProgressEvent] provides thread-safe communication

Visual Output (30-char progress bar):
    ⠹ mol_00002 ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░░░░░░░░░░░ 3/5 | CREST conformer search
"""

from __future__ import annotations

import logging
import threading
import time
from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from queue import Empty
from typing import TYPE_CHECKING, ClassVar

import pandas as pd
from rich.console import Console

from grimperium.crest_pm7.progress import (
    COMPLETION_STATUSES,
    CREST_STATUS_FAILED,
    CREST_STATUS_NOT_ATTEMPTED,
    CREST_STATUS_PREOPT,
    CREST_STATUS_SEARCH,
    CREST_STATUS_SUCCESS,
    MOPAC_STATUS_FAILED,
    MOPAC_STATUS_NOT_ATTEMPTED,
    MOPAC_STATUS_OK,
    MOPAC_STATUS_RUNNING,
    STATUS_OK,
    STATUS_PENDING,
    STATUS_RERUN,
    STATUS_RUNNING,
    STATUS_SELECTED,
    STATUS_SKIP,
)

if TYPE_CHECKING:
    from queue import Queue as QueueType


logger = logging.getLogger(__name__)

_STATUS_NORMALIZATION = {
    value.lower(): value
    for value in (
        STATUS_PENDING,
        STATUS_SELECTED,
        STATUS_RUNNING,
        STATUS_OK,
        STATUS_RERUN,
        STATUS_SKIP,
    )
}
_CREST_STATUS_NORMALIZATION = {
    value.lower(): value
    for value in (
        CREST_STATUS_NOT_ATTEMPTED,
        CREST_STATUS_PREOPT,
        CREST_STATUS_SEARCH,
        CREST_STATUS_SUCCESS,
        CREST_STATUS_FAILED,
    )
}
_MOPAC_STATUS_NORMALIZATION = {
    value.lower(): value
    for value in (
        MOPAC_STATUS_NOT_ATTEMPTED,
        MOPAC_STATUS_RUNNING,
        MOPAC_STATUS_OK,
        MOPAC_STATUS_FAILED,
    )
}


def _normalize_csv_value(
    value: object,
    *,
    default: str,
    mapping: dict[str, str],
) -> str:
    """Normalize a CSV value to a canonical string.

    Args:
        value: Raw CSV value
        default: Default value when missing
        mapping: Mapping of normalized strings

    Returns:
        Normalized string value
    """
    if value is None or pd.isna(value):
        return default

    text = str(value).strip()
    if not text:
        return default

    return mapping.get(text.lower(), text)


class ProcessingStage(IntEnum):
    """Processing stages mapped to CSV column transitions.

    Each stage represents a step in the molecular processing pipeline.
    Values are used for progress calculation (stage * 6 = filled chars).
    """

    NOT_STARTED = 0
    RDKIT_PARAMS = 1  # status: Pending/Selected -> Running
    XTB_PREOPT = 2  # crest_status: NOT_ATTEMPTED -> XTB_PREOPT
    CREST_SEARCH = 3  # crest_status: XTB_PREOPT -> CREST_SEARCH
    MOPAC_CALC = 4  # mopac_status: NOT_ATTEMPTED -> RUNNING
    FINAL_CALC = 5  # status: Running -> OK


@dataclass(frozen=True)
class StageTransition:
    """Defines a CSV column transition that triggers a stage change.

    Attributes:
        column: CSV column to monitor (status, crest_status, mopac_status)
        from_values: Previous values that trigger transition
        to_value: New value that triggers transition
        stage: ProcessingStage to transition to
        label: Human-readable label for display
    """

    column: str
    from_values: tuple[str, ...]
    to_value: str
    stage: ProcessingStage
    label: str


# Immutable tuple of 5 stage transitions (state machine definition)
EVENTS: tuple[StageTransition, ...] = (
    StageTransition(
        column="status",
        from_values=(STATUS_PENDING, STATUS_SELECTED),
        to_value=STATUS_RUNNING,
        stage=ProcessingStage.RDKIT_PARAMS,
        label="RDKit parameters",
    ),
    StageTransition(
        column="crest_status",
        from_values=(CREST_STATUS_NOT_ATTEMPTED,),
        to_value=CREST_STATUS_PREOPT,
        stage=ProcessingStage.XTB_PREOPT,
        label="Pre-optimization",
    ),
    StageTransition(
        column="crest_status",
        from_values=(CREST_STATUS_PREOPT, CREST_STATUS_NOT_ATTEMPTED),
        to_value=CREST_STATUS_SEARCH,
        stage=ProcessingStage.CREST_SEARCH,
        label="CREST conformer search",
    ),
    StageTransition(
        column="mopac_status",
        from_values=(MOPAC_STATUS_NOT_ATTEMPTED,),
        to_value=MOPAC_STATUS_RUNNING,
        stage=ProcessingStage.MOPAC_CALC,
        label="MOPAC PM7 calculation",
    ),
    StageTransition(
        column="status",
        from_values=(STATUS_RUNNING,),
        to_value=STATUS_OK,
        stage=ProcessingStage.FINAL_CALC,
        label="Final calculations",
    ),
)


@dataclass
class ProgressEvent:
    """Event emitted when molecule changes processing stage.

    Attributes:
        mol_id: Molecule identifier
        new_stage: New ProcessingStage after transition
        status: Optional status update (OK/Rerun/Skip)
        timestamp: Unix timestamp when event was created
    """

    mol_id: str
    new_stage: ProcessingStage | None = None
    status: str | None = None
    timestamp: float = field(default_factory=time.time)


@dataclass
class MoleculeProgress:
    """Track progress for a single molecule through processing stages.

    Uses last_csv_state to detect transitions by comparing previous
    and current CSV column values.

    Attributes:
        mol_id: Molecule identifier
        current_stage: Current ProcessingStage
        last_csv_state: Cache of previous CSV column values
        completed: Whether molecule has finished processing
        error: Error message if processing failed
    """

    mol_id: str
    current_stage: ProcessingStage = ProcessingStage.NOT_STARTED
    last_csv_state: dict[str, str] = field(
        default_factory=lambda: {
            "status": STATUS_PENDING,
            "crest_status": CREST_STATUS_NOT_ATTEMPTED,
            "mopac_status": MOPAC_STATUS_NOT_ATTEMPTED,
        }
    )
    completed: bool = False
    error: str | None = None

    def _detect_stage(self, current_row: dict[str, str]) -> ProcessingStage | None:
        """Compare last_csv_state with current_row to detect stage transition.

        Iterates through EVENTS and checks if any transition matches:
        - Previous value matches event.from_values
        - Current value matches event.to_value

        Args:
            current_row: Dict with current status, crest_status, mopac_status

        Returns:
            New ProcessingStage if transition detected, None otherwise
        """
        for event in EVENTS:
            prev_val = self.last_csv_state.get(event.column, "")
            curr_val = current_row.get(event.column, "")

            if prev_val in event.from_values and curr_val == event.to_value:
                return event.stage

        return None

    def _infer_stage(self, current_row: dict[str, str]) -> ProcessingStage | None:
        """Infer stage from current CSV values when transitions are missed.

        Args:
            current_row: Dict with current status, crest_status, mopac_status

        Returns:
            Inferred ProcessingStage or None if no inference applies
        """
        status = current_row.get("status", "")
        crest_status = current_row.get("crest_status", "")
        mopac_status = current_row.get("mopac_status", "")

        if status == STATUS_OK:
            return ProcessingStage.FINAL_CALC
        if mopac_status in {
            MOPAC_STATUS_RUNNING,
            MOPAC_STATUS_OK,
            MOPAC_STATUS_FAILED,
        }:
            return ProcessingStage.MOPAC_CALC
        if crest_status in {
            CREST_STATUS_SEARCH,
            CREST_STATUS_SUCCESS,
            CREST_STATUS_FAILED,
        }:
            return ProcessingStage.CREST_SEARCH
        if crest_status == CREST_STATUS_PREOPT:
            return ProcessingStage.XTB_PREOPT
        if status == STATUS_RUNNING:
            return ProcessingStage.RDKIT_PARAMS

        return None

    def update_from_csv_row(self, row: pd.Series) -> ProcessingStage | None:
        """Check CSV row for state transitions and update internal state.

        Args:
            row: Pandas Series with CSV column values

        Returns:
            New ProcessingStage if transition detected, None otherwise
        """
        current_row = {
            "status": _normalize_csv_value(
                row.get("status"),
                default=STATUS_PENDING,
                mapping=_STATUS_NORMALIZATION,
            ),
            "crest_status": _normalize_csv_value(
                row.get("crest_status"),
                default=CREST_STATUS_NOT_ATTEMPTED,
                mapping=_CREST_STATUS_NORMALIZATION,
            ),
            "mopac_status": _normalize_csv_value(
                row.get("mopac_status"),
                default=MOPAC_STATUS_NOT_ATTEMPTED,
                mapping=_MOPAC_STATUS_NORMALIZATION,
            ),
        }

        new_stage = self._detect_stage(current_row)

        if new_stage is None:
            inferred_stage = self._infer_stage(current_row)
            if inferred_stage is not None and inferred_stage > self.current_stage:
                new_stage = inferred_stage

        if new_stage is not None:
            self.current_stage = new_stage

        # Update last_csv_state for next comparison
        self.last_csv_state = current_row

        # Check completion (terminal states)
        status = current_row["status"]
        if status in COMPLETION_STATUSES:
            self.completed = True

        return new_stage


class ProgressTracker:
    """Main progress tracker for batch molecule processing.

    Operates in main thread only. Receives ProgressEvents from Queue
    and updates molecule states. Renders 30-character progress bars.

    Thread Safety:
        NOT thread-safe internally - all calls must be from main thread.
        Uses Queue to receive events from CSVMonitor daemon thread.

    Attributes:
        PROGRESS_BAR_WIDTH: Total width of progress bar (30 chars)
        STAGE_WIDTH: Width per stage (6 chars)
        FILLED_CHAR: Character for completed stages
        EMPTY_CHAR: Character for pending stages
    """

    PROGRESS_BAR_WIDTH: ClassVar[int] = 30
    STAGE_WIDTH: ClassVar[int] = 6  # 30 / 5 stages = 6 chars each
    FILLED_CHAR: ClassVar[str] = "▓"
    EMPTY_CHAR: ClassVar[str] = "░"

    SPINNER_FRAMES: ClassVar[tuple[str, ...]] = (
        "⠋",
        "⠙",
        "⠹",
        "⠸",
        "⠼",
        "⠴",
        "⠦",
        "⠧",
        "⠇",
        "⠏",
    )

    def __init__(
        self,
        console: Console,
        batch_size: int,
    ) -> None:
        """Initialize progress tracker.

        Args:
            console: Rich Console for output
            batch_size: Number of molecules in batch
        """
        self.console = console
        self.batch_size = batch_size
        self._molecules: dict[str, MoleculeProgress] = {}

        # Batch statistics
        self.total_processed: int = 0
        self.successful: int = 0
        self.failed: int = 0
        self.skipped: int = 0

    def register_molecule(self, mol_id: str) -> None:
        """Register a molecule for tracking.

        Args:
            mol_id: Molecule identifier to track
        """
        self._molecules[mol_id] = MoleculeProgress(mol_id=mol_id)

    def apply_event(self, event: ProgressEvent) -> None:
        """Apply a progress event to update molecule state.

        Args:
            event: ProgressEvent with mol_id and new_stage
        """
        if event.mol_id not in self._molecules:
            return

        progress = self._molecules[event.mol_id]
        if event.new_stage is not None:
            progress.current_stage = event.new_stage

        if event.status is not None:
            if event.status == STATUS_OK:
                self.mark_completed(event.mol_id, success=True)
            elif event.status == STATUS_RERUN:
                self.mark_completed(event.mol_id, success=False)
            elif event.status == STATUS_SKIP:
                self.mark_skipped(event.mol_id)

    def mark_completed(self, mol_id: str, *, success: bool) -> None:
        """Mark molecule as completed and update statistics.

        Removes molecule from active tracking after updating stats.

        Args:
            mol_id: Molecule identifier
            success: Whether processing succeeded
        """
        if mol_id in self._molecules:
            self._molecules[mol_id].completed = True
            del self._molecules[mol_id]  # Auto-cleanup

        self.total_processed += 1
        if success:
            self.successful += 1
        else:
            self.failed += 1

    def mark_skipped(self, mol_id: str) -> None:
        """Mark molecule as skipped.

        Args:
            mol_id: Molecule identifier
        """
        if mol_id in self._molecules:
            del self._molecules[mol_id]

        self.total_processed += 1
        self.skipped += 1

    def get_active_molecule_ids(self) -> list[str]:
        """Get list of currently active (not completed) molecule IDs.

        Returns:
            List of molecule IDs being tracked
        """
        return list(self._molecules.keys())

    def render_progress_bar(self, mol_id: str) -> str:
        """Render 30-character progress bar for a molecule.

        Each stage fills 6 characters. Stage 0 = 0 filled,
        Stage 5 = 30 filled.

        Args:
            mol_id: Molecule identifier

        Returns:
            30-character string with FILLED_CHAR and EMPTY_CHAR
        """
        progress = self._molecules.get(mol_id)
        if progress is None:
            return self.EMPTY_CHAR * self.PROGRESS_BAR_WIDTH

        filled = progress.current_stage.value * self.STAGE_WIDTH
        empty = self.PROGRESS_BAR_WIDTH - filled

        return self.FILLED_CHAR * filled + self.EMPTY_CHAR * empty

    def get_stage_label(self, mol_id: str) -> str:
        """Get current stage label for molecule.

        Args:
            mol_id: Molecule identifier

        Returns:
            Human-readable stage label
        """
        progress = self._molecules.get(mol_id)
        if progress is None or progress.current_stage == ProcessingStage.NOT_STARTED:
            return "Pending"

        event = EVENTS[progress.current_stage.value - 1]
        return event.label

    def get_spinner_frame(self, frame_idx: int) -> str:
        """Get spinner animation frame.

        Args:
            frame_idx: Frame index (will be wrapped)

        Returns:
            Single character spinner frame
        """
        return self.SPINNER_FRAMES[frame_idx % len(self.SPINNER_FRAMES)]

    def render_batch_header(self) -> str:
        """Render 7-line batch status header with bullets.

        Format:
            Batch Status:
              - Batch size: 6
              - Total processed: 2
              - Successful: 2
              - Failed: 0
              - Skipped: 0
              - Success rate: 100.0%
              - Total progress: 33%

        Returns:
            Formatted header string
        """
        success_rate = (
            (self.successful / self.total_processed * 100)
            if self.total_processed > 0
            else 0.0
        )
        total_progress = (
            int(self.total_processed / self.batch_size * 100)
            if self.batch_size > 0
            else 0
        )

        return (
            "Batch Status:\n"
            f"  - Batch size: {self.batch_size}\n"
            f"  - Total processed: {self.total_processed}\n"
            f"  - Successful: {self.successful}\n"
            f"  - Failed: {self.failed}\n"
            f"  - Skipped: {self.skipped}\n"
            f"  - Success rate: {success_rate:.1f}%\n"
            f"  - Total progress: {total_progress}%"
        )

    def render_molecule_line(self, mol_id: str, frame_idx: int) -> str:
        """Render single molecule progress line.

        Format:
            ⠹ mol_00002 ▓▓▓▓▓▓▓▓▓▓▓▓░░░░░░░░░░░░░░░░░░ 2/5 | Pre-optimization

        Args:
            mol_id: Molecule identifier
            frame_idx: Animation frame index

        Returns:
            Formatted progress line
        """
        progress = self._molecules.get(mol_id)
        if progress is None:
            return f"  ? {mol_id} (unknown)"

        spinner = self.get_spinner_frame(frame_idx)
        bar = self.render_progress_bar(mol_id)
        stage_num = progress.current_stage.value
        label = self.get_stage_label(mol_id)

        return f"  {spinner} {mol_id} {bar} {stage_num}/5 | {label}"


class CSVMonitor:
    """Daemon thread that polls CSV for state changes and enqueues events.

    Thread Safety:
        Only writes to Queue (thread-safe). Never touches ProgressTracker.
        Uses mtime optimization to skip unchanged files.

    Attributes:
        DEFAULT_POLL_INTERVAL_MS: Default polling interval (500ms)
    """

    DEFAULT_POLL_INTERVAL_MS: ClassVar[int] = 500

    def __init__(
        self,
        csv_path: Path,
        event_queue: QueueType[ProgressEvent],
        poll_interval_ms: int = DEFAULT_POLL_INTERVAL_MS,
    ) -> None:
        """Initialize CSV monitor.

        Args:
            csv_path: Path to batch tracking CSV
            event_queue: Queue to emit ProgressEvents
            poll_interval_ms: Polling interval in milliseconds
        """
        self.csv_path = Path(csv_path)
        self.event_queue = event_queue
        self.poll_interval_ms = poll_interval_ms

        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None
        self._last_mtime: float = 0

        # Cache previous state per molecule
        self._molecule_states: dict[str, MoleculeProgress] = {}

    def register_molecule(self, mol_id: str) -> None:
        """Register molecule for monitoring.

        Args:
            mol_id: Molecule identifier to monitor
        """
        self._molecule_states[mol_id] = MoleculeProgress(mol_id=mol_id)

    def start(self) -> None:
        """Start the monitoring thread (daemon=True)."""
        if self._thread is not None and self._thread.is_alive():
            return

        self._stop_event.clear()
        self._thread = threading.Thread(
            target=self._poll_loop,
            name="CSVMonitor",
            daemon=True,
        )
        self._thread.start()

    def stop(self, timeout: float = 1.0) -> None:
        """Stop the monitoring thread gracefully.

        Args:
            timeout: Maximum seconds to wait for thread to stop
        """
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=timeout)
            self._thread = None

    def _poll_loop(self) -> None:
        """Main polling loop (runs in daemon thread)."""
        while not self._stop_event.is_set():
            try:
                self._poll_csv()
            except Exception as e:
                logger.warning(f"CSV poll error: {e}")

            # Sleep with interruptible wait
            self._stop_event.wait(self.poll_interval_ms / 1000.0)

    def _poll_csv(self) -> None:
        """Read CSV and enqueue change events.

        Optimization: Only re-reads if file mtime changed.
        Safety: Catches all file/parse errors gracefully.
        """
        if not self.csv_path.exists():
            return

        # Check modification time (optimization)
        try:
            mtime = self.csv_path.stat().st_mtime
            if mtime == self._last_mtime:
                return
            self._last_mtime = mtime
        except OSError:
            return

        # Read CSV with explicit dtypes
        try:
            df = pd.read_csv(
                self.csv_path,
                dtype={
                    "mol_id": str,
                    "status": str,
                    "crest_status": str,
                    "mopac_status": str,
                },
            )
        except pd.errors.EmptyDataError:
            logger.debug("CSV empty, skipping poll")
            return
        except pd.errors.ParserError as e:
            logger.warning(f"CSV parse error: {e}")
            return
        except FileNotFoundError:
            return
        except PermissionError:
            logger.debug("CSV locked, skipping poll")
            return
        except UnicodeDecodeError:
            logger.warning("CSV encoding error")
            return

        # Detect changes and enqueue events
        try:
            for mol_id, progress in self._molecule_states.items():
                if progress.completed:
                    continue

                mask = df["mol_id"] == mol_id
                if not mask.any():
                    continue

                row = df[mask].iloc[0]
                prev_status = progress.last_csv_state.get("status", "")
                new_stage = progress.update_from_csv_row(row)
                current_status = progress.last_csv_state.get("status", "")

                status_changed = prev_status != current_status
                completion_status = (
                    current_status if current_status in COMPLETION_STATUSES else None
                )

                if new_stage is not None or (status_changed and completion_status):
                    self.event_queue.put(
                        ProgressEvent(
                            mol_id=mol_id,
                            new_stage=new_stage,
                            status=completion_status if status_changed else None,
                        )
                    )
                    if new_stage is not None:
                        logger.debug(f"[{mol_id}] Stage changed to {new_stage.name}")
        except KeyError as e:
            logger.error(
                f"Missing required column in CSV: {e}. "
                f"Available columns: {list(df.columns)}"
            )
            return


def consume_events(
    event_queue: QueueType[ProgressEvent],
    tracker: ProgressTracker,
) -> int:
    """Consume all pending events from queue and apply to tracker.

    Non-blocking: Returns immediately if queue is empty.

    Args:
        event_queue: Queue with ProgressEvents (stage and/or completion status)
        tracker: ProgressTracker to apply events to

    Returns:
        Number of events consumed
    """
    events_consumed = 0
    while True:
        try:
            event = event_queue.get_nowait()
            tracker.apply_event(event)
            events_consumed += 1
        except Empty:
            break
    return events_consumed
