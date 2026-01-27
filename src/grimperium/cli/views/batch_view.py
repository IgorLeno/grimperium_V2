"""
Batch processing view for GRIMPERIUM CLI.

Provides UI for:
- Running batch jobs
- Viewing batch status
- Configuring batch parameters

Progress Tracking:
    Uses a 5-stage progress bar (30 chars) with CSV-driven state machine.
    Daemon thread polls CSV every 500ms and communicates via Queue.
"""

import logging
import threading
import time
from pathlib import Path
from queue import Queue
from typing import TYPE_CHECKING, Any, ClassVar

from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from grimperium.cli.constants import DATA_DIR
from grimperium.cli.menu import MenuOption
from grimperium.cli.progress_tracker import (
    CSVMonitor,
    ProgressEvent,
    ProgressTracker,
    consume_events,
)
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController

logger = logging.getLogger(__name__)


class BatchView(BaseView):
    """View for batch processing operations."""

    name: ClassVar[str] = "batch"
    title: ClassVar[str] = "Batch Processing"
    icon: ClassVar[str] = "📦"
    color: ClassVar[str] = COLORS["batch"]

    # Default paths
    DEFAULT_CSV_PATH = DATA_DIR / "batch_tracking.csv"
    DEFAULT_DETAIL_DIR = DATA_DIR / "conformer_details"

    def __init__(self, controller: "CliController") -> None:
        """Initialize batch view.

        Args:
            controller: CLI controller
        """
        super().__init__(controller)
        self.csv_path: Path | None = None
        self.detail_dir: Path | None = None
        self.batch_size: int = 10
        self.crest_timeout: int = 30
        self.mopac_timeout: int = 60

    def render(self) -> None:
        """Render batch processing view."""
        self.show_header()

        # Show current configuration
        self._render_config()

        # Show status if CSV exists
        if self.csv_path and self.csv_path.exists():
            self._render_status()

    def _render_config(self) -> None:
        """Render current batch configuration."""
        table = Table(
            title="Batch Configuration",
            show_header=False,
            border_style=COLORS["muted"],
        )
        table.add_column("Setting", style=COLORS["batch"])
        table.add_column("Value", style=COLORS["highlight"])

        csv_status = (
            "✓ Found" if self.csv_path and self.csv_path.exists() else "✗ Not set"
        )
        detail_status = (
            "✓ Exists" if self.detail_dir and self.detail_dir.exists() else "✗ Not set"
        )

        table.add_row(
            "CSV Path", str(self.csv_path) if self.csv_path else "Not configured"
        )
        table.add_row("CSV Status", csv_status)
        table.add_row(
            "Detail Dir", str(self.detail_dir) if self.detail_dir else "Not configured"
        )
        table.add_row("Detail Status", detail_status)
        table.add_row("Batch Size", str(self.batch_size))
        table.add_row("CREST Timeout", f"{self.crest_timeout} min")
        table.add_row("MOPAC Timeout", f"{self.mopac_timeout} min")

        self.console.print(table)
        self.console.print()

    def _render_status(self) -> None:
        """Render current batch status from CSV."""
        try:
            from grimperium.crest_pm7.batch import BatchCSVManager

            if self.csv_path is None:
                self.console.print("[yellow]No CSV path configured[/]")
                return
            manager = BatchCSVManager(self.csv_path)
            manager.load_csv()
            counts = manager.get_status_counts()

            table = Table(
                title="Dataset Status",
                border_style=COLORS["muted"],
            )
            table.add_column("Status", style="bold")
            table.add_column("Count", justify="right")
            table.add_column("Percentage", justify="right")

            total = sum(counts.values())
            for status, count in sorted(counts.items()):
                pct = 100 * count / total if total > 0 else 0
                color = self._get_status_color(status)
                table.add_row(
                    f"[{color}]{status}[/{color}]",
                    str(count),
                    f"{pct:.1f}%",
                )

            table.add_row("─" * 10, "─" * 5, "─" * 8)
            table.add_row("[bold]Total[/bold]", str(total), "100%")

            self.console.print(table)
            self.console.print()

        except Exception as e:
            self.show_error(f"Failed to load status: {e}")

    def _get_status_color(self, status: str) -> str:
        """Get color for status display."""
        colors = {
            "Pending": COLORS["muted"],
            "Selected": COLORS["about"],
            "Running": COLORS["warning"],
            "OK": COLORS["success"],
            "Rerun": COLORS["warning"],
            "Skip": COLORS["error"],
        }
        return colors.get(status, COLORS["highlight"])

    def get_menu_options(self) -> list[MenuOption]:
        """Get menu options for batch view."""
        options = [
            MenuOption(
                label="Configure Paths",
                value="configure",
                icon="⚙️",
                description="Set CSV and detail directory paths",
            ),
            MenuOption(
                label="Set Batch Parameters",
                value="params",
                icon="📊",
                description="Configure batch size and timeouts",
            ),
        ]

        if self.csv_path and self.csv_path.exists():
            options.extend(
                [
                    MenuOption(
                        label="Run Batch",
                        value="run",
                        icon="▶️",
                        description="Execute next batch of molecules",
                    ),
                    MenuOption(
                        label="View Status",
                        value="status",
                        icon="📈",
                        description="View detailed processing status",
                    ),
                ]
            )

        options.append(
            MenuOption(
                label="Back to Main Menu",
                value="back",
                icon=ICONS["back"],
                description="Return to main menu",
            )
        )

        return options

    def handle_action(self, action: str) -> str | None:
        """Handle menu action."""
        if action == "back":
            return "main"
        elif action == "configure":
            self._configure_paths()
        elif action == "params":
            self._configure_params()
        elif action == "run":
            self._run_batch()
        elif action == "status":
            self._show_detailed_status()

        return None

    def _configure_paths(self) -> None:
        """Configure CSV and detail directory paths."""
        self.console.print()
        self.console.print("[bold]Configure Batch Paths[/bold]")
        self.console.print()

        # CSV Path
        default_csv = str(self.DEFAULT_CSV_PATH)
        csv_input = self.console.input(
            f"CSV Path [{COLORS['muted']}][{default_csv}][/{COLORS['muted']}]: "
        ).strip()
        self.csv_path = Path(csv_input) if csv_input else self.DEFAULT_CSV_PATH

        # Detail Directory
        default_detail = str(self.DEFAULT_DETAIL_DIR)
        detail_input = self.console.input(
            f"Detail Dir [{COLORS['muted']}][{default_detail}][/{COLORS['muted']}]: "
        ).strip()
        self.detail_dir = (
            Path(detail_input) if detail_input else self.DEFAULT_DETAIL_DIR
        )

        self.show_success("Paths configured successfully")

    def _configure_params(self) -> None:
        """Configure batch parameters."""
        self.console.print()
        self.console.print("[bold]Configure Batch Parameters[/bold]")
        self.console.print()

        try:
            # Batch size
            size_input = self.console.input(
                f"Batch Size [{COLORS['muted']}][{self.batch_size}][/{COLORS['muted']}]: "
            ).strip()
            if size_input:
                self.batch_size = max(1, int(size_input))

            # CREST timeout
            crest_input = self.console.input(
                f"CREST Timeout (min) [{COLORS['muted']}][{self.crest_timeout}][/{COLORS['muted']}]: "
            ).strip()
            if crest_input:
                self.crest_timeout = max(1, int(crest_input))

            # MOPAC timeout
            mopac_input = self.console.input(
                f"MOPAC Timeout (min) [{COLORS['muted']}][{self.mopac_timeout}][/{COLORS['muted']}]: "
            ).strip()
            if mopac_input:
                self.mopac_timeout = max(1, int(mopac_input))

            self.show_success("Parameters configured successfully")

        except ValueError as e:
            self.show_error(f"Invalid input: {e}")

    def _run_batch(self) -> None:
        """Run a batch of molecules with granular progress tracking.

        Uses a 5-stage progress bar with CSV-driven state machine.
        Falls back to legacy progress bar if new system fails.
        """
        if not self.csv_path or not self.csv_path.exists():
            self.show_error("CSV path not configured or file not found")
            return

        if not self.detail_dir:
            self.detail_dir = self.DEFAULT_DETAIL_DIR

        try:
            self._run_batch_with_tracker()
        except Exception as e:
            self.console.print(
                "[bold yellow]Warning: Progress tracker failed, "
                "using legacy mode[/bold yellow]"
            )
            logger.exception(f"Progress tracker error: {e}")
            self._run_batch_legacy()

    def _prepare_batch(self) -> tuple[Any, Any]:
        """Prepare batch components and create batch.

        Returns:
            Tuple of (BatchExecutionManager, Batch)

        Raises:
            ValueError: If CSV path is not configured
        """
        from grimperium.crest_pm7.batch import (
            BatchCSVManager,
            BatchExecutionManager,
            BatchSortingStrategy,
            ConformerDetailManager,
            FixedTimeoutProcessor,
        )
        from grimperium.crest_pm7.config import PM7Config

        if self.csv_path is None:
            raise ValueError("CSV path not configured")

        if self.detail_dir is None:
            self.detail_dir = self.DEFAULT_DETAIL_DIR

        # Initialize components
        csv_manager = BatchCSVManager(self.csv_path)
        csv_manager.load_csv()

        detail_manager = ConformerDetailManager(self.detail_dir)
        pm7_config = PM7Config()

        processor = FixedTimeoutProcessor(
            config=pm7_config,
            crest_timeout_minutes=self.crest_timeout,
            mopac_timeout_minutes=self.mopac_timeout,
        )

        exec_manager = BatchExecutionManager(
            csv_manager=csv_manager,
            detail_manager=detail_manager,
            pm7_config=pm7_config,
            processor_adapter=processor,
        )

        # Create batch
        batch_id = csv_manager.generate_batch_id()
        batch = csv_manager.select_batch(
            batch_id=batch_id,
            batch_size=self.batch_size,
            crest_timeout_minutes=self.crest_timeout,
            mopac_timeout_minutes=self.mopac_timeout,
            strategy=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
        )

        return exec_manager, batch

    def _run_batch_with_tracker(self) -> None:
        """Run batch with granular 5-stage progress tracking.

        Uses Queue pattern for thread-safe communication between
        CSVMonitor daemon thread and main thread Rich.Live updates.
        """
        # Prepare batch components
        exec_manager, batch = self._prepare_batch()

        if batch.is_empty:
            self.show_success("No molecules available for processing")
            return

        self.console.print()
        self.console.print(
            f"[bold]Starting batch {batch.batch_id}[/bold]: {batch.size} molecules"
        )
        self.console.print()

        # Create event queue for thread-safe communication
        event_queue: Queue[ProgressEvent] = Queue()

        # Initialize progress tracker (main thread only)
        tracker = ProgressTracker(
            console=self.console,
            batch_size=batch.size,
        )

        # Initialize CSV monitor (daemon thread)
        csv_monitor = CSVMonitor(
            csv_path=self.csv_path,
            event_queue=event_queue,
            poll_interval_ms=500,
        )

        # Register all molecules
        for mol in batch.molecules:
            tracker.register_molecule(mol.mol_id)
            csv_monitor.register_molecule(mol.mol_id)

        frame_idx = 0
        result = None

        csv_monitor.start()

        try:
            with Live(console=self.console, refresh_per_second=10) as live:
                # Start batch execution in a way that allows progress updates
                # We need to run execute_batch and update display concurrently

                # Use threading to run batch in background while updating display
                batch_complete = threading.Event()
                batch_error: Exception | None = None

                def run_batch() -> None:
                    nonlocal result, batch_error
                    try:
                        result = exec_manager.execute_batch(batch)
                    except Exception as e:
                        batch_error = e
                    finally:
                        batch_complete.set()

                batch_thread = threading.Thread(target=run_batch, daemon=True)
                batch_thread.start()

                # Main loop: consume events and update display
                while not batch_complete.is_set():
                    # Consume all pending events (non-blocking)
                    consume_events(event_queue, tracker)

                    # Render display
                    display = self._render_batch_display(tracker, frame_idx)
                    live.update(display)

                    frame_idx += 1
                    time.sleep(0.1)  # 10 FPS

                # Final update after batch completes
                consume_events(event_queue, tracker)
                display = self._render_batch_display(tracker, frame_idx)
                live.update(display)

                if batch_error is not None:
                    raise batch_error

        finally:
            csv_monitor.stop(timeout=1.0)

        # Update tracker stats from batch result
        if result is not None:
            # Sync tracker stats with actual batch results
            tracker.successful = result.success_count
            tracker.failed = result.rerun_count  # Rerun implies failure
            tracker.skipped = result.skip_count
            tracker.total_processed = result.total_count

        # Display final result
        if result is not None:
            self._display_batch_result(result)

    def _render_batch_display(self, tracker: ProgressTracker, frame_idx: int) -> Panel:
        """Render complete batch display with header and progress bars.

        Args:
            tracker: ProgressTracker with current state
            frame_idx: Animation frame index for spinner

        Returns:
            Rich Panel with formatted display
        """
        # Build header (7 lines)
        header = tracker.render_batch_header()

        # Build molecule progress lines
        lines = [header, ""]  # Header + blank line

        for mol_id in tracker.get_active_molecule_ids():
            line = tracker.render_molecule_line(mol_id, frame_idx)
            lines.append(line)

        content = "\n".join(lines)

        return Panel(
            content,
            title=f"[bold {COLORS['batch']}]Batch Processing[/bold {COLORS['batch']}]",
            border_style=COLORS["batch"],
        )

    def _run_batch_legacy(self) -> None:
        """Run batch with legacy simple progress bar.

        Fallback when new progress tracker fails.
        """
        # Prepare batch components
        exec_manager, batch = self._prepare_batch()

        if batch.is_empty:
            self.show_success("No molecules available for processing")
            return

        self.console.print()
        self.console.print(
            f"[bold]Starting batch {batch.batch_id}[/bold]: {batch.size} molecules"
        )
        self.console.print()

        # Run with legacy progress bar
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            console=self.console,
        ) as progress:
            task = progress.add_task("Processing", total=batch.size)

            def update_progress(mol_id: str, current: int, total: int) -> None:
                progress.update(
                    task, completed=current, description=f"Processing {mol_id}"
                )

            result = exec_manager.execute_batch(
                batch, progress_callback=update_progress
            )

        # Display result
        self._display_batch_result(result)

    def _display_batch_result(self, result: Any) -> None:
        """Display batch execution result."""
        self.console.print()

        panel_content = (
            f"[bold]Batch:[/bold] {result.batch_id}\n"
            f"[bold]Total:[/bold] {result.total_count}\n"
            f"[{COLORS['success']}]OK:[/{COLORS['success']}] {result.success_count}\n"
            f"[{COLORS['warning']}]Rerun:[/{COLORS['warning']}] {result.rerun_count}\n"
            f"[{COLORS['error']}]Skip:[/{COLORS['error']}] {result.skip_count}\n"
            f"[bold]Time:[/bold] {result.total_time:.1f}s\n"
            f"[bold]Success Rate:[/bold] {result.success_rate:.1f}%"
        )

        if result.min_hof is not None:
            panel_content += (
                f"\n\n[bold]Energy Range:[/bold]\n"
                f"  Min: {result.min_hof:.2f} kcal/mol ({result.min_hof_mol_id})\n"
                f"  Max: {result.max_hof:.2f} kcal/mol ({result.max_hof_mol_id})"
            )

        color = COLORS["success"] if result.success_rate >= 80 else COLORS["warning"]
        self.console.print(
            Panel(
                panel_content,
                title=f"[bold {color}]Batch Complete[/bold {color}]",
                border_style=color,
            )
        )

        self.wait_for_enter()

    def _show_detailed_status(self) -> None:
        """Show detailed processing status."""
        self._render_status()
        self.wait_for_enter()
