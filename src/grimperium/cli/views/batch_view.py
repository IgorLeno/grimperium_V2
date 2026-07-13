"""
Batch processing view for GRIMPERIUM CLI.

Provides UI for:
- Running batch jobs
- Viewing batch status
- Configuring batch parameters

Progress Tracking:
    Uses a 5-stage progress bar (60 chars) with CSV-driven state machine.
    Daemon thread polls CSV every 500ms and communicates via Queue.
"""

import logging
import threading
import time
from dataclasses import asdict
from pathlib import Path
from queue import Queue
from typing import TYPE_CHECKING, Any, ClassVar

from rich.live import Live
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.constants import DATA_DIR
from grimperium.cli.menu import MenuOption
from grimperium.cli.progress_tracker import (
    CSVMonitor,
    ProgressEvent,
    ProgressTracker,
    consume_events,
)
from grimperium.cli.session import MethodExecutionOverrides, RunRef
from grimperium.cli.settings_manager import SettingsManager
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.crest_pm7.batch import (
    BatchCSVManager,
    BatchExecutionManager,
    BatchOutputLayout,
    BatchSortingStrategy,
    BatchStateManager,
    ConformerDetailManager,
    FixedTimeoutProcessor,
)
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.config import PM7Config
from grimperium.runs.models import RunManifest
from grimperium.runs.service import RunService

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

    def _run_service(self) -> RunService:
        service = getattr(self.controller, "__dict__", {}).get("run_service")
        return (
            service
            if isinstance(service, RunService)
            else RunService.from_environment()
        )

    def _method_run_fields(self) -> tuple[str, str, dict[str, Any]]:
        # BatchExecutionManager executa exclusivamente crest_pm7 (PM7-only).
        # Não herdar o método da sessão (ex.: pm7_delta_learning).
        return "crest_pm7", "1.0.0", {"method_id": "crest_pm7"}

    def _execution_overrides_snapshot(self, batch_size: int) -> dict[str, Any]:
        session = getattr(self.controller, "session", None)
        overrides_obj = (
            session.overrides if session is not None else MethodExecutionOverrides()
        )
        overrides = asdict(overrides_obj)
        active = {key: value for key, value in overrides.items() if value is not None}
        active.update(
            {
                "batch_size": batch_size,
                "crest_timeout_minutes": self.crest_timeout,
                "mopac_timeout_minutes": self.mopac_timeout,
            }
        )
        return active

    def _dataset_ref_snapshot(self, csv_path: Path) -> dict[str, str]:
        return {"path": str(csv_path)}

    def _attach_existing_outputs(
        self,
        manifest: RunManifest,
        output_layout: BatchOutputLayout,
    ) -> RunManifest:
        paths = {
            "calculation_results_csv": output_layout.calculation_results_csv,
            "calculation_results_xlsx": output_layout.calculation_results_xlsx,
            "batch_state_csv": output_layout.batch_state_csv,
        }
        existing_paths = {key: path for key, path in paths.items() if path.exists()}
        if not existing_paths:
            return manifest
        return self._run_service().attach_output_paths(manifest.run_id, existing_paths)

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
        """
        if not self.csv_path or not self.csv_path.exists():
            self.show_error("CSV path not configured or file not found")
            return

        if not self.detail_dir:
            self.detail_dir = self.DEFAULT_DETAIL_DIR

        try:
            self._run_batch_with_tracker()
        except Exception as e:
            self.show_error(f"Batch processing failed: {e}")
            logger.exception("Batch processing error")

    def _prepare_batch(self) -> tuple[Any, Any, RunManifest]:
        """Select batch, create RunManifest, then wire outputs under runs/<run_id>/.

        Order (locked): select/create batch → create RunManifest →
        BatchOutputLayout(runs_root / run_id) → BatchExecutionManager with
        output_layout + canonical_run_id.

        Returns:
            Tuple of (BatchExecutionManager, Batch, RunManifest)

        Raises:
            ValueError: If CSV path is not configured
        """
        if self.csv_path is None:
            raise ValueError("CSV path not configured")

        if self.detail_dir is None:
            self.detail_dir = self.DEFAULT_DETAIL_DIR

        # Initialize components
        csv_manager = BatchCSVManager(self.csv_path)
        csv_manager.load_csv()

        session = getattr(self.controller, "session", None)
        overrides = (
            session.overrides if session is not None else MethodExecutionOverrides()
        )
        batch_size = overrides.batch_size or self.batch_size
        crest_timeout = int(overrides.crest_timeout_minutes or self.crest_timeout)
        mopac_timeout = int(overrides.mopac_timeout_minutes or self.mopac_timeout)

        pm7_config = PM7Config()
        if overrides.n_conformers is not None:
            pm7_config.max_conformers = overrides.n_conformers
        state_manager = BatchStateManager(
            self.csv_path.parent / "batch_state.csv",
            pm7_config,
        )
        state_manager.reconcile_molecules(csv_manager.state_seed_rows())
        stuck_mol_ids = state_manager.reset_stuck_running_molecules()
        for mol_id in stuck_mol_ids:
            csv_manager.apply_operational_status(mol_id, MoleculeStatus.PENDING.value)
        if stuck_mol_ids:
            self.console.print(
                f"[yellow]Recovered {len(stuck_mol_ids)} stuck RUNNING molecule(s) -> PENDING[/yellow]"
            )

        detail_manager = ConformerDetailManager(self.detail_dir)
        settings_manager: SettingsManager | None = getattr(
            self.controller, "settings_manager", None
        )
        if settings_manager is not None:
            settings_manager.apply_to_pm7_config(pm7_config)
        xtb_timeout_seconds = (
            int(overrides.xtb_timeout_seconds)
            if overrides.xtb_timeout_seconds is not None
            else (
                int(settings_manager.xtb.timeout_seconds)
                if settings_manager is not None
                and settings_manager.xtb.timeout_seconds is not None
                else None
            )
        )

        processor = FixedTimeoutProcessor(
            config=pm7_config,
            crest_timeout_minutes=crest_timeout,
            mopac_timeout_minutes=mopac_timeout,
            enable_xtb_preopt=pm7_config.xtb_preopt,
            xtb_timeout_seconds=xtb_timeout_seconds,
        )

        # Create batch
        batch_id = csv_manager.generate_batch_id()
        batch = csv_manager.select_batch(
            batch_id=batch_id,
            batch_size=batch_size,
            crest_timeout_minutes=crest_timeout,
            mopac_timeout_minutes=mopac_timeout,
            strategy=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
        )
        state_manager.mark_selected_from_batch(batch)

        method_id, method_version, method_snapshot = self._method_run_fields()
        manifest = self._run_service().create_run(
            property_id="standard_enthalpy_of_formation",
            method_id=method_id,
            method_version=method_version,
            method_snapshot=method_snapshot,
            execution_overrides=self._execution_overrides_snapshot(batch.size),
            dataset_ref=self._dataset_ref_snapshot(self.csv_path),
            model_ref=None,
            molecule_count=batch.size,
        )
        # Outputs autoritativos vivem em runs/<run_id>/ (portáveis).
        # batch_outputs/ legado não é mais a fonte científica.
        output_layout = BatchOutputLayout(self._batch_output_dir(manifest.run_id))
        exec_manager = BatchExecutionManager(
            csv_manager=csv_manager,
            state_manager=state_manager,
            detail_manager=detail_manager,
            pm7_config=pm7_config,
            processor_adapter=processor,
            output_layout=output_layout,
            canonical_run_id=manifest.run_id,
        )

        return exec_manager, batch, manifest

    def _batch_output_dir(self, run_id: str) -> Path:
        """Canonical scientific outputs live under runs/<run_id>/ (portable)."""
        return self._run_service().run_dir(run_id)

    def _run_batch_with_tracker(self) -> None:
        """Run batch with granular 5-stage progress tracking.

        Uses Queue pattern for thread-safe communication between
        CSVMonitor daemon thread and main thread Rich.Live updates.
        """
        # Prepare batch components (select → create_run → layout under runs/)
        exec_manager, batch, manifest = self._prepare_batch()
        output_layout = getattr(exec_manager, "_output_layout", None)
        if not isinstance(output_layout, BatchOutputLayout):
            raise ValueError("Batch output layout not configured")

        csv_path = self.csv_path
        if csv_path is None:
            raise ValueError("CSV path not configured")

        manifest = self._run_service().start_run(manifest.run_id)
        self.controller.session.run = RunRef(
            run_id=manifest.run_id,
            status=manifest.status.value,
        )

        if batch.is_empty:
            cancelled = self._run_service().cancel_run(
                manifest.run_id,
                error="No molecules available for processing",
            )
            self.controller.session.run = RunRef(
                run_id=cancelled.run_id,
                status=cancelled.status.value,
            )
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

        monitor_path = exec_manager.state_manager.state_csv_path

        # Initialize CSV monitor (daemon thread)
        csv_monitor = CSVMonitor(
            csv_path=monitor_path if monitor_path.exists() else csv_path,
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

        # Suppress noisy logs during Rich.Live rendering to avoid display corruption.
        previous_disable = logging.root.manager.disable
        logging.disable(logging.INFO)

        try:
            try:
                with Live(console=self.console, refresh_per_second=3) as live:
                    # Reduced to 3 FPS to avoid display glitches on resize.
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

                    last_size_check = 0.0
                    last_terminal_size: tuple[int, int] | None = None
                    resize_cooldown_until = 0.0
                    size_check_interval = 0.5
                    resize_debounce_seconds = 0.6
                    frame_delay = 1 / 3
                    while not batch_complete.is_set():
                        consume_events(event_queue, tracker)

                        now = time.monotonic()
                        if now - last_size_check >= size_check_interval:
                            size = self.console.size
                            size_tuple = (size.width, size.height)
                            if last_terminal_size is None:
                                last_terminal_size = size_tuple
                            elif size_tuple != last_terminal_size:
                                last_terminal_size = size_tuple
                                resize_cooldown_until = now + resize_debounce_seconds
                            last_size_check = now

                        if now < resize_cooldown_until:
                            time.sleep(frame_delay)
                            continue

                        # Render display
                        display = self._render_batch_display(tracker, frame_idx)
                        live.update(display)

                        frame_idx += 1
                        time.sleep(frame_delay)  # 3 FPS

                    # Final update after batch completes
                    consume_events(event_queue, tracker)
                    if result is not None:
                        tracker.successful = result.success_count
                        tracker.failed = result.rerun_count  # Rerun implies failure
                        tracker.skipped = result.skip_count
                        tracker.total_processed = result.total_count
                    display = self._render_batch_display(tracker, frame_idx)
                    live.update(display)

                    if batch_error is not None:
                        raise batch_error
            except Exception as exc:
                failed = self._run_service().fail_run(manifest.run_id, error=str(exc))
                self.controller.session.run = RunRef(
                    run_id=failed.run_id,
                    status=failed.status.value,
                )
                raise

        finally:
            logging.disable(previous_disable)
            csv_monitor.stop(timeout=1.0)

        # Display final result
        if result is not None:
            self._display_batch_result(result)
            manifest = self._attach_existing_outputs(manifest, output_layout)
            failure_count = min(
                result.failed_count + result.rerun_count + result.skip_count,
                max(result.total_count - result.success_count, 0),
            )
            if result.invalidated:
                finalized = self._run_service().invalidate_run(
                    manifest.run_id,
                    error="ALL_OR_NOTHING reset invalidated scientific completion",
                    success_count=result.success_count,
                    failure_count=failure_count,
                )
            else:
                finalized = self._run_service().complete_run(
                    manifest.run_id,
                    success_count=result.success_count,
                    failure_count=failure_count,
                )
            self.controller.session.run = RunRef(
                run_id=finalized.run_id,
                status=finalized.status.value,
            )

    def _render_batch_display(self, tracker: ProgressTracker, frame_idx: int) -> Panel:
        """Render batch display: header + current molecule + completion history.

        Args:
            tracker: ProgressTracker with current state
            frame_idx: Animation frame index for spinner

        Returns:
            Rich Panel with formatted display
        """
        header = tracker.render_batch_header()
        lines = [header, ""]

        current_mol = tracker.get_current_molecule_id()
        if current_mol:
            lines.append(tracker.render_molecule_line(current_mol, frame_idx))
        else:
            spinner = tracker.get_spinner_frame(frame_idx)
            lines.append(f"{spinner} Aguardando próxima molécula...")
        lines.append("")

        completed = tracker.get_completed_molecules()
        if completed:
            lines.append("Recent Completions:")
            for mol_id, success, elapsed_minutes in completed[-5:]:
                status = "✓" if success else "✗"
                color = "green" if success else "red"
                label = "Completed successfully" if success else "Failed / Rerun"
                colored_status = f"[{color}]{status}[/{color}]"
                mins = int(elapsed_minutes)
                if mins < 1:
                    time_str = "< 1 min"
                elif mins == 1:
                    time_str = "in 1 minute"
                else:
                    time_str = f"in {mins} minutes"
                lines.append(f"  {colored_status} {mol_id} {label}  {time_str}")

        content = "\n".join(lines)

        return Panel(
            content,
            title=f"[bold {COLORS['batch']}]Batch Processing[/bold {COLORS['batch']}]",
            border_style=COLORS["batch"],
        )

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
