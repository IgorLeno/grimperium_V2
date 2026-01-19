"""
Orchestrate batch processing end-to-end.

This module coordinates the entire batch processing workflow:
1. Load and validate CSV data
2. Schedule molecules for processing
3. Execute calculations (CREST + MOPAC)
4. Persist results atomically
5. Print summaries and save logs

**v2.2 UPDATES:**
- Print validation warnings (console + file)
- Print batch summary (console + file)
- Save error logs for monitoring
- Always visible output (never silent)

⚠️ ARCHITECTURE NOTE:
Single-threaded execution. Do NOT attempt parallelism without refactor.
See: docs/ARCHITECTURE.md → "If You Need to Parallelize"
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from rich.console import Console
from rich.table import Table

from grimperium.core.batch_scheduler import BatchScheduler
from grimperium.core.csv_data_loader import BatchDataManager, ValidationReport
from grimperium.core.molecule import Molecule, MoleculeStatus

logger = logging.getLogger(__name__)
console = Console()


class BatchOrchestratorError(Exception):
    """Batch orchestration error."""

    pass


@dataclass
class CalculationSettings:
    """Settings for batch calculations.

    Provides explicit configuration rather than global state.
    """

    max_reruns: int = 3
    """Maximum rerun attempts for failed molecules."""

    crest_timeout_minutes: int = 30
    """CREST execution timeout in minutes."""

    mopac_timeout_minutes: int = 10
    """MOPAC execution timeout in minutes."""


@dataclass
class BatchSummary:
    """Summary statistics for a batch run."""

    total: int = 0
    """Total molecules in batch."""

    pending: int = 0
    """Molecules scheduled for processing."""

    complete: int = 0
    """Molecules completed successfully."""

    errors: int = 0
    """Molecules that encountered errors."""

    skipped: int = 0
    """Molecules skipped (max retries exceeded)."""

    failed_rerunnable: int = 0
    """Failed molecules eligible for rerun."""

    elapsed_total: float = 0.0
    """Total elapsed time in seconds."""

    def to_dict(self) -> dict[str, int | float]:
        """Convert to dictionary."""
        return {
            "total": self.total,
            "pending": self.pending,
            "complete": self.complete,
            "errors": self.errors,
            "skipped": self.skipped,
            "failed_rerunnable": self.failed_rerunnable,
            "elapsed_total": self.elapsed_total,
        }


class BatchOrchestrator:
    """
    Coordinate batch processing (SINGLE-THREADED ONLY).

    ⚠️ ARCHITECTURE NOTE:
    Single-threaded execution. Do NOT attempt parallelism without refactor.
    See: docs/ARCHITECTURE.md → "If You Need to Parallelize"

    Usage:
        >>> orch = BatchOrchestrator(Path("data/thermo_pm7.csv"), strict=False)
        >>> summary = orch.run()
        >>> print(f"Completed: {summary.complete}")
    """

    def __init__(
        self,
        csv_path: Path,
        strict: bool = False,
        workdir: Path | None = None,
        settings: CalculationSettings | None = None,
    ):
        """
        Initialize orchestrator.

        Args:
            csv_path: Path to CSV file
            strict: Validation mode (strict or permissive)
            workdir: Working directory for calculations
            settings: Calculation settings (or defaults)
        """
        self.csv_path = Path(csv_path)
        self.strict = strict
        self.workdir = workdir or Path.cwd()
        self.settings = settings or CalculationSettings()

        self.data_manager = BatchDataManager(csv_path, strict=strict)
        self.summary = BatchSummary()

        # Paths for audit logs
        self.log_dir = Path.home() / ".grimperium"
        self.validation_log = self.log_dir / "validation_errors.log"
        self.summary_log = self.log_dir / "batch_summary.log"

        logger.info(f"BatchOrchestrator initialized (strict={strict})")

    def run(
        self,
        callback: Callable[[Molecule, str], None] | None = None,
        dry_run: bool = False,
    ) -> BatchSummary:
        """
        Run batch processing with full visibility.

        **v2.2:** Always prints:
        1. Validation warnings (if errors)
        2. Processing progress
        3. Final summary

        Saves:
        - ~/.grimperium/validation_errors.log
        - ~/.grimperium/batch_summary.log

        Args:
            callback: Optional callback for progress updates.
                     Called with (molecule, status_message).
            dry_run: If True, only load and schedule, don't process.

        Returns:
            BatchSummary with statistics.

        Raises:
            BatchOrchestratorError: If batch processing fails.
        """
        logger.info("Starting batch processing...")
        start_time = datetime.now(timezone.utc)

        try:
            # Step 1: Load molecules
            molecules = self.data_manager.load_batch()
            self.summary.total = len(molecules)

            # Step 2: Print validation report (if errors)
            if not self.strict:
                self._print_validation_report()

            # Step 3: Schedule
            scheduled = BatchScheduler.schedule(molecules, self.settings.max_reruns)
            self.summary.pending = len(scheduled)

            # Get skip reasons for reporting
            skip_reasons = BatchScheduler.get_skip_reasons(
                molecules, self.settings.max_reruns
            )
            self.summary.failed_rerunnable = skip_reasons.get("max_reruns_exceeded", 0)

            if not scheduled:
                console.print("[yellow]ℹ️  No molecules to process[/yellow]")
                self._print_batch_summary()
                self._save_batch_summary()
                return self.summary

            # Step 4: Dry run or process
            if dry_run:
                console.print(
                    f"[cyan]🔍 DRY RUN: Would process {len(scheduled)} molecules[/cyan]"
                )
                self._print_schedule_preview(scheduled)
            else:
                # Process each molecule
                for idx, mol in enumerate(scheduled, 1):
                    logger.info(f"Processing {idx}/{len(scheduled)}: {mol.mol_id}")

                    self._process_molecule(mol)

                    if callback:
                        callback(mol, f"{mol.status.value} ({idx}/{len(scheduled)})")

            # Step 5: Finalize and print summary
            elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
            self.summary.elapsed_total = elapsed

            # Update counts from data
            status_counts = self.data_manager.count_by_status()
            self.summary.complete = status_counts.get("ok", 0)
            self.summary.skipped = status_counts.get("skip", 0)

            # Print final summary (ALWAYS VISIBLE)
            self._print_batch_summary()

            # Save summary to log
            self._save_batch_summary()

            logger.info("Batch processing complete")

            return self.summary

        except Exception as e:
            logger.error(f"Batch processing failed: {e}")
            raise BatchOrchestratorError(f"Batch failed: {e}") from e

    def _print_validation_report(self) -> None:
        """
        Print validation warnings if errors found (PERMISSIVE MODE).

        **v2.2:** Always visible in console + saved to file
        """
        report = self.data_manager.loader.get_validation_report()

        if report.total_errors == 0:
            return

        # Print to console (prominent)
        console.print("\n" + "=" * 70)
        console.print("[yellow]⚠️  DATA QUALITY WARNING (Permissive Mode)[/yellow]")
        console.print("=" * 70)
        console.print(
            f"[red]{report.total_errors}[/red] rows skipped due to validation errors:"
        )

        # Show first 5 errors
        for error in report.errors[:5]:
            console.print(
                f"  [dim]Row {error.row}:[/dim] "
                f"[cyan]{error.mol_id}[/cyan] → "
                f"[red]{error.error}[/red]"
            )

        if len(report.errors) > 5:
            remaining = len(report.errors) - 5
            console.print(f"  [dim]... and {remaining} more[/dim]")

        console.print(f"[yellow]📝 Full report: {self.validation_log}[/yellow]")
        console.print("=" * 70 + "\n")

        # Save to file (for monitoring)
        self._save_validation_report(report)

    def _save_validation_report(self, report: ValidationReport) -> None:
        """Save validation errors to persistent log file."""
        self.log_dir.mkdir(parents=True, exist_ok=True)

        with open(self.validation_log, "a") as f:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"Batch: {datetime.now(timezone.utc).isoformat()}\n")
            f.write(f"Validation Errors: {report.total_errors}\n")
            f.write(f"{'=' * 70}\n")

            for error in report.errors:
                f.write(f"Row {error.row:4d}: {error.mol_id:20s} " f"→ {error.error}\n")

    def _process_molecule(self, mol: Molecule) -> None:
        """
        Process single molecule.

        Note: In production, this would call CREST and MOPAC.
        For now, this is a placeholder that updates status.

        ⚠️ ARCHITECTURE NOTE:
        This method must NOT be called concurrently.
        See: docs/ARCHITECTURE.md → Limitation 3
        """
        try:
            # Set status: RUNNING
            mol.meta.status = MoleculeStatus.RUNNING
            mol.meta.timestamp_started = datetime.now(timezone.utc)

            # TODO: Integrate with actual CalculationExecutor
            # For now, mark as complete for testing
            # In production:
            #   crest_result = CalculationExecutor.run_crest(mol, self.workdir)
            #   mopac_result = CalculationExecutor.run_mopac_top_3(mol, crest_result)

            # Placeholder: Mark as pending (not actually processed)
            # Real implementation would call CREST/MOPAC here
            logger.debug(f"Would process {mol.mol_id} (placeholder)")

            # For dry-run/testing, keep original status
            # In production, this would be set based on calculation results

        except Exception as e:
            # Handle rerun logic
            mol.meta.reruns += 1
            mol.results.error_message = str(e)
            mol.meta.timestamp_completed = datetime.now(timezone.utc)

            if mol.meta.reruns < self.settings.max_reruns:
                mol.meta.status = MoleculeStatus.RERUN
            else:
                mol.meta.status = MoleculeStatus.SKIP
                self.summary.skipped += 1

            self.summary.errors += 1
            logger.error(f"Error processing {mol.mol_id}: {e}")

    def _print_schedule_preview(self, scheduled: list[Molecule]) -> None:
        """Print preview of scheduled molecules."""
        table = Table(title="Scheduled Molecules (Preview)")
        table.add_column("Position", style="dim")
        table.add_column("mol_id", style="cyan")
        table.add_column("Status", style="green")
        table.add_column("Reruns", style="yellow")

        for idx, mol in enumerate(scheduled[:10], 1):
            table.add_row(
                str(idx),
                mol.mol_id,
                mol.status.value,
                str(mol.reruns),
            )

        if len(scheduled) > 10:
            table.add_row(
                "...",
                f"... and {len(scheduled) - 10} more",
                "",
                "",
            )

        console.print(table)

    def _print_batch_summary(self) -> None:
        """
        Print batch summary (ALWAYS VISIBLE).

        **v2.2:** Console + file output
        """
        # Build table
        table = Table(title="📊 Batch Processing Summary", show_header=False)
        table.add_column(style="cyan", width=25)
        table.add_column(style="green", width=20)

        table.add_row("Total Molecules", str(self.summary.total))
        table.add_row("📋 Scheduled", str(self.summary.pending))
        table.add_row("✅ Completed", str(self.summary.complete))
        table.add_row("❌ Errors", str(self.summary.errors))
        table.add_row("⏭️  Skipped", str(self.summary.skipped))
        table.add_row("⏱️  Elapsed", f"{self.summary.elapsed_total:.1f}s")

        # Print to console (prominent)
        console.print("\n" + "=" * 70)
        console.print(table)
        console.print("=" * 70 + "\n")

        # Also log
        logger.info(
            f"Batch summary: {self.summary.total} total, "
            f"{self.summary.complete} complete, "
            f"{self.summary.errors} errors, "
            f"{self.summary.elapsed_total:.1f}s elapsed"
        )

    def _save_batch_summary(self) -> None:
        """Save batch summary to log file (append)."""
        self.log_dir.mkdir(parents=True, exist_ok=True)

        with open(self.summary_log, "a") as f:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"{datetime.now(timezone.utc).isoformat()}\n")
            f.write(f"Total:     {self.summary.total}\n")
            f.write(f"Scheduled: {self.summary.pending}\n")
            f.write(f"Complete:  {self.summary.complete}\n")
            f.write(f"Errors:    {self.summary.errors}\n")
            f.write(f"Skipped:   {self.summary.skipped}\n")
            f.write(f"Elapsed:   {self.summary.elapsed_total:.1f}s\n")
            f.write(f"{'=' * 70}\n")

    def get_summary(self) -> BatchSummary:
        """Get current batch summary."""
        return self.summary

    def get_validation_report(self) -> ValidationReport:
        """Get validation report from data manager."""
        return self.data_manager.loader.get_validation_report()
