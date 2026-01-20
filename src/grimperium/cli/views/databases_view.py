"""
Databases view for GRIMPERIUM CLI.

Displays and manages molecular databases.
"""

import json
import os
from datetime import date, datetime
from typing import TYPE_CHECKING, Any

import pandas as pd
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

from grimperium.cli.constants import DATA_DIR, PHASE_A_RESULTS_FILE
from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.mock_data import DATABASES, Database
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController


class DatabasesView(BaseView):
    """View for managing molecular databases."""

    name = "databases"
    title = "Databases"
    icon = ICONS["databases"]
    color = COLORS["databases"]

    def __init__(self, controller: "CliController") -> None:
        """Initialize the databases view."""
        super().__init__(controller)
        self.selected_db: Database | None = None

    @staticmethod
    def load_real_phase_a_results() -> Database | None:
        """Load real Phase A CREST PM7 results with CSV fallback.

        Priority:
        1. Load from phase_a_results.json if exists and has n_molecules > 0
        2. Fallback to counting rows in data/thermo_pm7.csv
        3. Return None if no data source available

        Returns:
            Database object with real data, or None if no data available.
        """
        molecules_count = 0
        last_updated: date | None = None

        # Strategy 1: Try loading from JSON
        if PHASE_A_RESULTS_FILE.exists():
            try:
                with open(PHASE_A_RESULTS_FILE, encoding="utf-8") as f:
                    data: dict[str, Any] = json.load(f)

                molecules_count = data.get("n_molecules", 0)
                results = data.get("results", [])

                # Extract timestamp from results
                last_updated_str = None
                if results:
                    last_result = results[-1]
                    timestamp = last_result.get("timestamp", "")
                    if timestamp:
                        last_updated_str = timestamp[:10]

                if last_updated_str:
                    last_updated = date.fromisoformat(last_updated_str)
                else:
                    # Infer from file modification time
                    try:
                        file_mtime = os.path.getmtime(PHASE_A_RESULTS_FILE)
                        last_updated = date.fromtimestamp(file_mtime)
                    except OSError:
                        last_updated = None

            except (json.JSONDecodeError, OSError, ValueError):
                # JSON exists but is corrupted, will fallback to CSV
                molecules_count = 0

        # Strategy 2: Fallback to CSV if JSON missing or has 0 molecules
        if molecules_count == 0:
            csv_path = DATA_DIR / "thermo_pm7.csv"
            if csv_path.exists():
                try:
                    df = pd.read_csv(csv_path, low_memory=False)

                    # ✅ FIX: Count ONLY molecules with status="OK"
                    # (calculated molecules, not PENDING)
                    if "status" in df.columns:
                        molecules_count = len(df[df["status"] == "OK"])
                    else:
                        # Legacy CSV without status column
                        molecules_count = len(df)

                    # Use CSV modification time if no JSON timestamp
                    if last_updated is None:
                        try:
                            file_mtime = os.path.getmtime(csv_path)
                            last_updated = date.fromtimestamp(file_mtime)
                        except OSError:
                            last_updated = None

                except (pd.errors.EmptyDataError, OSError, ValueError):
                    # CSV exists but is empty or corrupted
                    molecules_count = 0

        # Strategy 3: Return None if no valid data source
        if (
            molecules_count == 0
            and not PHASE_A_RESULTS_FILE.exists()
            and not (DATA_DIR / "thermo_pm7.csv").exists()
        ):
            return None

        # Return Database object with molecule count from JSON or CSV
        return Database(
            name="CREST PM7",
            description="CREST conformer search with PM7 optimization",
            molecules=molecules_count,
            last_updated=last_updated,
            status="ready" if molecules_count > 0 else "in_development",
            properties=["H298_pm7", "conformers", "smiles", "quality_grade"],
        )

    def get_databases(self) -> list[Database]:
        """Get list of databases, replacing CREST PM7 with real data if available.

        Returns:
            List of Database objects with CREST PM7 from real calculations.
        """
        real_pm7 = self.load_real_phase_a_results()
        databases = []

        for db in DATABASES:
            if db.name == "CREST PM7" and real_pm7 is not None:
                databases.append(real_pm7)
            else:
                databases.append(db)

        return databases

    def _format_last_updated(self, value: datetime | date | None) -> str:
        """Format database last_updated timestamp consistently.

        Args:
            value: A datetime, date, or None.

        Returns:
            Formatted string or "-" if value is None.

        Note:
            If the Database.last_updated is timezone-aware, the timezone info
            is preserved in the formatted output.
        """
        if value is None:
            return "-"
        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d %H:%M:%S")
        if isinstance(value, date):
            return value.strftime("%Y-%m-%d")
        return str(value)

    def render(self) -> None:
        """Render the databases overview."""
        self.clear_screen()
        self.show_header()

        # Databases table
        table = Table(
            title="Available Databases",
            show_header=True,
            header_style=f"bold {COLORS['databases']}",
            border_style=COLORS["border"],
        )
        table.add_column("Name", style="bold")
        table.add_column("Molecules", justify="right")
        table.add_column("Last Updated")
        table.add_column("Status")

        for db in self.get_databases():
            if db.status == "ready":
                status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            else:
                status = (
                    f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Dev[/{COLORS['in_dev']}]"
                )

            # Format last_updated consistently
            if isinstance(db.last_updated, datetime):
                last_updated_str = db.last_updated.strftime("%Y-%m-%d %H:%M:%S")
            elif isinstance(db.last_updated, date):
                last_updated_str = db.last_updated.strftime("%Y-%m-%d")
            elif db.last_updated is None:
                last_updated_str = "-"
            else:
                last_updated_str = str(db.last_updated)

            table.add_row(
                db.name,
                f"{db.molecules:,}" if db.molecules > 0 else "-",
                last_updated_str,
                status,
            )

        self.console.print(table)
        self.console.print()

    def render_database_detail(self, db: Database) -> None:
        """Render detailed view for a specific database."""
        self.clear_screen()
        self.show_header()

        # Database info panel
        status_text = (
            f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            if db.status == "ready"
            else f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Development[/{COLORS['in_dev']}]"
        )

        last_updated_str = self._format_last_updated(db.last_updated)

        info = f"""
[bold]Name:[/bold]         {db.name}
[bold]Description:[/bold]  {db.description}
[bold]Molecules:[/bold]    {db.molecules:,}
[bold]Last Updated:[/bold] {last_updated_str}
[bold]Status:[/bold]       {status_text}

[bold]Properties:[/bold]
"""
        for prop in db.properties:
            info += f"  • {prop}\n"

        self.console.print(
            Panel(
                info,
                title=f"[bold {COLORS['databases']}]{db.name}[/bold {COLORS['databases']}]",
                border_style=COLORS["databases"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the databases view."""
        options = []
        for db in self.get_databases():
            options.append(
                MenuOption(
                    label=f"View {db.name}",
                    value=f"view_{db.name}",
                    icon=ICONS["databases"],
                )
            )

        options.extend(
            [
                MenuOption(
                    label="Calculate PM7 Values",
                    value="calculate",
                    disabled=False,
                ),
                MenuOption(
                    label="Refresh Database",
                    value="refresh",
                    icon="",
                    description="Resync working CSV from source",
                ),
                MenuOption(
                    label="Add New Database",
                    value="add",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )

        return options

    def get_detail_menu_options(self) -> list[MenuOption]:
        """Return menu options for database detail view."""
        return [
            MenuOption(
                label="Calculate PM7 Values",
                value="calculate",
                disabled=False,
            ),
            MenuOption(
                label="Edit Database",
                value="edit",
                disabled=True,
                disabled_reason="In Development",
            ),
            MenuOption(
                label="Delete Database",
                value="delete",
                disabled=True,
                disabled_reason="In Development",
            ),
        ]

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions."""
        if action == "back":
            if self.selected_db:
                self.selected_db = None
                return None  # Stay in databases view
            return "main"

        if action and action.startswith("view_"):
            db_name = action.removeprefix("view_")
            found = False
            for db in self.get_databases():
                if db.name == db_name:
                    self.selected_db = db
                    found = True
                    break

            if not found:
                self.selected_db = None

            return None

        if action == "calculate":
            self.handle_calculate_pm7()
            return None

        if action == "refresh":
            self.refresh_databases_from_filesystem()
            self.console.input("[dim]Press Enter to continue...[/dim]")
            return None

        if action in ["add", "edit", "delete"]:
            self.show_in_development(action.title())
            return None

        return None

    def _refresh_database(self) -> None:
        """Trigger database refresh from source.

        Resyncs the working CSV from the source-of-truth CSV,
        resetting all status fields to PENDING.
        """
        from grimperium.cli.dataset_manager import DatasetManager

        manager = DatasetManager(
            source_csv=DATA_DIR / "thermo_cbs_chon.csv",
            working_csv=DATA_DIR / "thermo_pm7.csv",
        )
        try:
            manager.refresh_database()
            self.console.print("[green]Database refreshed[/green]")
        except FileNotFoundError as e:
            self.console.print(f"[red]Refresh failed: {e}[/red]")
        except Exception as e:
            self.console.print(f"[red]Refresh failed: {e}[/red]")
        self.console.input("[dim]Press Enter to continue...[/dim]")

    def refresh_databases_from_filesystem(self) -> int:
        """Scan data/ directory for CSV files and display database info.

        Discovers available database CSV files in the data/ directory,
        counts molecules in each file, and displays summary.

        Returns:
            Number of database CSV files discovered.
        """
        try:
            if not DATA_DIR.exists():
                self.console.print(
                    f"[{COLORS['warning']}]⚠️  Data directory not found: {DATA_DIR}[/{COLORS['warning']}]"
                )
                return 0

            # Find all CSV files in data/ directory
            csv_files = list(DATA_DIR.glob("*.csv"))

            if not csv_files:
                self.console.print(
                    f"[{COLORS['muted']}]No CSV files found in {DATA_DIR}[/{COLORS['muted']}]"
                )
                return 0

            self.console.print()
            self.console.print(
                f"[bold {COLORS['databases']}]Discovered {len(csv_files)} database(s):[/bold {COLORS['databases']}]"
            )
            self.console.print()

            # Display each database with molecule count
            for csv_file in sorted(csv_files):
                try:
                    df = pd.read_csv(csv_file)

                    # ✅ FIX: Show breakdown of calculated (OK) vs pending
                    if "status" in df.columns:
                        ok_count = len(df[df["status"] == "OK"])
                        total_count = len(df)
                        pending_count = total_count - ok_count

                        if ok_count > 0:
                            msg = f"{ok_count} calculated"
                            if pending_count > 0:
                                msg += f" ({pending_count} pending)"
                        else:
                            msg = f"{total_count} pending (none calculated)"
                    else:
                        # Legacy CSV without status - just show total
                        msg = f"{len(df)} molecules"

                    self.console.print(
                        f"  [{COLORS['success']}]✓[/{COLORS['success']}] "
                        f"[bold]{csv_file.name}[/bold]: "
                        f"[{COLORS['databases']}]{msg}[/{COLORS['databases']}]"
                    )
                except Exception as e:
                    self.console.print(
                        f"  [{COLORS['error']}]✗[/{COLORS['error']}] "
                        f"[bold]{csv_file.name}[/bold]: "
                        f"[{COLORS['muted']}]Error reading file: {e}[/{COLORS['muted']}]"
                    )

            self.console.print()
            self.console.print(
                f"[{COLORS['success']}]✓ Refresh complete: {len(csv_files)} database(s) found[/{COLORS['success']}]"
            )

            return len(csv_files)

        except Exception as e:
            self.console.print(
                f"[{COLORS['error']}]❌ Refresh failed: {e}[/{COLORS['error']}]"
            )
            return 0

    def handle_calculate_pm7(self) -> None:
        """Interactive handler for CREST PM7 batch calculations.

        Collects user inputs, runs batch processor with progress bar,
        updates results and returns to database view.
        """
        self.console.print()
        self.console.print(
            Panel(
                "[bold cyan]CREST PM7 Batch Calculation Configuration[/bold cyan]",
                border_style="cyan",
            )
        )
        self.console.print()

        try:
            num_molecules = self._prompt_positive_int(
                "How many molecules to calculate?",
                default=10,
                max_value=30026,
            )
            if num_molecules is None:
                return

            crest_timeout = self._prompt_positive_int(
                "CREST timeout per molecule (minutes)?",
                default=30,
            )
            if crest_timeout is None:
                return

            mopac_timeout = self._prompt_positive_int(
                "MOPAC/PM7 timeout per molecule (minutes)?",
                default=60,
            )
            if mopac_timeout is None:
                return

            self.console.print()
            self.console.print("[bold]Configuration Summary:[/bold]")
            self.console.print(f"  • Molecules to calculate: {num_molecules}")
            self.console.print(f"  • CREST timeout: {crest_timeout} min")
            self.console.print(f"  • MOPAC timeout: {mopac_timeout} min")
            self.console.print()

            confirm = (
                self.console.input(
                    "[yellow]Proceed with calculation? (yes/no) [/yellow]"
                )
                .strip()
                .lower()
            )
            if confirm not in ("yes", "y"):
                self.console.print("[yellow]Calculation cancelled.[/yellow]")
                self.console.input("[dim]Press Enter to continue...[/dim]")
                return

            self._run_pm7_batch(
                num_molecules=num_molecules,
                crest_timeout_minutes=crest_timeout,
                mopac_timeout_minutes=mopac_timeout,
            )

        except KeyboardInterrupt:
            self.console.print("\n[yellow]Calculation interrupted.[/yellow]")
            self.console.input("[dim]Press Enter to continue...[/dim]")

    def _prompt_positive_int(
        self,
        prompt: str,
        default: int,
        max_value: int | None = None,
    ) -> int | None:
        """Prompt user for a positive integer with validation.

        Args:
            prompt: Question to display
            default: Default value if user presses Enter
            max_value: Maximum allowed value (optional)

        Returns:
            Validated integer or None if user cancelled
        """
        while True:
            try:
                value_str = self.console.input(
                    f"[yellow]{prompt} [default: {default}][/yellow] "
                ).strip()

                if not value_str:
                    return default

                value = int(value_str)
                if value <= 0:
                    self.console.print("[red]✗ Must be a positive number[/red]")
                    continue
                if max_value and value > max_value:
                    self.console.print(f"[red]✗ Cannot exceed {max_value}[/red]")
                    continue
                return value

            except ValueError:
                self.console.print("[red]✗ Invalid number. Please try again.[/red]")
            except KeyboardInterrupt:
                return None

    def _run_pm7_batch(
        self,
        num_molecules: int,
        crest_timeout_minutes: int,
        mopac_timeout_minutes: int,
    ) -> None:
        """Execute PM7 batch processing with progress display.

        Args:
            num_molecules: Total molecules to process (also used as batch_size)
            crest_timeout_minutes: CREST timeout per molecule
            mopac_timeout_minutes: MOPAC timeout per molecule
        """
        from grimperium.crest_pm7.batch import (
            BatchCSVManager,
            BatchExecutionManager,
            BatchSortingStrategy,
            ConformerDetailManager,
            FixedTimeoutProcessor,
        )
        from grimperium.crest_pm7.config import PM7Config

        csv_path = DATA_DIR / "thermo_pm7.csv"
        detail_dir = DATA_DIR / "molecules_pm7" / "conformer_details"

        if not csv_path.exists():
            self.console.print(
                f"[red]✗ Batch CSV not found: {csv_path}[/red]\n"
                "[dim]Create a batch CSV with columns: mol_id, smiles, nheavy, status[/dim]"
            )
            self.console.input("[dim]Press Enter to continue...[/dim]")
            return

        try:
            self.console.print()
            self.console.print("[bold cyan]Initializing batch processor...[/bold cyan]")

            pm7_config = PM7Config()
            csv_manager = BatchCSVManager(csv_path)
            csv_manager.load_csv()

            detail_dir.mkdir(parents=True, exist_ok=True)
            detail_manager = ConformerDetailManager(detail_dir)

            processor = FixedTimeoutProcessor(
                config=pm7_config,
                crest_timeout_minutes=float(crest_timeout_minutes),
                mopac_timeout_minutes=float(mopac_timeout_minutes),
            )

            exec_manager = BatchExecutionManager(
                csv_manager=csv_manager,
                detail_manager=detail_manager,
                pm7_config=pm7_config,
                processor_adapter=processor,
            )

            batch_id = csv_manager.generate_batch_id()
            batch = csv_manager.select_batch(
                batch_id=batch_id,
                batch_size=num_molecules,
                crest_timeout_minutes=crest_timeout_minutes,
                mopac_timeout_minutes=mopac_timeout_minutes,
                strategy=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
            )

            if batch.is_empty:
                self.console.print(
                    "[yellow]No molecules available for processing.[/yellow]\n"
                    "[dim]All molecules may already be processed or none are PENDING.[/dim]"
                )
                self.console.input("[dim]Press Enter to continue...[/dim]")
                return

            self.console.print(
                f"[bold cyan]Starting batch {batch_id}: "
                f"{batch.size} molecules[/bold cyan]"
            )
            self.console.print()

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
                        task,
                        completed=current,
                        description=f"Processing {mol_id}",
                    )

                result = exec_manager.execute_batch(
                    batch, progress_callback=update_progress
                )

            self.console.print()
            self.console.print("[bold green]✓ Batch completed![/bold green]")
            self.console.print()
            self.console.print("[bold]Results Summary:[/bold]")
            self.console.print(f"  • Total processed: {result.total_count}")
            self.console.print(f"  • Successful: {result.success_count}")
            self.console.print(f"  • Failed: {result.failed_count}")
            self.console.print(f"  • Skipped: {result.skip_count}")

            if result.total_count > 0:
                rate = result.success_count / result.total_count * 100
                self.console.print(f"  • Success rate: {rate:.1f}%")

            self.console.print()
            self.console.print(
                f"[cyan]Results saved to:[/cyan] {PHASE_A_RESULTS_FILE.parent}"
            )

        except FileNotFoundError as e:
            self.console.print(f"[red]✗ File not found: {e}[/red]")
        except RuntimeError as e:
            self.console.print(f"[red]✗ Runtime error: {e}[/red]")
        except Exception as e:
            self.console.print(f"[red]✗ Unexpected error: {e}[/red]")

        self.console.print()
        self.console.input("[dim]Press Enter to continue...[/dim]")

    def run(self) -> str | None:
        """Run the databases view interaction loop."""
        while True:
            if self.selected_db:
                self.render_database_detail(self.selected_db)
                result = show_back_menu(
                    options=self.get_detail_menu_options(),
                    title="Actions",
                )
            else:
                self.render()
                result = show_back_menu(
                    options=self.get_menu_options(),
                    title="Select Database",
                )

            if result is None or result == "back":  # Ctrl+C
                if self.selected_db:
                    self.selected_db = None
                else:
                    return "main"
            else:
                next_view = self.handle_action(result)
                if next_view:
                    return next_view
