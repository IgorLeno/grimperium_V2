"""Rich progress tracker for grimperium_mini TUI."""

from __future__ import annotations

from rich.panel import Panel
from rich.table import Table

from .styles import COLORS


class MiniProgressTracker:
    """Tracks pipeline progress and renders a Rich Table for Live display."""

    def __init__(self, total_tasks: int) -> None:
        self.total_tasks = total_tasks
        self.rows: list[dict[str, object]] = []
        self._current_mol_id: str = ""
        self.success_count = 0
        self.failed_count = 0
        self.skipped_count = 0

    def on_task_start(self, mol_id: str, _method: str) -> None:
        self._current_mol_id = mol_id

    def on_task_done(self, row: dict[str, object]) -> None:
        self.rows.append(row)
        status = str(row.get("status", ""))
        if status in {"success", "dry_run"}:
            self.success_count += 1
        elif status == "skipped":
            self.skipped_count += 1
        else:
            self.failed_count += 1

    def render(self) -> Table:
        table = Table(
            title="Pipeline Progress",
            border_style=COLORS["border"],
            show_lines=False,
            expand=True,
        )
        table.add_column("mol_id", style="bold", no_wrap=True)
        table.add_column("conformers", justify="right")
        table.add_column("status")
        table.add_column("hf_PM7_kJ/mol", justify="right")
        table.add_column("tempo(s)", justify="right")
        table.add_column("crest")
        table.add_column("mopac")

        if self._current_mol_id:
            done_ids = {str(r.get("mol_id", "")) for r in self.rows}
            if self._current_mol_id not in done_ids:
                table.add_row(
                    f"[{COLORS['primary']}]⟳ {self._current_mol_id}[/{COLORS['primary']}]",
                    "...",
                    f"[{COLORS['primary']}]calculando...[/{COLORS['primary']}]",
                    "...",
                    "...",
                    "...",
                    "...",
                )

        for row in self.rows:
            status = str(row.get("status", ""))
            if status == "skipped":
                table.add_row(
                    str(row.get("mol_id", "")),
                    "—",
                    f"[{COLORS['muted']}]✓ já calculada[/{COLORS['muted']}]",
                    "—",
                    "—",
                    "—",
                    "—",
                )
                continue

            if status == "success":
                status_style = COLORS["success"]
            elif status == "dry_run":
                status_style = COLORS["muted"]
            else:
                status_style = COLORS["error"]

            hof = row.get("hf_pm7_kJmol", "")
            hof_str = f"{hof:.2f}" if isinstance(hof, float) else str(hof)

            elapsed = row.get("total_time_s", "")
            elapsed_str = (
                f"{elapsed:.1f}" if isinstance(elapsed, int | float) else str(elapsed)
            )

            table.add_row(
                str(row.get("mol_id", "")),
                str(row.get("n_conformers_founded", "")),
                f"[{status_style}]{status}[/{status_style}]",
                hof_str,
                elapsed_str,
                str(row.get("crest_status", "")),
                str(row.get("mopac_status", "")),
            )

        done = self.success_count + self.failed_count
        remaining = self.total_tasks - done
        table.caption = (
            f"{done}/{self.total_tasks} concluídas  |  "
            f"[{COLORS['primary']}]{remaining} restantes[/{COLORS['primary']}]"
        )
        return table

    def summary(self) -> Panel:
        total = self.success_count + self.failed_count + self.skipped_count
        body = (
            f"[{COLORS['success']}]{self.success_count} ✓[/{COLORS['success']}]  "
            f"[{COLORS['error']}]{self.failed_count} ✗[/{COLORS['error']}]  "
            f"[{COLORS['muted']}]{self.skipped_count} ↷ puladas[/{COLORS['muted']}]  "
            f"[{COLORS['muted']}]{total} total[/{COLORS['muted']}]"
        )
        return Panel(body, title="[bold]Summary[/bold]", border_style=COLORS["primary"])
