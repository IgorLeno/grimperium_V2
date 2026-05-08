"""Rich TUI application for grimperium_mini."""

from __future__ import annotations

import sys
from pathlib import Path

from rich.live import Live
from rich.panel import Panel
from rich.prompt import Prompt

from .config import MiniConfig
from .io import load_pipeline_tasks
from .pipeline import run_pipeline
from .progress import MiniProgressTracker
from .styles import COLORS, MINI_BANNER, MINI_SUBTITLE, MINI_VERSION, console


class MiniApp:
    """Interactive TUI application that mirrors GrimperiumCLI structure."""

    def __init__(self) -> None:
        self.console = console
        self._running = True

    def show_welcome(self) -> None:
        self.console.clear()
        self.console.print()
        banner_content = f"{MINI_BANNER}\n{MINI_SUBTITLE}\n{MINI_VERSION}"
        self.console.print(
            Panel(
                banner_content,
                border_style=COLORS["primary"],
                padding=(1, 4),
            )
        )
        self.console.print()

    def display_main_menu(self) -> str | None:
        self.console.print(
            Panel(
                "[bold]MAIN MENU[/bold]\n\n"
                "  [cyan]1[/cyan] · Executar pipeline\n"
                "  [cyan]2[/cyan] · Validar dados (--xlsx)\n"
                "  [cyan]3[/cyan] · Exportar resumo\n"
                "  [cyan]Q[/cyan] · Sair",
                border_style=COLORS["border"],
                padding=(0, 2),
            )
        )
        self.console.print()
        choice = Prompt.ask(
            "[bold cyan]Escolha[/bold cyan]",
            choices=["1", "2", "3", "q", "Q"],
            default="Q",
        ).lower()
        if choice == "q":
            return None
        mapping = {"1": "run", "2": "validate", "3": "export"}
        return mapping.get(choice)

    def run_pipeline_interactive(self) -> None:
        default_xlsx = "data/grimperium_mini_pipeline_tcc.xlsx"
        xlsx_str = Prompt.ask("Caminho do xlsx", default=default_xlsx)
        xlsx = Path(xlsx_str)

        methods_str = Prompt.ask("Métodos (espaço separado)", default="AM1 PM3 PM7")
        methods = [m.upper() for m in methods_str.split()]

        limit_str = Prompt.ask("Limit (Enter para processar tudo)", default="")
        limit: int | None = int(limit_str) if limit_str.strip() else None

        dry_str = Prompt.ask("Dry-run?", choices=["s", "n"], default="n")
        dry_run = dry_str.lower() == "s"

        config = MiniConfig()
        tasks = load_pipeline_tasks(xlsx, methods=methods, limit=limit)
        total = len(tasks)

        tracker = MiniProgressTracker(total_tasks=total)

        def on_progress(event: str, data: object) -> None:
            if event == "start" and isinstance(data, dict):
                tracker.on_task_start(
                    str(data.get("mol_id", "")), str(data.get("method", ""))
                )
            elif event == "done" and isinstance(data, dict):
                tracker.on_task_done(data)

        self.console.print()
        with Live(tracker.render(), refresh_per_second=4, console=self.console) as live:

            def _progress(event: str, data: object) -> None:
                on_progress(event, data)
                live.update(tracker.render())

            run_pipeline(
                xlsx,
                methods=methods,
                config=config,
                limit=limit,
                dry_run=dry_run,
                on_progress=_progress,
            )

        self.console.print(tracker.summary())

    def run_validate_interactive(self) -> None:
        from .io import validate_workbook

        default_xlsx = "data/grimperium_mini_pipeline_tcc.xlsx"
        xlsx_str = Prompt.ask("Caminho do xlsx", default=default_xlsx)
        xlsx = Path(xlsx_str)
        issues = validate_workbook(xlsx)
        if issues:
            self.console.print(f"[{COLORS['error']}]INVALID[/{COLORS['error']}]")
            for issue in issues:
                self.console.print(
                    f"  [{COLORS['error']}]- {issue}[/{COLORS['error']}]"
                )
        else:
            self.console.print(f"[{COLORS['success']}]OK: {xlsx}[/{COLORS['success']}]")

    def run_export_interactive(self) -> None:
        from .io import read_csv_rows, write_simple_xlsx

        input_str = Prompt.ask("CSV de formação (input)")
        reactions_str = Prompt.ask("CSV de reações")
        output_str = Prompt.ask(
            "Arquivo xlsx de saída", default="grimperium_mini_summary.xlsx"
        )

        input_path = Path(input_str)
        reactions_path = Path(reactions_str)
        output_path = Path(output_str)

        formation = read_csv_rows(input_path)
        reactions = read_csv_rows(reactions_path)
        write_simple_xlsx(
            output_path,
            {
                "formacao_grimperium_mini": (
                    list(formation[0]) if formation else [],
                    [dict(row) for row in formation],
                ),
                "reacoes_grimperium_mini": (
                    list(reactions[0]) if reactions else [],
                    [dict(row) for row in reactions],
                ),
            },
        )
        self.console.print(
            f"[{COLORS['success']}]Wrote {output_path}[/{COLORS['success']}]"
        )

    def run(self) -> int:
        try:
            while self._running:
                self.show_welcome()
                selection = self.display_main_menu()

                if selection is None:
                    break

                if selection == "run":
                    self.run_pipeline_interactive()
                    Prompt.ask("\nPressione Enter para continuar", default="")
                elif selection == "validate":
                    self.run_validate_interactive()
                    Prompt.ask("\nPressione Enter para continuar", default="")
                elif selection == "export":
                    self.run_export_interactive()
                    Prompt.ask("\nPressione Enter para continuar", default="")

        except KeyboardInterrupt:
            self.console.print()
            self.console.print(
                f"[{COLORS['muted']}]Interrupted by user[/{COLORS['muted']}]"
            )

        except Exception as exc:
            self.console.print()
            self.console.print(
                f"[{COLORS['error']}]Unexpected error: {exc}[/{COLORS['error']}]"
            )
            return 1

        return 0


def main() -> int:
    return MiniApp().run()


if __name__ == "__main__":
    sys.exit(main())
