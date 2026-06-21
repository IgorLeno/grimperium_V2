"""Rich TUI application for grimperium_mini."""

from __future__ import annotations

import sys
from pathlib import Path

from rich.live import Live
from rich.panel import Panel
from rich.prompt import Prompt

from .config import MiniConfig
from .io import load_molecules
from .multi_conformer import run_multi_conformer_pipeline
from .pipeline import run_pipeline
from .progress import MiniProgressTracker
from .styles import COLORS, MINI_BANNER, MINI_SUBTITLE, MINI_VERSION, console

# Bundled TCC dataset shipped with the package (packages/grimperium-mini/data/).
# Resolved relative to this file so the default works regardless of CWD.
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_XLSX_PATH = str(_PACKAGE_ROOT / "data" / "grimperium_mini_pipeline_tcc.xlsx")


class MiniApp:
    """Interactive TUI application that mirrors GrimperiumCLI structure."""

    def __init__(self) -> None:
        self.console = console
        self._running = True
        self._config = MiniConfig()

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
                "  [cyan]4[/cyan] · Configurações\n"
                "  [cyan]5[/cyan] · Multi-conformer mode\n"
                "  [cyan]Q[/cyan] · Sair",
                border_style=COLORS["border"],
                padding=(0, 2),
            )
        )
        self.console.print()
        choice = Prompt.ask(
            "[bold cyan]Escolha[/bold cyan]",
            choices=["1", "2", "3", "4", "5", "q", "Q"],
            default="Q",
        ).lower()
        if choice == "q":
            return None
        mapping = {
            "1": "run",
            "2": "validate",
            "3": "export",
            "4": "settings",
            "5": "multiconf",
        }
        return mapping.get(choice)

    def run_pipeline_interactive(self) -> None:
        default_xlsx = DEFAULT_XLSX_PATH
        xlsx_str = Prompt.ask("Caminho do xlsx", default=default_xlsx)
        xlsx = Path(xlsx_str)

        limit_str = Prompt.ask("Limit (Enter para processar tudo)", default="")
        limit: int | None = int(limit_str) if limit_str.strip() else None

        dry_str = Prompt.ask("Dry-run?", choices=["s", "n"], default="n")
        dry_run = dry_str.lower() == "s"

        config = self._config
        molecules = load_molecules(xlsx, limit=limit)
        total = len(molecules)

        tracker = MiniProgressTracker(total_tasks=total)

        def on_progress(event: str, data: object) -> None:
            if event == "start" and isinstance(data, dict):
                tracker.on_task_start(str(data.get("mol_id", "")), "")
            elif event == "done" and isinstance(data, dict):
                tracker.on_task_done(data)

        self.console.print()
        with Live(tracker.render(), refresh_per_second=4, console=self.console) as live:

            def _progress(event: str, data: object) -> None:
                on_progress(event, data)
                live.update(tracker.render())

            run_pipeline(
                xlsx,
                config=config,
                limit=limit,
                dry_run=dry_run,
                on_progress=_progress,
            )

        self.console.print(tracker.summary())

    def run_settings_interactive(self) -> None:
        """Tela de configurações: editar timeout e max conformers."""
        self.console.print()
        self.console.print(
            Panel(
                "[bold]CONFIGURAÇÕES DO CREST[/bold]",
                border_style=COLORS["primary"],
                padding=(0, 2),
            )
        )
        self.console.print()

        current_timeout_min = self._config.timeout_crest_s // 60
        self.console.print(
            f"  Timeout atual: [{COLORS['success']}]{current_timeout_min} minutos"
            f"[/{COLORS['success']}]  "
            f"([{COLORS['muted']}]{self._config.timeout_crest_s}s[/{COLORS['muted']}])"
        )
        timeout_str = Prompt.ask(
            "  Novo timeout [bold](minutos)[/bold] — Enter para manter",
            default=str(current_timeout_min),
        )
        try:
            new_timeout_min = int(timeout_str.strip())
            if new_timeout_min <= 0:
                raise ValueError
            self._config.timeout_crest_s = new_timeout_min * 60
            self.console.print(
                f"  [{COLORS['success']}]✓ Timeout definido: "
                f"{new_timeout_min} min ({self._config.timeout_crest_s}s)"
                f"[/{COLORS['success']}]"
            )
        except ValueError:
            self.console.print(
                f"  [{COLORS['muted']}]Valor inválido — timeout não alterado."
                f"[/{COLORS['muted']}]"
            )

        self.console.print()

        self.console.print(
            f"  Max conformers atual: [{COLORS['success']}]"
            f"{self._config.crest_max_structures}[/{COLORS['success']}]"
        )
        mstruct_str = Prompt.ask(
            "  Novo máximo de conformers — Enter para manter",
            default=str(self._config.crest_max_structures),
        )
        try:
            new_mstruct = int(mstruct_str.strip())
            if new_mstruct <= 0:
                raise ValueError
            self._config.crest_max_structures = new_mstruct
            self.console.print(
                f"  [{COLORS['success']}]✓ Max conformers definido: "
                f"{new_mstruct}[/{COLORS['success']}]"
            )
        except ValueError:
            self.console.print(
                f"  [{COLORS['muted']}]Valor inválido — max conformers não alterado."
                f"[/{COLORS['muted']}]"
            )

        self.console.print()
        self.console.print(
            Panel(
                f"  Timeout:       [{COLORS['primary']}]"
                f"{self._config.timeout_crest_s // 60} min[/{COLORS['primary']}]\n"
                f"  Max conformers:[{COLORS['primary']}]"
                f" {self._config.crest_max_structures}[/{COLORS['primary']}]",
                title="[bold]Configurações ativas[/bold]",
                border_style=COLORS["border"],
            )
        )

    def run_multiconf_interactive(self) -> None:
        """Interactive multi-conformer mode: reuses CREST files, runs MOPAC×3."""
        self.console.print()
        self.console.print(
            Panel(
                "[bold]MULTI-CONFORMER MODE[/bold]\n\n"
                "Reutiliza os arquivos CREST já calculados e roda MOPAC\n"
                "(AM1 + PM3 + PM7) para os N melhores conformeros.",
                border_style=COLORS["primary"],
                padding=(0, 2),
            )
        )
        self.console.print()

        default_xlsx = DEFAULT_XLSX_PATH
        xlsx_str = Prompt.ask("Caminho do xlsx", default=default_xlsx)
        xlsx = Path(xlsx_str)

        limit_str = Prompt.ask("Limit (Enter para processar tudo)", default="")
        limit: int | None = int(limit_str) if limit_str.strip() else None

        max_conf_str = Prompt.ask("Máximo de conformeros por molécula", default="10")
        try:
            max_conformers = int(max_conf_str.strip())
            if max_conformers <= 0:
                raise ValueError
        except ValueError:
            self.console.print(
                f"  [{COLORS['muted']}]Valor inválido — usando 10.[/{COLORS['muted']}]"
            )
            max_conformers = 10

        config = self._config
        molecules = load_molecules(xlsx, limit=limit)
        total = len(molecules)

        tracker = MiniProgressTracker(total_tasks=total)

        self.console.print()
        with Live(tracker.render(), refresh_per_second=4, console=self.console) as live:

            def _progress(event: str, data: object) -> None:
                if event == "start" and isinstance(data, dict):
                    tracker.on_task_start(str(data.get("mol_id", "")), "")
                elif event == "done" and isinstance(data, dict):
                    tracker.on_task_done(data)
                live.update(tracker.render())

            run_multi_conformer_pipeline(
                xlsx,
                config=config,
                limit=limit,
                max_conformers=max_conformers,
                on_progress=_progress,
            )

        self.console.print(tracker.summary())
        output = config.results_dir / "grimperium_mini_multiconf_summary.csv"
        self.console.print(
            f"\n[{COLORS['success']}]✓ Resultados gravados em {output}[/{COLORS['success']}]"
        )

    def run_validate_interactive(self) -> None:
        from .io import validate_workbook

        default_xlsx = DEFAULT_XLSX_PATH
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
                elif selection == "settings":
                    self.run_settings_interactive()
                    Prompt.ask("\nPressione Enter para continuar", default="")
                elif selection == "multiconf":
                    self.run_multiconf_interactive()
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
