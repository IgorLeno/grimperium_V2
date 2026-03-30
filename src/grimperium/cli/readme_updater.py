"""Utilities for updating the project README with live project metrics."""

from __future__ import annotations

import difflib
import pickle
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TypedDict, cast

import joblib
import numpy as np
import pandas as pd
import questionary
from prompt_toolkit.styles import Style
from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table

from grimperium.cli.styles import COLORS
from grimperium.core.metrics import mae, r2_score
from grimperium.ml.gate import evaluate_gate


class PM7Stats(TypedDict):
    """Structured PM7 baseline stats for README rendering."""

    total_ok: int
    total_rows: int
    r2: float | None
    mare: float | None
    p50: float | None
    p90: float | None
    bias: float | None


class ModelStats(TypedDict):
    """Structured trained-model stats for README rendering."""

    model_found: bool
    model_name: str | None
    model_path: str | None
    mae: float | None
    rmse: float | None
    r2: float | None
    gate_pass: bool | None
    n_train: int | None
    n_test: int | None
    date: str | None


README_UPDATER_STYLE = Style.from_dict(
    {
        "qmark": COLORS["settings"],
        "question": "bold",
        "answer": COLORS["success"],
        "pointer": f"{COLORS['settings']} bold",
        "highlighted": f"{COLORS['settings']} bold",
        "selected": COLORS["success"],
        "separator": COLORS["muted"],
        "instruction": COLORS["muted"],
        "text": "",
        "disabled": COLORS["in_dev"],
    }
)


@dataclass
class ReadmeUpdater:
    """Collect live project stats and refresh the dynamic README section."""

    console: Console = field(default_factory=Console)
    readme_path: Path = field(default_factory=lambda: Path.cwd() / "README.md")
    pm7_csv_path: Path = field(default_factory=lambda: Path("data/thermo_pm7.csv"))
    model_bundle_path: Path = field(default_factory=lambda: Path("models/"))

    def _resolve_path(self, path: Path) -> Path:
        """Resolve *path* relative to the current working directory."""
        return path if path.is_absolute() else Path.cwd() / path

    def _empty_pm7_stats(self) -> PM7Stats:
        """Return a PM7 stats payload for missing or unreadable datasets."""
        return {
            "total_ok": 0,
            "total_rows": 0,
            "r2": None,
            "mare": None,
            "p50": None,
            "p90": None,
            "bias": None,
        }

    def _empty_model_stats(self) -> ModelStats:
        """Return a model stats payload for missing model bundles."""
        return {
            "model_found": False,
            "model_name": None,
            "model_path": None,
            "mae": None,
            "rmse": None,
            "r2": None,
            "gate_pass": None,
            "n_train": None,
            "n_test": None,
            "date": None,
        }

    def _get_series(
        self,
        df: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> pd.Series | None:
        """Return the first matching column from *df*, if any."""
        for column in candidates:
            if column in df.columns:
                return cast(pd.Series, df[column])
        return None

    def _get_numeric_series(
        self,
        df: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> pd.Series | None:
        """Return a numeric Series for the first matching column."""
        series = self._get_series(df, candidates)
        if series is None:
            return None
        return pd.to_numeric(series, errors="coerce")

    def _coerce_bool(self, value: object) -> bool | None:
        """Convert heterogeneous truthy values to bool when possible."""
        if isinstance(value, bool):
            return value
        if value is None:
            return None

        if isinstance(value, int | float) and not isinstance(value, bool):
            return bool(value)

        text = str(value).strip().lower()
        if text in {"true", "1", "yes", "y", "sim"}:
            return True
        if text in {"false", "0", "no", "n", "nao", "não"}:
            return False
        return None

    def _coerce_int(self, value: object) -> int | None:
        """Convert *value* to int when possible."""
        if isinstance(value, bool) or value is None:
            return None
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            return int(value)

        text = str(value).strip()
        if not text:
            return None
        try:
            return int(text)
        except ValueError:
            return None

    def _coerce_float(self, value: object) -> float | None:
        """Convert *value* to float when possible."""
        if isinstance(value, bool) or value is None:
            return None
        if isinstance(value, int | float):
            return float(value)

        text = str(value).strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None

    def _format_trained_at(self, value: object) -> str | None:
        """Format a stored training timestamp into a short human-readable string."""
        if value is None:
            return None

        text = str(value).strip()
        if not text:
            return None

        try:
            trained_at = datetime.fromisoformat(text)
        except ValueError:
            return text

        if trained_at.tzinfo is None:
            return trained_at.strftime("%Y-%m-%d %H:%M")
        return trained_at.strftime("%Y-%m-%d %H:%M %Z").strip()

    def _load_model_bundle(self, path: Path) -> object:
        """Load a serialized model bundle from *path*."""
        try:
            return joblib.load(path)
        except Exception:
            with path.open("rb") as handle:
                return pickle.load(handle)

    def _extract_model_metric(
        self,
        metrics: object,
        section: str,
        key: str,
    ) -> object:
        """Safely extract nested metrics from a bundle payload."""
        if not isinstance(metrics, dict):
            return None

        section_value = metrics.get(section)
        if not isinstance(section_value, dict):
            return None

        return section_value.get(key)

    def _discover_model_candidates(self) -> list[Path]:
        """Return model files ordered from newest to oldest."""
        base_path = self._resolve_path(self.model_bundle_path)
        if base_path.is_file():
            return [base_path]
        if not base_path.exists():
            return []

        patterns = ("*.joblib", "*.pkl")
        candidates: list[Path] = []
        for pattern in patterns:
            candidates.extend(base_path.rglob(pattern))

        return sorted(
            [candidate for candidate in candidates if candidate.is_file()],
            key=lambda candidate: candidate.stat().st_mtime,
            reverse=True,
        )

    def _format_metric(self, value: float | None, precision: int) -> str:
        """Render a float metric or placeholder for Markdown tables."""
        if value is None:
            return "-"
        return f"{value:.{precision}f}"

    def _format_signed_metric(self, value: float | None, precision: int) -> str:
        """Render a signed float metric or placeholder for Markdown tables."""
        if value is None:
            return "-"
        return f"{value:+.{precision}f}"

    def _format_int(self, value: int | None) -> str:
        """Render an integer metric or placeholder for Markdown tables."""
        if value is None:
            return "-"
        return f"{value:,}"

    def _pm7_variation_text(self, r2_value: float | None) -> str:
        """Return the PM7 interpretation text for the R² row."""
        if r2_value is None:
            return "Correlação indisponível por falta de dados válidos"
        return f"PM7 captura {r2_value * 100:.1f}% da variação real"

    def _bias_text(self, bias_value: float | None) -> str:
        """Return a human-readable bias interpretation."""
        if bias_value is None:
            return "Viés indisponível por falta de dados válidos"
        if bias_value < 0:
            return "PM7 subestima sistematicamente"
        if bias_value > 0:
            return "PM7 superestima sistematicamente"
        return "Sem viés sistemático relevante"

    def _build_diff(self, before: str, after: str) -> str:
        """Return a unified diff string between two text blobs."""
        diff_lines = difflib.unified_diff(
            before.splitlines(),
            after.splitlines(),
            fromfile="README.md (atual)",
            tofile="README.md (novo)",
            lineterm="",
        )
        return "\n".join(diff_lines)

    def collect_pm7_stats(self) -> PM7Stats:
        """Collect PM7 baseline metrics from ``data/thermo_pm7.csv``.

        Returns:
            PM7Stats with row counts and absolute-error metrics. If the CSV is
            missing, empty, or lacks required metric columns, count fields are
            still returned when possible and metric fields remain ``None``.
        """
        csv_path = self._resolve_path(self.pm7_csv_path)
        if not csv_path.exists():
            return self._empty_pm7_stats()

        try:
            df = pd.read_csv(csv_path, low_memory=False)
        except (FileNotFoundError, pd.errors.EmptyDataError):
            return self._empty_pm7_stats()

        if df.empty:
            return self._empty_pm7_stats()

        stats = self._empty_pm7_stats()
        stats["total_rows"] = int(len(df.index))

        status_series = self._get_series(df, ("status",))
        if status_series is not None:
            status_mask = (
                status_series.fillna("").astype(str).str.strip().str.upper().eq("OK")
            )
            stats["total_ok"] = int(status_mask.sum())
            metric_df = df.loc[status_mask].copy()
        else:
            metric_df = df.copy()

        flag_series = self._get_series(metric_df, ("cbs_quality_flag",))
        if flag_series is not None:
            suspect_mask = (
                flag_series.fillna("").astype(str).str.strip().str.upper().eq("SUSPECT")
            )
            metric_df = metric_df.loc[~suspect_mask].copy()

        pm7_series = self._get_numeric_series(metric_df, ("H298_PM7", "H298_pm7"))
        cbs_series = self._get_numeric_series(metric_df, ("H298_cbs", "H298_CBS"))
        if pm7_series is None or cbs_series is None:
            return stats

        valid_mask = pm7_series.notna() & cbs_series.notna()
        if not bool(valid_mask.any()):
            return stats

        pm7_values = pm7_series.loc[valid_mask].to_numpy(dtype=float)
        cbs_values = cbs_series.loc[valid_mask].to_numpy(dtype=float)

        abs_err = np.abs(pm7_values - cbs_values)
        stats["r2"] = float(r2_score(cbs_values.tolist(), pm7_values.tolist()))
        stats["mare"] = float(mae(cbs_values.tolist(), pm7_values.tolist()))
        stats["p50"] = float(np.percentile(abs_err, 50))
        stats["p90"] = float(np.percentile(abs_err, 90))
        stats["bias"] = float(np.mean(pm7_values - cbs_values))
        return stats

    def collect_model_stats(self) -> ModelStats:
        """Collect metrics from the newest trained-model bundle, if present."""
        for candidate in self._discover_model_candidates():
            try:
                bundle = self._load_model_bundle(candidate)
            except Exception:
                continue

            if not isinstance(bundle, dict):
                continue

            metrics = bundle.get("metrics")
            model_r2 = self._coerce_float(
                self._extract_model_metric(metrics, "test", "r2")
            )
            model_mae = self._coerce_float(
                self._extract_model_metric(metrics, "test", "mae")
            )
            model_rmse = self._coerce_float(
                self._extract_model_metric(metrics, "test", "rmse")
            )
            gate_pass = self._coerce_bool(
                self._extract_model_metric(metrics, "test", "gate_pass")
            )

            if (
                gate_pass is None
                and model_r2 is not None
                and model_mae is not None
                and model_rmse is not None
            ):
                gate_pass = evaluate_gate(
                    {
                        "r2": model_r2,
                        "mae": model_mae,
                        "rmse": model_rmse,
                    }
                )

            return {
                "model_found": True,
                "model_name": candidate.name,
                "model_path": str(candidate),
                "mae": model_mae,
                "rmse": model_rmse,
                "r2": model_r2,
                "gate_pass": gate_pass,
                "n_train": self._coerce_int(bundle.get("n_train")),
                "n_test": self._coerce_int(bundle.get("n_test")),
                "date": self._format_trained_at(bundle.get("trained_at")),
            }

        return self._empty_model_stats()

    def render_results_table(self, pm7: PM7Stats, model: ModelStats) -> str:
        """Render the full Markdown block for ``## Resultados Atuais``."""
        pm7_lines = [
            "## Resultados Atuais",
            "",
            "### Baseline PM7",
            "",
            "| Métrica | Valor | O que significa |",
            "|---|---|---|",
            (
                "| Moléculas calculadas (status OK) "
                f"| {self._format_int(pm7['total_ok'])} "
                "| Base de dados em crescimento ativo para treinar a correção ML |"
            ),
            (
                "| Correlação PM7 vs CBS-QB3 (R²) "
                f"| {self._format_metric(pm7['r2'], 4)} "
                f"| {self._pm7_variation_text(pm7['r2'])} |"
            ),
            (
                "| Erro absoluto mediano (P50) "
                f"| {self._format_metric(pm7['p50'], 2)} kcal/mol "
                "| Erro típico que o modelo ML precisará corrigir |"
            ),
            (
                "| Erro no percentil 90 (P90) "
                f"| {self._format_metric(pm7['p90'], 2)} kcal/mol "
                "| Pior caso frequente |"
            ),
            (
                "| Erro absoluto médio (MARE) "
                f"| {self._format_metric(pm7['mare'], 2)} kcal/mol "
                "| Média geral dos erros |"
            ),
            (
                "| Viés médio (bias) "
                f"| {self._format_signed_metric(pm7['bias'], 2)} kcal/mol "
                f"| {self._bias_text(pm7['bias'])} |"
            ),
            "",
        ]

        if not model["model_found"]:
            pm7_lines.append(
                "> Nenhum modelo ML treinado ainda. Execute 'Treinar Modelo' na CLI."
            )
            return "\n".join(pm7_lines)

        gate_text = "-"
        if model["gate_pass"] is True:
            gate_text = "✅ Sim"
        elif model["gate_pass"] is False:
            gate_text = "❌ Não"

        pm7_lines.extend(
            [
                "### Modelo ML Treinado",
                "",
                "| Métrica | Valor | O que significa |",
                "|---|---|---|",
                (
                    "| R² do modelo "
                    f"| {self._format_metric(model['r2'], 4)} "
                    "| Quanto da variação do CBS-QB3 é explicada após a correção ML |"
                ),
                (
                    "| MAE do modelo "
                    f"| {self._format_metric(model['mae'], 2)} kcal/mol "
                    "| Erro médio absoluto após aplicar a correção delta |"
                ),
                (
                    "| RMSE do modelo "
                    f"| {self._format_metric(model['rmse'], 2)} kcal/mol "
                    "| Penaliza erros grandes e resume a robustez global |"
                ),
                (
                    "| Gate pass "
                    f"| {gate_text} "
                    "| Modelo aprovado nos critérios de qualidade |"
                ),
                (
                    "| Moléculas de treino "
                    f"| {self._format_int(model['n_train'])} "
                    "| Tamanho efetivo do conjunto usado no último treinamento |"
                ),
                (
                    "| Data do treino "
                    f"| {model['date'] or '-'} "
                    "| Momento em que o bundle atual foi gerado |"
                ),
            ]
        )
        return "\n".join(pm7_lines)

    def update_readme(self, dry_run: bool = False) -> tuple[bool, str]:
        """Update only the dynamic ``## Resultados Atuais`` block in the README.

        Args:
            dry_run: When ``True``, return the updated content without writing it.

        Returns:
            Tuple of ``(success, payload)`` where payload is either the updated
            README content (dry-run) or a status/error message.
        """
        readme_path = self._resolve_path(self.readme_path)
        if not readme_path.exists():
            return False, f"README não encontrado: {readme_path}"

        try:
            original = readme_path.read_text(encoding="utf-8")
        except OSError as exc:
            return False, f"Falha ao ler README.md: {exc}"

        section_pattern = re.compile(
            r"(^## Resultados Atuais\n)(.*?)(\n---\n)",
            flags=re.MULTILINE | re.DOTALL,
        )
        match = section_pattern.search(original)
        if match is None:
            return False, "Seção '## Resultados Atuais' não encontrada no README.md"

        pm7_stats = self.collect_pm7_stats()
        model_stats = self.collect_model_stats()
        new_section = self.render_results_table(pm7_stats, model_stats)
        separator = "\n---\n"
        updated = (
            original[: match.start()]
            + new_section
            + original[match.end() - len(separator) :]
        )

        if dry_run:
            return True, updated

        try:
            readme_path.write_text(updated, encoding="utf-8")
        except OSError as exc:
            return False, f"Falha ao escrever README.md: {exc}"

        return True, "README atualizado"

    def _show_pm7_table(self, stats: PM7Stats) -> None:
        """Display current PM7 statistics in a Rich table."""
        table = Table(
            title="Estatísticas PM7",
            show_header=True,
            header_style=f"bold {COLORS['settings']}",
            border_style=COLORS["settings"],
        )
        table.add_column("Métrica", style="bold")
        table.add_column("Valor", justify="right")

        table.add_row("Linhas lidas", self._format_int(stats["total_rows"]))
        table.add_row("Status OK", self._format_int(stats["total_ok"]))
        table.add_row("R²", self._format_metric(stats["r2"], 4))
        table.add_row("MARE", self._format_metric(stats["mare"], 2))
        table.add_row("P50", self._format_metric(stats["p50"], 2))
        table.add_row("P90", self._format_metric(stats["p90"], 2))
        bias_text = "-"
        if stats["bias"] is not None:
            bias_text = f"{stats['bias']:+.2f}"
        table.add_row("Bias", bias_text)
        self.console.print(table)

    def _show_model_table(self, stats: ModelStats) -> None:
        """Display current trained-model statistics in a Rich table."""
        if not stats["model_found"]:
            self.console.print(
                Panel(
                    "Nenhum bundle de modelo treinado foi encontrado em models/.",
                    title="[bold]Modelo ML[/bold]",
                    border_style=COLORS["settings"],
                )
            )
            return

        table = Table(
            title=f"Modelo treinado: {stats['model_name']}",
            show_header=True,
            header_style=f"bold {COLORS['settings']}",
            border_style=COLORS["settings"],
        )
        table.add_column("Métrica", style="bold")
        table.add_column("Valor", justify="right")

        table.add_row("R²", self._format_metric(stats["r2"], 4))
        table.add_row("MAE", self._format_metric(stats["mae"], 2))
        table.add_row("RMSE", self._format_metric(stats["rmse"], 2))
        if stats["gate_pass"] is True:
            gate_text = "✅ Sim"
        elif stats["gate_pass"] is False:
            gate_text = "❌ Não"
        else:
            gate_text = "-"
        table.add_row("Gate pass", gate_text)
        table.add_row("Moléculas de treino", self._format_int(stats["n_train"]))
        table.add_row("Moléculas de teste", self._format_int(stats["n_test"]))
        table.add_row("Data do treino", stats["date"] or "-")
        self.console.print(table)

    def display_menu(self) -> None:
        """Display the README updater menu and handle user actions."""
        while True:
            self.console.clear()
            self.console.print()
            self.console.print(
                Panel(
                    "[dim]Atualiza o README.md com os dados reais[/dim]",
                    title="[bold]📊 Atualizar Documentação do Projeto[/bold]",
                    border_style=COLORS["settings"],
                )
            )

            choice = questionary.select(
                "Selecione uma ação:",
                choices=[
                    questionary.Choice(
                        "📈 Ver estatísticas PM7 atuais", value="pm7_stats"
                    ),
                    questionary.Choice(
                        "🤖 Ver métricas do modelo treinado", value="model_stats"
                    ),
                    questionary.Choice(
                        "📝 Prévia das mudanças (dry run)", value="dry_run"
                    ),
                    questionary.Choice("✅ Atualizar README.md agora", value="write"),
                    questionary.Separator(),
                    questionary.Choice("◀ Voltar", value="back"),
                ],
                style=README_UPDATER_STYLE,
                qmark="",
                pointer="❯",
            ).ask()

            if choice is None or choice == "back":
                return

            self.console.clear()
            self.console.print()

            if choice == "pm7_stats":
                self._show_pm7_table(self.collect_pm7_stats())
                self.console.print()
                self.console.input("[dim]Press Enter to continue...[/dim]")
                continue

            if choice == "model_stats":
                self._show_model_table(self.collect_model_stats())
                self.console.print()
                self.console.input("[dim]Press Enter to continue...[/dim]")
                continue

            if choice == "dry_run":
                readme_path = self._resolve_path(self.readme_path)
                original = (
                    readme_path.read_text(encoding="utf-8")
                    if readme_path.exists()
                    else ""
                )
                success, payload = self.update_readme(dry_run=True)
                if not success:
                    self.console.print(
                        f"[{COLORS['error']}]{payload}[/{COLORS['error']}]"
                    )
                else:
                    diff_text = self._build_diff(original, payload)
                    if not diff_text:
                        self.console.print(
                            f"[{COLORS['warning']}]Nenhuma mudança detectada no README.[/{COLORS['warning']}]"
                        )
                    else:
                        self.console.print(
                            Syntax(
                                diff_text, "diff", theme="ansi_dark", word_wrap=False
                            )
                        )
                self.console.print()
                self.console.input("[dim]Press Enter to continue...[/dim]")
                continue

            confirmed = questionary.confirm(
                "Isso vai reescrever o README.md. Confirmar?",
                default=False,
                style=README_UPDATER_STYLE,
                qmark="",
            ).ask()
            if not confirmed:
                self.console.print(
                    f"[{COLORS['muted']}]Atualização cancelada.[/{COLORS['muted']}]"
                )
                self.console.input("[dim]Press Enter to continue...[/dim]")
                continue

            success, payload = self.update_readme(dry_run=False)
            color = COLORS["success"] if success else COLORS["error"]
            self.console.print(f"[{color}]{payload}[/{color}]")
            self.console.print()
            self.console.input("[dim]Press Enter to continue...[/dim]")
