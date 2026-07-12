"""Componentes Rich reutilizaveis para a CLI."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from grimperium.cli.styles import COLORS


@dataclass(frozen=True)
class StatusBadge:
    """Badge textual para estados de execucao e sessao."""

    status: str

    @property
    def color(self) -> str:
        normalized = self.status.lower()
        if normalized in {"completed", "ready", "available"}:
            return COLORS["success"]
        if normalized in {"running", "partial", "created"}:
            return COLORS["warning"]
        if normalized in {"failed", "invalidated", "missing", "incompatible"}:
            return COLORS["error"]
        if normalized == "cancelled":
            return COLORS["muted"]
        return COLORS["highlight"]

    @property
    def text(self) -> str:
        return f"[{self.color}]{self.status}[/{self.color}]"

    def render(self) -> Text:
        return Text(self.status, style=self.color)


@dataclass(frozen=True)
class DetailsTable:
    """Tabela simples de pares label/valor."""

    title: str
    rows: list[tuple[str, object]]
    border_style: str = COLORS["border"]

    def render(self) -> Table:
        table = Table(
            title=self.title,
            show_header=False,
            border_style=self.border_style,
        )
        table.add_column("Field", style="bold")
        table.add_column("Value")
        for label, value in self.rows:
            table.add_row(label, str(value))
        return table


@dataclass(frozen=True)
class EmptyStatePanel:
    """Painel para estados vazios com proximo passo claro."""

    title: str
    message: str
    border_style: str = COLORS["border"]

    def render(self) -> Panel:
        return Panel(
            f"[{COLORS['muted']}]{self.message}[/{COLORS['muted']}]",
            title=self.title,
            border_style=self.border_style,
            padding=(1, 2),
        )


@dataclass(frozen=True)
class ConfirmationSummary:
    """Resumo de confirmacao antes de acoes destrutivas ou longas."""

    title: str
    rows: list[tuple[str, object]]
    border_style: str = COLORS["warning"]

    def render(self) -> Panel:
        return Panel(
            DetailsTable(
                title="",
                rows=self.rows,
                border_style=self.border_style,
            ).render(),
            title=self.title,
            border_style=self.border_style,
            padding=(1, 1),
        )


@dataclass(frozen=True)
class SessionContextPanel:
    """Painel compacto com o contexto cientifico ativo."""

    summary: Mapping[str, str]
    title: str = "Session"

    def render(self) -> Panel:
        table = DetailsTable(
            title="",
            rows=[
                ("Property", self.summary.get("property", "Not selected")),
                ("Method", self.summary.get("method", "Not selected")),
                ("Dataset", self.summary.get("dataset", "Not selected")),
                ("Model", self.summary.get("model", "No model selected")),
                ("Status", StatusBadge(self.summary.get("status", "Unknown")).text),
            ],
        ).render()
        return Panel(
            table,
            title=self.title,
            border_style=COLORS["border"],
            padding=(1, 1),
        )
