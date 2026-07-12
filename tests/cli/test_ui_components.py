from __future__ import annotations

from grimperium.cli.components import (
    ConfirmationSummary,
    DetailsTable,
    EmptyStatePanel,
    SessionContextPanel,
    StatusBadge,
)
from grimperium.cli.menu import format_session_header


def test_status_badge_formats_known_statuses() -> None:
    assert "completed" in StatusBadge("completed").text
    assert "failed" in StatusBadge("failed").text
    assert "unknown" in StatusBadge("unknown").text


def test_details_table_preserves_labels() -> None:
    table = DetailsTable(
        title="Run",
        rows=[("Run ID", "run_1"), ("Status", "completed")],
    ).render()

    assert table.title == "Run"
    assert len(table.rows) == 2


def test_empty_state_and_confirmation_summary_are_renderable() -> None:
    empty = EmptyStatePanel(
        title="No Results",
        message="Select a dataset or run a calculation.",
    ).render()
    summary = ConfirmationSummary(
        title="Confirm Run",
        rows=[("Method", "PM7"), ("Molecules", "2")],
    ).render()

    assert empty.title == "No Results"
    assert summary.title == "Confirm Run"


def test_session_context_panel_accepts_controller_summary() -> None:
    panel = SessionContextPanel(
        {
            "property": "Standard enthalpy of formation",
            "method": "PM7 + Delta",
            "dataset": "PM7",
            "model": "Model A",
            "status": "Ready",
        }
    ).render()

    assert panel.title == "Session"


def test_format_session_header_respects_terminal_widths() -> None:
    long_value = "Standard enthalpy of formation with a very long method label"

    for width in (80, 100, 120):
        header = format_session_header(
            property_label=long_value,
            method_label=long_value,
            dataset_label=long_value,
            model_label=long_value,
            status="Ready",
            width=width,
        )
        assert len(header) <= width + 20
        assert "Property:" in header
        assert "Method:" in header
