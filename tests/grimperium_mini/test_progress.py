"""Tests for MiniProgressTracker."""

from grimperium_mini.progress import MiniProgressTracker


def _make_row(mol_id: str, status: str = "success") -> dict[str, object]:
    return {
        "mol_id": mol_id,
        "method": "PM7",
        "status": status,
        "hof_selected_kJmol": -100.5,
        "total_time_s": 1.23,
        "crest_status": "success",
        "mopac_status": "success",
    }


def test_on_task_done_registers_row() -> None:
    tracker = MiniProgressTracker(total_tasks=3)
    row = _make_row("mol_001")
    tracker.on_task_done(row)
    assert tracker.rows[-1]["mol_id"] == row["mol_id"]


def test_summary_counts_success_and_failed() -> None:
    tracker = MiniProgressTracker(total_tasks=3)
    tracker.on_task_done(_make_row("m1", "success"))
    tracker.on_task_done(_make_row("m2", "success"))
    tracker.on_task_done(_make_row("m3", "failed"))
    panel = tracker.summary()
    rendered = str(panel.renderable)
    assert "2" in rendered
    assert "1" in rendered


def test_render_returns_table() -> None:
    from rich.table import Table

    tracker = MiniProgressTracker(total_tasks=2)
    tracker.on_task_done(_make_row("mol_a"))
    table = tracker.render()
    assert isinstance(table, Table)
    assert table.row_count == 1


def test_on_task_start_does_not_raise() -> None:
    tracker = MiniProgressTracker(total_tasks=1)
    tracker.on_task_start("mol_001", "AM1")
