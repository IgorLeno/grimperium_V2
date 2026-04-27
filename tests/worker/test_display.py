"""Tests for worker/display.py — MoleculeProgressDisplay."""

from io import StringIO

from rich.console import Console

from grimperium.worker.display import (
    MoleculeProgressDisplay,
    _format_elapsed,
    _status_style,
)


def _headless_display() -> MoleculeProgressDisplay:
    """Display with a non-TTY console so Live renders synchronously in tests."""
    console = Console(file=StringIO(), force_terminal=False, highlight=False)
    return MoleculeProgressDisplay(console=console)


class TestFormatElapsed:
    def test_seconds_only(self) -> None:
        assert _format_elapsed(45) == "00:45"

    def test_minutes_and_seconds(self) -> None:
        assert _format_elapsed(125) == "02:05"

    def test_hours(self) -> None:
        assert _format_elapsed(3661) == "01:01:01"

    def test_zero(self) -> None:
        assert _format_elapsed(0) == "00:00"


class TestStatusStyle:
    def test_done_is_green(self) -> None:
        assert _status_style("DONE") == "green"

    def test_ok_is_green(self) -> None:
        assert _status_style("ok") == "green"

    def test_failed_is_red(self) -> None:
        assert _status_style("FAILED") == "red"

    def test_crest_is_yellow(self) -> None:
        assert _status_style("CREST") == "yellow"

    def test_mopac_is_yellow(self) -> None:
        assert _status_style("mopac") == "yellow"

    def test_idle_is_dim(self) -> None:
        assert _status_style("IDLE") == "dim"

    def test_unknown_is_cyan(self) -> None:
        assert _status_style("COMPUTING") == "cyan"


class TestMoleculeProgressDisplay:
    def test_update_sets_state(self) -> None:
        display = _headless_display()
        display.update("mol_001", "CCO", "CREST")
        assert display._mol_id == "mol_001"
        assert display._smiles == "CCO"
        assert display._status == "CREST"

    def test_update_with_result(self) -> None:
        display = _headless_display()
        display.update("mol_001", "CCO", "DONE", result="OK")
        assert display._result == "OK"

    def test_initial_state(self) -> None:
        display = _headless_display()
        assert display._mol_id == ""
        assert display._status == "IDLE"

    def test_render_snapshot_contains_mol_id(self) -> None:
        display = _headless_display()
        display.update("mol_042", "CCO", "CREST")
        snapshot = display.render_snapshot()
        assert "mol_042" in snapshot

    def test_render_snapshot_contains_status(self) -> None:
        display = _headless_display()
        display.update("mol_001", "CCO", "MOPAC")
        snapshot = display.render_snapshot()
        assert "MOPAC" in snapshot

    def test_render_snapshot_truncates_long_smiles(self) -> None:
        long_smiles = "C" * 50
        display = _headless_display()
        display.update("mol_001", long_smiles, "CREST")
        snapshot = display.render_snapshot()
        assert "..." in snapshot
        assert long_smiles not in snapshot

    def test_render_snapshot_keeps_short_smiles(self) -> None:
        display = _headless_display()
        display.update("mol_001", "CCO", "CREST")
        snapshot = display.render_snapshot()
        assert "CCO" in snapshot

    def test_render_snapshot_contains_elapsed(self) -> None:
        display = _headless_display()
        display.update("mol_001", "CCO", "CREST")
        snapshot = display.render_snapshot()
        # Should contain at least MM:SS format
        import re

        assert re.search(r"\d{2}:\d{2}", snapshot)

    def test_start_stop_does_not_crash(self) -> None:
        display = _headless_display()
        display.start()
        display.update("mol_001", "CCO", "CREST")
        display.stop()

    def test_stop_without_start_is_safe(self) -> None:
        display = _headless_display()
        display.stop()

    def test_start_twice_is_safe(self) -> None:
        display = _headless_display()
        display.start()
        display.start()  # Should not raise
        display.stop()

    def test_new_mol_resets_timer(self) -> None:
        import time

        display = _headless_display()
        display.update("mol_001", "CCO", "CREST")
        t1 = display._start_time
        time.sleep(0.01)
        display.update("mol_002", "CCCO", "CREST")
        t2 = display._start_time
        assert t2 > t1

    def test_same_mol_keeps_timer(self) -> None:
        display = _headless_display()
        display.update("mol_001", "CCO", "CREST")
        t1 = display._start_time
        display.update("mol_001", "CCO", "MOPAC")
        assert display._start_time == t1
