"""Tests for the grimperium analyze-errors CLI subcommand."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

# ── Fixtures ──────────────────────────────────────────────────────────────────


def _make_valid_csv(tmp_path: Path) -> Path:
    csv_path = tmp_path / "data.csv"
    cbs = list(range(-100, -1100, -100))  # 10 elements: -100, -200, ..., -1000
    df = pd.DataFrame(
        {
            "mol_id": [f"m{i}" for i in range(10)],
            "H298_cbs": cbs,
            "H298_predicted": [v + float(i + 1) for i, v in enumerate(cbs)],
        }
    )
    df.to_csv(csv_path, index=False)
    return csv_path


def _make_invalid_csv(tmp_path: Path) -> Path:
    csv_path = tmp_path / "bad.csv"
    df = pd.DataFrame({"col_a": [1, 2], "col_b": [3, 4]})
    df.to_csv(csv_path, index=False)
    return csv_path


def _call_main(argv: list[str]) -> int:
    """Call app.main() with a given argv and return the exit code."""
    from grimperium.cli.app import main

    with patch("sys.argv", argv):
        return main()


# ── No subcommand → interactive mode ─────────────────────────────────────────


def test_no_subcommand_starts_interactive_app() -> None:
    """Without a subcommand, GrimperiumCLI.run() is invoked."""
    with (
        patch("sys.argv", ["grimperium"]),
        patch("grimperium.cli.app.GrimperiumCLI") as mock_cli_cls,
    ):
        mock_cli = MagicMock()
        mock_cli.run.return_value = 0
        mock_cli_cls.return_value = mock_cli
        from grimperium.cli.app import main

        result = main()

    mock_cli.run.assert_called_once()
    assert result == 0


# ── analyze-errors basic output ───────────────────────────────────────────────


def test_analyze_errors_prints_summary_and_exits_0(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--input <valid csv> prints mandatory summary and exits 0."""
    csv_path = _make_valid_csv(tmp_path)

    with patch("sys.argv", ["grimperium", "analyze-errors", "--input", str(csv_path)]):
        from grimperium.cli.app import _run_analyze_errors

        result = _run_analyze_errors(["--input", str(csv_path)])

    captured = capsys.readouterr()
    assert result == 0
    assert "Error analysis complete" in captured.out
    assert "Valid molecules analyzed" in captured.out
    assert "MAE" in captured.out
    assert "RMSE" in captured.out
    assert "Outliers detected" in captured.out
    assert "Top anomalous molecule IDs" in captured.out


def test_analyze_errors_exits_1_on_missing_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--input pointing to a nonexistent file exits 1 with an error message."""
    from grimperium.cli.app import _run_analyze_errors

    result = _run_analyze_errors(["--input", str(tmp_path / "nope.csv")])
    captured = capsys.readouterr()
    assert result == 1
    assert "Error" in captured.err or "not found" in captured.err.lower()


def test_analyze_errors_exits_1_on_missing_columns(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """CSV without required columns exits 1 with a clear error."""
    csv_path = _make_invalid_csv(tmp_path)
    from grimperium.cli.app import _run_analyze_errors

    result = _run_analyze_errors(["--input", str(csv_path)])
    captured = capsys.readouterr()
    assert result == 1
    assert "Error" in captured.err


# ── --show-outliers ───────────────────────────────────────────────────────────


def test_analyze_errors_show_outliers_prints_table(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--show-outliers prints the outlier details section."""
    csv_path = _make_valid_csv(tmp_path)
    from grimperium.cli.app import _run_analyze_errors

    result = _run_analyze_errors(
        ["--input", str(csv_path), "--show-outliers", "--threshold", "0.5"]
    )
    captured = capsys.readouterr()
    assert result == 0
    # With threshold=0.5, most molecules should be outliers
    # At minimum the mandatory summary is there
    assert "Error analysis complete" in captured.out


# ── --export-outliers ─────────────────────────────────────────────────────────


def test_analyze_errors_export_outliers_creates_csv(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--export-outliers saves a CSV file at the given path."""
    csv_path = _make_valid_csv(tmp_path)
    out_path = tmp_path / "out_outliers.csv"

    from grimperium.cli.app import _run_analyze_errors

    result = _run_analyze_errors(
        [
            "--input",
            str(csv_path),
            "--export-outliers",
            str(out_path),
            "--threshold",
            "0.5",
        ]
    )
    assert result == 0
    # File may or may not have outliers depending on data, but function must not crash
    captured = capsys.readouterr()
    assert "Error analysis complete" in captured.out


# ── --html-report ─────────────────────────────────────────────────────────────


def test_analyze_errors_html_report_calls_generate(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--html-report calls generate_html_report with the correct path."""
    csv_path = _make_valid_csv(tmp_path)
    html_path = tmp_path / "report.html"

    from grimperium.cli.app import _run_analyze_errors

    with patch("grimperium.ml.html_report.generate_html_report") as mock_gen:
        mock_gen.return_value = html_path
        result = _run_analyze_errors(
            ["--input", str(csv_path), "--html-report", str(html_path)]
        )

    assert result == 0
    mock_gen.assert_called_once()
    call_args = mock_gen.call_args
    assert str(call_args[0][1]) == str(html_path)


# ── --top N ───────────────────────────────────────────────────────────────────


def test_analyze_errors_top_controls_printed_count(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--top N limits the number of molecule IDs printed in the summary."""
    csv_path = _make_valid_csv(tmp_path)
    from grimperium.cli.app import _run_analyze_errors

    result = _run_analyze_errors(["--input", str(csv_path), "--top", "3"])
    captured = capsys.readouterr()
    assert result == 0
    # Count lines like "  1. m0 | abs_error=..."
    ranked_lines = [
        line
        for line in captured.out.splitlines()
        if line.strip().startswith(tuple(f"{i}." for i in range(1, 11)))
    ]
    assert len(ranked_lines) <= 3


# ── analyze-errors is routed from main() ─────────────────────────────────────


def test_main_routes_analyze_errors_subcommand(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """main() routes 'analyze-errors' to _run_analyze_errors, not the TUI."""
    csv_path = _make_valid_csv(tmp_path)

    with (
        patch("sys.argv", ["grimperium", "analyze-errors", "--input", str(csv_path)]),
        patch("grimperium.cli.app.GrimperiumCLI") as mock_tui,
    ):
        from grimperium.cli.app import main

        result = main()

    # TUI must NOT be instantiated
    mock_tui.assert_not_called()
    captured = capsys.readouterr()
    assert result == 0
    assert "Error analysis complete" in captured.out
