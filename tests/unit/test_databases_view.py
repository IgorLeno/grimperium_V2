"""Tests for databases view module."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

from grimperium.cli.views.databases_view import DatabasesView


def test_load_real_phase_a_results_from_json(tmp_path: Path) -> None:
    """
    Bug #2: Test loading molecule count from phase_a_results.json.

    When JSON exists with n_molecules, that value should be used.
    """
    json_file = tmp_path / "phase_a_results.json"
    json_data = {
        "n_molecules": 29568,
        "results": [
            {"smiles": "CCO", "H298_pm7": -100.5, "timestamp": "2026-01-15T10:30:00"}
        ],
    }
    json_file.write_text(json.dumps(json_data), encoding="utf-8")

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        result = DatabasesView.load_real_phase_a_results()

    assert result is not None
    assert result.molecules == 29568
    assert result.name == "CREST PM7"
    assert result.status == "ready"  # Because molecules > 0


def test_load_real_phase_a_results_from_csv_fallback(tmp_path: Path) -> None:
    """
    Bug #2: Test fallback to CSV row count when JSON doesn't exist.

    When phase_a_results.json is missing, should count rows in thermo_pm7.csv.
    """
    # JSON doesn't exist
    json_file = tmp_path / "phase_a_results.json"

    # CSV exists with data
    csv_file = tmp_path / "thermo_pm7.csv"
    csv_content = "smiles,H298_pm7,charge,multiplicity\nCCO,-100.5,0,1\nCC,-50.2,0,1\nC,-25.1,0,1\n"
    csv_file.write_text(csv_content, encoding="utf-8")

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
            result = DatabasesView.load_real_phase_a_results()

    assert result is not None
    assert result.molecules == 3  # 3 data rows (excluding header)
    assert result.name == "CREST PM7"
    assert result.status == "ready"


def test_load_real_phase_a_results_empty_csv_fallback(tmp_path: Path) -> None:
    """
    Bug #2: Test CSV fallback with empty CSV (only header).

    Should return Database with 0 molecules and in_development status.
    """
    json_file = tmp_path / "phase_a_results.json"
    csv_file = tmp_path / "thermo_pm7.csv"
    csv_file.write_text("smiles,H298_pm7,charge,multiplicity\n", encoding="utf-8")

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
            result = DatabasesView.load_real_phase_a_results()

    assert result is not None
    assert result.molecules == 0
    assert result.status == "in_development"


def test_load_real_phase_a_results_no_files(tmp_path: Path) -> None:
    """
    Bug #2: Test when neither JSON nor CSV exists.

    Should return None to use mock data fallback.
    """
    json_file = tmp_path / "phase_a_results.json"
    # Neither file exists

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
            result = DatabasesView.load_real_phase_a_results()

    assert result is None


def test_load_real_phase_a_results_invalid_json(tmp_path: Path) -> None:
    """
    Bug #2: Test handling of corrupted JSON file.

    Should fallback to CSV count when JSON is invalid.
    """
    json_file = tmp_path / "phase_a_results.json"
    json_file.write_text("{ invalid json }", encoding="utf-8")

    csv_file = tmp_path / "thermo_pm7.csv"
    csv_content = "smiles,H298_pm7\nCCO,-100.5\n"
    csv_file.write_text(csv_content, encoding="utf-8")

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
            result = DatabasesView.load_real_phase_a_results()

    assert result is not None
    assert result.molecules == 1  # Fallback to CSV count


def test_load_real_phase_a_results_json_with_zero_molecules(tmp_path: Path) -> None:
    """
    Bug #2: Test JSON with n_molecules=0 falls back to CSV.

    Even if JSON says 0 molecules, check CSV for actual data.
    """
    json_file = tmp_path / "phase_a_results.json"
    json_data = {"n_molecules": 0, "results": []}
    json_file.write_text(json.dumps(json_data), encoding="utf-8")

    csv_file = tmp_path / "thermo_pm7.csv"
    csv_content = "smiles,H298_pm7\nCCO,-100.5\nCC,-50.2\n"
    csv_file.write_text(csv_content, encoding="utf-8")

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
            result = DatabasesView.load_real_phase_a_results()

    assert result is not None
    assert result.molecules == 2  # Use CSV count, not JSON's 0


def test_refresh_databases_from_filesystem_finds_csvs(tmp_path: Path) -> None:
    """
    Bug #4: Test that refresh discovers CSV files in data/ directory.

    Should scan filesystem and return count of discovered databases.
    """
    # Create mock controller
    controller = MagicMock()
    view = DatabasesView(controller)

    # Create test CSV files
    (tmp_path / "thermo_cbs_chon.csv").write_text(
        "smiles,H298_cbs\nCCO,-100.5\nCC,-50.2\n", encoding="utf-8"
    )
    (tmp_path / "thermo_pm7.csv").write_text(
        "smiles,H298_pm7\nCCO,-100.3\n", encoding="utf-8"
    )

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        count = view.refresh_databases_from_filesystem()

    assert count == 2, "Should discover 2 CSV files"


def test_refresh_databases_from_filesystem_missing_directory(tmp_path: Path) -> None:
    """
    Bug #4: Test refresh handles missing data/ directory gracefully.

    Should return 0 and not crash when data/ doesn't exist.
    """
    controller = MagicMock()
    view = DatabasesView(controller)

    nonexistent_dir = tmp_path / "does_not_exist"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", nonexistent_dir):
        count = view.refresh_databases_from_filesystem()

    assert count == 0, "Should return 0 when directory doesn't exist"


def test_refresh_databases_from_filesystem_empty_directory(tmp_path: Path) -> None:
    """
    Bug #4: Test refresh handles empty data/ directory.

    Should return 0 when no CSV files exist.
    """
    controller = MagicMock()
    view = DatabasesView(controller)

    # Create empty directory
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    with patch("grimperium.cli.views.databases_view.DATA_DIR", empty_dir):
        count = view.refresh_databases_from_filesystem()

    assert count == 0, "Should return 0 when directory is empty"


def test_refresh_databases_from_filesystem_counts_molecules(tmp_path: Path) -> None:
    """
    Bug #4: Test that refresh counts actual CSV rows for each database.

    Should display molecule count for each discovered database.
    """
    controller = MagicMock()
    view = DatabasesView(controller)

    # Create CSV with known row count
    (tmp_path / "test_database.csv").write_text(
        "smiles,energy\nCCO,-100\nCC,-50\nC,-25\n", encoding="utf-8"
    )

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        count = view.refresh_databases_from_filesystem()

    assert count == 1, "Should discover 1 CSV file"
