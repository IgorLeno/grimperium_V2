"""Tests for DatabasesView catalog management actions."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock

from rich.console import Console

from grimperium.cli.controller import CliController
from grimperium.cli.views.databases_view import DatabasesView


def _view(tmp_path: Path, monkeypatch) -> DatabasesView:
    monkeypatch.setenv("GRIMPERIUM_CONFIG_DIR", str(tmp_path / "config"))
    data_dir = tmp_path / "data"
    data_dir.mkdir(exist_ok=True)
    controller = CliController()
    controller.console = Console(file=io.StringIO(), force_terminal=True, width=120)
    monkeypatch.setattr("grimperium.cli.views.databases_view.DATA_DIR", data_dir)
    view = DatabasesView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view


def test_add_database_wizard_persists_user_overlay(
    tmp_path: Path,
    monkeypatch,
) -> None:
    csv_path = tmp_path / "custom.csv"
    csv_path.write_text("smiles,H298_cbs\nC,-10\n", encoding="utf-8")
    view = _view(tmp_path, monkeypatch)
    answers = iter(
        [
            str(csv_path),
            "Custom DB",
            "CUSTOM",
            "Custom description",
            "readable,analysis_input",
        ]
    )

    monkeypatch.setattr(
        "grimperium.cli.views.databases_view.text_input",
        lambda *args, **kwargs: next(answers),
    )
    monkeypatch.setattr(
        "grimperium.cli.views.databases_view.show_menu",
        lambda *args, **kwargs: "analysis",
    )

    view._handle_add_database_wizard()

    added = view.registry.get_by_alias("CUSTOM")
    assert added is not None
    assert added.origin == "user"
    assert added.database_id.startswith("user.")
    assert added.capabilities == frozenset({"readable", "analysis_input"})


def test_use_as_session_dataset_updates_controller_session(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "thermo_pm7.csv").write_text("smiles,status\nC,OK\n", encoding="utf-8")
    view = _view(tmp_path, monkeypatch)
    view.selected_db = view.registry.get_by_id("official.crest_pm7")

    view._handle_use_as_session_dataset()

    assert view.controller.session.dataset is not None
    assert view.controller.session.dataset.database_id == "official.crest_pm7"
    assert view.controller.current_csv_path == data_dir / "thermo_pm7.csv"


def test_detail_actions_validate_capabilities(
    tmp_path: Path,
    monkeypatch,
) -> None:
    view = _view(tmp_path, monkeypatch)
    view.selected_db = view.registry.get_by_id("official.nist_experimental")
    view.show_error = MagicMock()  # type: ignore[method-assign]

    assert view.handle_action("analyze") is None

    view.show_error.assert_called_once()
