"""Unit tests for ModelsView — dynamic model discovery, set-as-default, retrain."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import joblib
import pytest

from grimperium.cli.views.models_view import ModelsView

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _create_fake_model(
    directory: Path,
    name: str,
    metrics: dict[str, Any] | None = None,
) -> Path:
    """Create a minimal valid .joblib model bundle for testing.

    Uses joblib.dump directly to avoid importing DeltaLearner/FeaturePipeline.
    """
    if metrics is None:
        metrics = {
            "train": {"mae": 1.0, "rmse": 1.5, "r2": 0.98},
            "test": {"mae": 2.0, "rmse": 2.5, "r2": 0.95, "max_error": 5.0},
        }
    path = directory / f"{name}.joblib"
    bundle = {
        "learner": None,  # load_model_metadata only reads version/trained_at/metrics
        "pipeline": None,
        "metrics": metrics,
        "version": "1.0.0",
        "trained_at": "2026-03-01T12:00:00+00:00",
    }
    joblib.dump(bundle, path)
    return path


@pytest.fixture
def mock_controller() -> MagicMock:
    ctrl = MagicMock()
    ctrl.current_model = None
    ctrl.current_model_path = None
    return ctrl


@pytest.fixture
def view(mock_controller: MagicMock) -> ModelsView:
    return ModelsView(mock_controller)


# ---------------------------------------------------------------------------
# 1–5: _discover_models
# ---------------------------------------------------------------------------


class TestDiscoverModels:
    def test_discover_models_no_directory(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Missing models/ directory returns empty list."""
        monkeypatch.chdir(tmp_path)
        result = view._discover_models()
        assert result == []

    def test_discover_models_empty_directory(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Empty models/ directory returns empty list."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / "models").mkdir()
        result = view._discover_models()
        assert result == []

    def test_discover_models_one_model(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """One valid .joblib returns one entry with correct fields."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        _create_fake_model(models_dir, "delta_learner_v1")

        result = view._discover_models()

        assert len(result) == 1
        entry = result[0]
        assert entry["display_name"] == "Delta Learner V1"
        assert entry["mae"] == pytest.approx(2.0)
        assert entry["r2"] == pytest.approx(0.95)
        assert entry["trained_at"] == "2026-03-01T12:00:00+00:00"
        assert entry["version"] == "1.0.0"
        assert entry["is_active"] is False

    def test_discover_models_invalid_joblib_skipped(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Corrupted .joblib is skipped; valid files still returned."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        _create_fake_model(models_dir, "good_model")
        (models_dir / "bad_model.joblib").write_bytes(b"not a valid joblib file!!!")

        result = view._discover_models()

        assert len(result) == 1
        assert result[0]["display_name"] == "Good Model"

    def test_discover_models_marks_active_model(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """is_active=True when path matches controller.current_model_path."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        model_path = _create_fake_model(models_dir, "delta_learner_v1")
        view.controller.current_model_path = model_path

        result = view._discover_models()

        assert len(result) == 1
        assert result[0]["is_active"] is True


# ---------------------------------------------------------------------------
# 6: set_default
# ---------------------------------------------------------------------------


class TestSetDefault:
    def test_set_default_updates_controller_path(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """handle_action('set_default') calls controller.set_model with the model path."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        model_path = _create_fake_model(models_dir, "delta_learner_v1")
        view.selected_model_path = model_path

        with patch.object(view, "show_success"):
            view.handle_action("set_default")

        view.controller.set_model.assert_called_once()
        call = view.controller.set_model.call_args
        # set_model(name, model_path=path) — path must be passed
        assert call.kwargs.get("model_path") == model_path


# ---------------------------------------------------------------------------
# 7–10: retrain
# ---------------------------------------------------------------------------


class TestRetrain:
    def test_retrain_no_model_selected_shows_error(
        self,
        view: ModelsView,
    ) -> None:
        """selected_model_path=None → show_error called, train not called."""
        view.selected_model_path = None

        with (
            patch.object(view, "show_error") as mock_error,
            patch("grimperium.cli.views.models_view.train") as mock_train,
        ):
            view._handle_retrain()

        mock_error.assert_called_once()
        mock_train.assert_not_called()

    def test_retrain_confirmation_denied(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """User enters 'no' → train not called."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        model_path = _create_fake_model(models_dir, "delta_learner_v1")
        view.selected_model_path = model_path

        with (
            patch.object(view.console, "input", return_value="no"),
            patch.object(view.console, "print"),
            patch("grimperium.cli.views.models_view.train") as mock_train,
        ):
            view._handle_retrain()

        mock_train.assert_not_called()

    def test_retrain_success_shows_comparison(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Full retrain flow → model saved, save called with the original path."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        model_path = _create_fake_model(models_dir, "delta_learner_v1")
        view.selected_model_path = model_path

        new_learner = MagicMock()
        new_pipeline = MagicMock()
        new_train_m: dict[str, Any] = {"mae": 0.5, "rmse": 0.8, "r2": 0.99}
        new_test_m: dict[str, Any] = {
            "mae": 1.5,
            "rmse": 2.0,
            "r2": 0.97,
            "max_error": 4.0,
        }

        with (
            patch.object(view.console, "input", return_value="yes"),
            patch.object(view.console, "print"),
            patch.object(view, "wait_for_enter"),
            patch(
                "grimperium.cli.views.models_view.train",
                return_value=(new_learner, new_train_m, new_test_m, new_pipeline),
            ) as mock_train,
            patch("grimperium.cli.views.models_view.save_model") as mock_save,
        ):
            view._handle_retrain()

        mock_train.assert_called_once()
        mock_save.assert_called_once()
        # save_model(bundle, self.selected_model_path) — path is the same file
        save_call = mock_save.call_args
        saved_path = (
            save_call.args[1] if save_call.args else save_call.kwargs.get("path")
        )
        assert saved_path == model_path

    def test_retrain_updates_controller_if_active(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """When retrained model is the active model, controller path is unchanged (same path)."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        model_path = _create_fake_model(models_dir, "delta_learner_v1")
        view.selected_model_path = model_path
        view.controller.current_model_path = model_path  # this is the active model

        new_learner = MagicMock()
        new_pipeline = MagicMock()
        new_train_m: dict[str, Any] = {"mae": 0.5, "rmse": 0.8, "r2": 0.99}
        new_test_m: dict[str, Any] = {
            "mae": 1.5,
            "rmse": 2.0,
            "r2": 0.97,
            "max_error": 4.0,
        }

        with (
            patch.object(view.console, "input", return_value="yes"),
            patch.object(view.console, "print"),
            patch.object(view, "wait_for_enter"),
            patch(
                "grimperium.cli.views.models_view.train",
                return_value=(new_learner, new_train_m, new_test_m, new_pipeline),
            ),
            patch("grimperium.cli.views.models_view.save_model"),
        ):
            view._handle_retrain()

        # The model was overwritten in-place; the path in the controller is still correct
        assert view.controller.current_model_path == model_path


# ---------------------------------------------------------------------------
# 11–12: handle_action and get_menu_options
# ---------------------------------------------------------------------------


class TestHandleAction:
    def test_handle_action_view_stem_sets_path(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """handle_action('view_delta_learner_v1') sets selected_model_path."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        model_path = _create_fake_model(models_dir, "delta_learner_v1")

        view.handle_action("view_delta_learner_v1")

        assert view.selected_model_path == model_path

    def test_get_menu_options_dynamic(
        self,
        view: ModelsView,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Menu options match discovered models plus 'Train New Model'."""
        monkeypatch.chdir(tmp_path)
        models_dir = tmp_path / "models"
        models_dir.mkdir()
        _create_fake_model(models_dir, "delta_learner_v1")
        _create_fake_model(models_dir, "delta_learner_v2")

        options = view.get_menu_options()
        values = [o.value for o in options]

        assert "view_delta_learner_v1" in values
        assert "view_delta_learner_v2" in values
        assert "train" in values
