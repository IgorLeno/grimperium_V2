"""Tests for split metadata persisted with trained model bundles."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.ml.persistence import load_model_metadata, save_model
from grimperium.ml.trainer import train


def _write_fifty_row_fixture(synthetic_df: pd.DataFrame, tmp_path: Path) -> Path:
    """Create a 50-row OK/OK CSV compatible with the ML loader."""
    eligible = synthetic_df.loc[
        (synthetic_df["status"] == "OK") & (synthetic_df["cbs_quality_flag"] == "OK")
    ].copy()

    replicas: list[pd.DataFrame] = []
    for batch_idx in range(5):
        chunk = eligible.copy()
        chunk["mol_id"] = [
            f"mol_{batch_idx:02d}_{row_idx:03d}" for row_idx in range(len(chunk.index))
        ]
        chunk["H298_cbs"] = chunk["H298_cbs"] + (batch_idx * 0.1)
        chunk["H298_pm7"] = chunk["H298_pm7"] + (batch_idx * 0.1)
        replicas.append(chunk)

    fixture_df = pd.concat(replicas, ignore_index=True)
    csv_path = tmp_path / "thermo_pm7_split_fixture.csv"
    fixture_df.to_csv(csv_path, index=False)
    return csv_path


def test_train_metrics_contain_n_samples(
    synthetic_df: pd.DataFrame,
    tmp_path: Path,
) -> None:
    """train() exposes split sizes in both metric dictionaries."""
    csv_path = _write_fifty_row_fixture(synthetic_df, tmp_path)

    _, train_metrics, test_metrics = train(csv_path, test_size=0.2, random_state=42)

    assert train_metrics["n_samples"] == 40
    assert test_metrics["n_samples"] == 10
    assert train_metrics["n_total"] == 50
    assert test_metrics["n_total"] == 50


def test_load_model_metadata_with_n_samples(
    synthetic_df: pd.DataFrame,
    tmp_path: Path,
) -> None:
    """Saved bundles expose split metadata through load_model_metadata()."""
    csv_path = _write_fifty_row_fixture(synthetic_df, tmp_path)
    learner, train_metrics, test_metrics, pipeline = train(
        csv_path,
        test_size=0.2,
        random_state=42,
        return_pipeline=True,
    )
    bundle: dict[str, Any] = {
        "learner": learner,
        "pipeline": pipeline,
        "metrics": {"train": train_metrics, "test": test_metrics},
    }
    model_path = tmp_path / "split_metadata_model.joblib"

    save_model(bundle, model_path)
    meta = load_model_metadata(model_path)

    assert meta["n_train"] is not None
    assert meta["n_test"] is not None
    assert meta["n_total"] is not None
    assert meta["test_size"] is not None
    assert isinstance(meta["n_train"], int)
    assert isinstance(meta["n_test"], int)
    assert isinstance(meta["n_total"], int)
    assert isinstance(meta["test_size"], float)


def test_load_model_metadata_old_bundle_graceful(
    trained_model_fixture: dict[str, Any],
    tmp_path: Path,
) -> None:
    """Bundles without split metadata still load without raising."""
    legacy_bundle = {
        "learner": trained_model_fixture["learner"],
        "pipeline": trained_model_fixture["pipeline"],
        "metrics": {
            "train": {
                key: value
                for key, value in trained_model_fixture["metrics"]["train"].items()
                if key not in {"n_samples", "n_total", "test_size"}
            },
            "test": {
                key: value
                for key, value in trained_model_fixture["metrics"]["test"].items()
                if key not in {"n_samples", "n_total", "test_size"}
            },
        },
    }
    model_path = tmp_path / "legacy_model.joblib"

    save_model(legacy_bundle, model_path)
    meta = load_model_metadata(model_path)

    assert meta["n_train"] is None
    assert meta["n_test"] is None
    assert meta["n_total"] is None
    assert meta["test_size"] is None
