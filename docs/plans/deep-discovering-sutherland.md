# Phase E: Model Serialization + CLI Integration

## Context

The ML pipeline (Phase D) is complete — `trainer.py` can train a `DeltaLearner` and return metrics. But the trained model lives only in memory. Phase E bridges this gap by:
1. **Serializing** the trained model + pipeline + metrics to disk (`joblib`)
2. **Integrating** real model data into the CLI's Models view (replacing mock data)

This enables the user to train once, save, and reuse the model across sessions.

## Scope

### Files to CREATE
- `src/grimperium/ml/persistence.py` — save/load/metadata functions
- `tests/ml/test_persistence.py` — 7 TDD tests

### Files to MODIFY
- `tests/ml/conftest.py` — add `trained_model_fixture`
- `src/grimperium/cli/views/models_view.py` — replace mock data with real model

### Files NOT to touch
- `trainer.py`, `evaluator.py`, `data_loader.py`, `features.py`
- `mock_data.py`, `controller.py`, `app.py`

## Implementation Plan

### Step 1: RED Phase — Write Tests

**File: `tests/ml/conftest.py`** (append only)
- Add `trained_model_fixture(synthetic_csv_path)` that calls `train()` with `return_pipeline=True` and returns a bundle dict with keys: `learner`, `pipeline`, `metrics`

**File: `tests/ml/test_persistence.py`**
- `TestSaveModel`: 3 tests (creates file, creates parent dirs, includes metadata keys)
- `TestLoadModel`: 4 tests (returns correct types, model is fitted, raises on missing file, roundtrip predictions identical)

Run and confirm all 7 fail with `ModuleNotFoundError`.

### Step 2: GREEN Phase — Implement persistence.py

**File: `src/grimperium/ml/persistence.py`**

Three functions:

1. `save_model(bundle, path)` — validates required keys, creates parent dirs, adds `version` + `trained_at` metadata, writes via `joblib.dump`
2. `load_model(path)` — checks existence, loads bundle, returns `(DeltaLearner, FeaturePipeline)` tuple
3. `load_model_metadata(path)` — loads bundle but returns only metadata dict (version, trained_at, metrics) without learner/pipeline. Useful for CLI display without full deserialization cost.

Key details:
- `MODEL_VERSION = "1.0.0"`
- Bundle format: `{learner, pipeline, metrics, version, trained_at}`
- `trained_at` stored as ISO 8601 UTC string
- Uses `from grimperium import DictStrAny` (exists in `__init__.py`)

Run tests — all 7 should pass.

### Step 3: Modify models_view.py for Real Data

**File: `src/grimperium/cli/views/models_view.py`**

Changes:
1. Remove `from grimperium.cli.mock_data import MODELS, Model`
2. Add imports: `persistence.load_model_metadata`, `persistence.save_model`, `trainer.train`, `Path`
3. Define `_MODEL_PATH = Path("models/delta_learner_v1.joblib")`
4. Add `_load_model_info(self) -> dict | None` helper method
5. Rewrite `render()` — single row for DeltaLearner v1, showing real metrics from metadata or "Not trained" if file missing
6. Enable "Train New Model" menu option (remove `disabled=True`):
   - In `handle_action("train")`: call `train()`, build bundle, call `save_model()`, show Rich panel with results
7. Enable "Test Model" in detail menu (remove `disabled=True`):
   - In `handle_action("test")`: load metadata and display metrics panel
8. Keep "compare", "retrain", "export" as In Development
9. Remove `selected_model` / `Model` dataclass usage — simplify since we have only one real model
10. Keep `render_model_detail` working with metadata dict instead of `Model` dataclass

Existing references to `Model` from `mock_data.py`:
- `results_view.py` imports `MODELS` from `mock_data` — unaffected since we don't touch `mock_data.py`
- `models_view.py` is the ONLY file importing `Model` class — safe to remove

### Step 4: Verification

1. `pytest tests/ml/test_persistence.py -v --tb=short` → 7 passed
2. `pytest tests/ -q --tb=short` → all passing (was 572+20 skipped)
3. `mypy --strict` on both new/modified files
4. `ruff check` + `black --check`
5. Integration script: train real model, save, load metadata, roundtrip predictions

## Key Reusable Existing Code
- `grimperium.ml.trainer.train()` — `trainer.py:41` (with `return_pipeline=True`)
- `grimperium.ml.data_loader.load_ml_data()` — `data_loader.py`
- `grimperium.core.delta_learning.DeltaLearner` — has `is_fitted` property, `predict(X, y_pm7)`
- `grimperium.ml.features.FeaturePipeline` — has `train()`, `transform()`, serializable via joblib
- `grimperium.cli.styles.COLORS` / `ICONS` — use existing status colors/icons
- `grimperium.cli.views.base_view.BaseView` — `show_success()`, `show_error()`, `show_in_development()`, `wait_for_enter()`
- `grimperium.DictStrAny` — type alias from `__init__.py`
