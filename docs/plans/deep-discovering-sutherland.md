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
- `TestLoadModel`: 4 base tests + 4 robustness tests
  - Base: returns correct types, model is fitted, raises on missing file, roundtrip predictions identical
  - Robustness (use `trained_model_fixture` where applicable):
    - Corrupted joblib file (write random bytes, expect clear error from `load_model`)
    - Invalid bundle structure (missing required key like `learner` or `metrics`, expect `ValueError`)
    - `MODEL_VERSION` mismatch (save with altered version, verify `load_model` / `load_model_metadata` behaviour)
    - Wrong-type entries (`learner` value is not a `DeltaLearner`, expect `TypeError` or `ValueError`)

Run and confirm all 11 fail with `ModuleNotFoundError`.

### Step 2: GREEN Phase — Implement persistence.py

**File: `src/grimperium/ml/persistence.py`**

Three functions (all accept `path: Path | str`):

1. `save_model(bundle, path)` — validates required keys, creates parent dirs, adds `version` + `trained_at` metadata, writes atomically via temp file + rename (see **Atomic Write Procedure** below), then validates the written file
2. `load_model(path)` — checks existence, loads bundle, validates required keys and types, returns `(DeltaLearner, FeaturePipeline)` tuple
3. `load_model_metadata(path)` — loads bundle but returns only metadata dict (version, trained_at, metrics) without learner/pipeline. Useful for CLI display without full deserialization cost.

#### Bundle Schema

```
{
  "learner":    DeltaLearner   (fitted, required)
  "pipeline":   FeaturePipeline (fitted, required)
  "metrics":    {              (required, dict[str, dict[str, float]])
      "train": {"rmse": float, "mae": float, "r2": float, ...},
      "test":  {"rmse": float, "mae": float, "r2": float,
                "mape": float, "max_error": float, "gate_pass": bool, ...}
  }
  "version":    str            (required, e.g. "1.0.0")
  "trained_at": str            (required, ISO 8601 UTC, e.g. "2026-03-15T08:00:00+00:00")
}
```

- **Required keys:** `learner`, `pipeline`, `metrics`, `version`, `trained_at`.
- **Canonical metric names:** `rmse`, `mae`, `r2` (both train and test); `mape`, `max_error`, `gate_pass` (test only). Additional numeric metrics are optional.
- **metrics** values come directly from `trainer.train()` output; if the trainer API changes, the caller must transform results to this schema before calling `save_model`.
- `MODEL_VERSION = "1.0.0"` — default constant.
- `trained_at` stored as ISO 8601 UTC string.
- Uses `from grimperium import DictStrAny` (exists in `__init__.py`).

#### Version Compatibility and Migration Policy

- **Major version mismatch** (e.g. loaded "2.0.0" vs current "1.0.0"): `load_model` and `load_model_metadata` must raise `ValueError` with a message describing the mismatch.
- **Minor / patch differences** (e.g. "1.1.0" vs "1.0.0"): log a warning but proceed normally.
- **Backward compatibility:** bundles saved with `MODEL_VERSION = "1.x.y"` must remain loadable by any "1.z.w" loader.
- **Breaking changes** (new required key, changed type of existing key, removed key) require a major version bump.
- **`trained_at` format:** always ISO 8601 UTC. Loaders should parse with `datetime.fromisoformat`.
- **Migration procedure:** when bumping the major version, add a migration function in `persistence.py` that reads the old bundle, transforms it to the new schema, and re-saves. Update `MODEL_VERSION`, document the change in `CHANGELOG.md`, and adjust `load_model` to call the migrator when it detects the old major version.

#### Atomic Write Procedure (save_model)

1. Write to a temp file in the same directory (`path.with_suffix(".joblib.tmp")`).
2. If the target file already exists, copy it to a backup (`path.with_suffix(".joblib.bak")`).
3. Call `joblib.dump(payload, temp_path)` and flush.
4. Rename temp file to final path (atomic on POSIX).
5. Validate by loading the file back (`joblib.load`); on failure, restore backup and raise.

#### Security Notes

- `joblib.load` deserializes objects and can execute arbitrary code. Only load files from trusted sources.
- Document this warning in the `load_model` and `load_model_metadata` docstrings.
- Future: consider adding a SHA-256 checksum written alongside the bundle (e.g. `.joblib.sha256`) and verifying on load.

#### Validation in load_model / load_model_metadata

- Check all required bundle keys (`learner`, `pipeline`, `metrics`, `version`, `trained_at`); raise `ValueError` on missing keys.
- Check types: `learner` must be `DeltaLearner`, `pipeline` must be `FeaturePipeline`, `metrics` must be `dict`.
- Wrap `joblib.load` errors to surface corrupted or incomplete files with a clear message.

Run tests — all 11 should pass.

### Step 3: Modify models_view.py for Real Data

**File: `src/grimperium/cli/views/models_view.py`**

Changes:
1. Remove `from grimperium.cli.mock_data import MODELS, Model`
2. Add imports: `persistence.load_model_metadata`, `persistence.save_model`, `trainer.train`, `Path`, `os`
3. Define configurable paths via environment variables with sensible defaults:
   - `_MODEL_PATH = Path(os.environ.get("GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib"))`
   - `_DATA_PATH = Path(os.environ.get("GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"))`
   - Future multi-model support: allow the env var to point to a directory; treat filenames like `delta_learner_v1.joblib` or timestamped files within it; add an optional `model_path` parameter to `render_model_detail` and `_load_model_info` to override the global default.
4. Add `_safe_metric(value, default=0.0) -> float` helper to coerce `None` values before formatting.
5. Add `_load_model_info(self) -> dict | None` helper method
6. Rewrite `render()` — single row for DeltaLearner v1, showing real metrics from metadata or "Not trained" if file missing
7. Enable "Train New Model" menu option (remove `disabled=True`):
   - In `handle_action("train")`: call `train()` using `_DATA_PATH`, validate the returned tuple structure (check `isinstance(result, tuple)` and length), build bundle, call `save_model()`, show Rich panel with results; wrap unpacking in `try/except (ValueError, TypeError)` with a clear error message referencing `train_model` and `return_pipeline`.
8. Enable "Test Model" in detail menu (remove `disabled=True`):
   - In `handle_action("test")`: load metadata and display metrics panel; use `_safe_metric` for all `.get()` values before format specifiers
9. Keep "compare", "retrain", "export" as In Development
10. Remove `selected_model` / `Model` dataclass usage — simplify since we have only one real model
11. Keep `render_model_detail` working with metadata dict instead of `Model` dataclass; use `_safe_metric` for all metric formatting

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
