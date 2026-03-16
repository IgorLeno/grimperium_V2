# Phase F: Batch Predict — Implementation Plan

## Context

Phase F adds batch prediction capability to Grimperium. Given a trained DeltaLearner model (`models/delta_learner_v1.joblib`), we apply it to all eligible molecules in `data/thermo_pm7.csv` and write `H298_predicted` + `delta_correction` columns back to the CSV. The results_view.py is updated to show real model data instead of mocks.

**TDD approach:** tests first (RED), implementation second (GREEN).

---

## Files to Create/Modify

| Action | File | Purpose |
|--------|------|---------|
| CREATE | `tests/ml/test_predictor.py` | 8 tests for batch prediction |
| CREATE | `src/grimperium/ml/predictor.py` | Batch prediction logic |
| MODIFY | `tests/ml/conftest.py` | Add `csv_with_predictions` fixture (append only) |
| MODIFY | `src/grimperium/cli/views/models_view.py` | Fix env var bug: module constants → instance methods |
| MODIFY | `src/grimperium/cli/views/results_view.py` | Replace mock data with real model/CSV data |

---

## Step 1: RED Phase — Write Tests

### 1a. Append fixture to `tests/ml/conftest.py`

Add `csv_with_predictions` fixture at end of file. It copies `synthetic_csv_path` (12-row CSV: 10 OK/OK + 2 non-eligible) to tmp_path for write testing.

### 1b. Create `tests/ml/test_predictor.py`

8 tests in `TestPredictBatch`:
1. `test_returns_dataframe_with_prediction_columns` — result has `H298_predicted`, `delta_correction`
2. `test_predictions_only_for_eligible_rows` — only status==OK & cbs_quality_flag==OK get predictions; others NaN
3. `test_delta_correction_equals_predicted_minus_pm7` — math check
4. `test_csv_is_written_with_new_columns` — CSV on disk updated
5. `test_original_columns_preserved` — existing columns untouched
6. `test_raises_on_missing_model` — `FileNotFoundError("Model file not found")`
7. `test_raises_on_missing_csv` — `FileNotFoundError("CSV file not found")`
8. `test_returns_summary_stats` — `return_stats=True` gives dict with n_predicted, mean_delta, std_delta

**Verify RED:** `pytest tests/ml/test_predictor.py -v --tb=short` → 8 failures (ImportError)

---

## Step 2: GREEN Phase — Implement `predictor.py`

### File: `src/grimperium/ml/predictor.py`

**Key function:** `predict_batch(csv_path, model_path, *, return_stats=False)`

**Algorithm:**
1. Validate paths exist (raise `FileNotFoundError` with specific messages)
2. `load_model(model_path)` → `(learner, pipeline)`
3. `pd.read_csv(csv_path)` → full DataFrame
4. Init `H298_predicted` and `delta_correction` columns as NaN
5. Build eligible mask: `status=="OK" & cbs_quality_flag=="OK"`
6. Extract eligible subset, run `pipeline.transform()` → X
7. Extract `y_pm7` from eligible rows
8. `learner.predict(X, y_pm7)` → `H298_predicted` array
9. Compute `delta_correction = H298_predicted - y_pm7`
10. Assign back to df using `.loc[eligible_mask]`
11. Overwrite CSV: `df.to_csv(csv_path, index=False)`
12. If `return_stats=True`: return `(df, stats_dict)` with n_predicted, mean_delta, std_delta, min_predicted, max_predicted

**Reused components:**
- `grimperium.ml.persistence.load_model` — loads (learner, pipeline) tuple
- `grimperium.DictStrAny` — type alias for dict[str, Any]
- `DeltaLearner.predict(X, y_pm7)` — returns H298_predicted array
- `FeaturePipeline.transform(df)` — returns feature matrix

**Verify GREEN:** `pytest tests/ml/test_predictor.py -v --tb=short` → 8 passed

---

## Step 2.5: Fix env var bug in `models_view.py`

**Bug:** `_MODEL_PATH` and `_DATA_PATH` (lines 26-29) are module-level constants that call `os.environ.get()` at **import time**. Any env var set after import is silently ignored.

**Fix:** Convert to instance methods on `ModelsView`:

1. **Remove** module-level `_MODEL_PATH` and `_DATA_PATH` (lines 26-29)
2. **Add** two instance methods to `ModelsView`:
   ```python
   def _get_model_path(self) -> Path:
       return Path(os.environ.get(
           "GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib"
       ))

   def _get_data_path(self) -> Path:
       return Path(os.environ.get(
           "GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"
       ))
   ```
3. **Replace all references** in the class:
   - `_MODEL_PATH` → `self._get_model_path()` (lines 53, 288, 311)
   - `_DATA_PATH` → `self._get_data_path()` (line 266)

---

## Step 3: Update `results_view.py`

### Changes:
1. **Replace imports:** Remove `mock_data.DIVERGENCE_STATS, MODELS` → add `os, Path, np, pd, load_model_metadata`
2. **Use instance methods for paths (NOT module constants):**
   ```python
   # Instance methods on ResultsView:
   def _get_model_path(self) -> Path:
       return Path(os.environ.get(
           "GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib"
       ))

   def _get_csv_path(self) -> Path:
       return Path(os.environ.get(
           "GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"
       ))
   ```
   `_DIVERGENCE_THRESHOLDS` stays as module constant (static data, no env dependency).
3. **New method `_load_real_model_row()`:** Loads metadata via `load_model_metadata(self._get_model_path())`, returns dict or None
4. **New method `_compute_divergence_stats()`:** Reads CSV, computes `abs(predicted - cbs)/abs(cbs) * 100`, bins into severity thresholds
5. **Rewrite `_render_model_comparison()`:** Use `_load_real_model_row()` — single row or "no model" message
6. **Rewrite `_render_divergence_analysis()`:** Use `_compute_divergence_stats()` — real stats or "run predict first" message
7. **Replace "Export Report" MenuOption** with `MenuOption(label="Predict Batch", value="predict_batch", icon=ICONS.get("calc", "⚡"))`
8. **Add `handle_action` for "predict_batch":** Import `predict_batch`, run it, display Rich Panel with stats
9. **Keep** "Detailed Metrics" and "Visualization Charts" as In Development

---

## Step 4: Verification

```bash
# 1. Full test suite
python -m pytest tests/ -q --tb=short

# 2. Type checking
mypy src/grimperium/ml/predictor.py --strict
mypy src/grimperium/cli/views/results_view.py --strict
mypy src/grimperium/cli/views/models_view.py --strict

# 3. Linting & formatting
ruff check src/grimperium/ml/predictor.py src/grimperium/cli/views/results_view.py src/grimperium/cli/views/models_view.py
black src/grimperium/ml/predictor.py src/grimperium/cli/views/results_view.py src/grimperium/cli/views/models_view.py --check

# 4. Integration test on real data
python -c "
from pathlib import Path
from grimperium.ml.predictor import predict_batch
df, stats = predict_batch(Path('data/thermo_pm7.csv'), Path('models/delta_learner_v1.joblib'), return_stats=True)
print(f'N predicted: {stats[\"n_predicted\"]}')
print(f'Mean delta: {stats[\"mean_delta\"]:.3f}')
print(f'Std delta: {stats[\"std_delta\"]:.3f}')
"
```

---

## Expected Final Report Format

```
TESTES PREDICTOR:    8 passed, 0 failed
TESTES TOTAL:        N passed, N skipped, 0 failed
MYPY predictor:      OK
MYPY results_view:   OK
RUFF:                OK
BLACK:               OK
INTEGRAÇÃO REAL:     N predições, stats printed
```
