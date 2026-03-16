# Phase G: Visualization Charts — Implementation Plan

## Context

Phase G adds 3 analysis charts (Parity Plot, Delta Histogram, Residuals Plot) to the Grimperium CLI, accessible via Results → Visualization Charts. Currently this menu option is disabled with "In Development". The charts visualize ML model predictions vs CBS reference values from the CSV dataset.

**Goal:** Generate publication-quality dark-themed PNG charts from CLI, following TDD.

---

## Files to Create

| File | Purpose |
|------|---------|
| `src/grimperium/ml/charts.py` | Chart generation logic (3 plots + dataclass result) |
| `tests/ml/test_charts.py` | 9 tests for chart generation |

## Files to Modify

| File | Change |
|------|--------|
| `tests/ml/conftest.py` | Add `csv_with_predictions_for_charts` fixture |
| `src/grimperium/ml/__init__.py` | Export `generate_charts`, `ChartGenerationResult`, `predict_batch`, `load_model_metadata`, `save_model`, `load_model` |
| `src/grimperium/cli/views/results_view.py` | Enable charts menu, implement `_handle_charts()`, update `handle_action()` dispatch |
| `.gitignore` | Add `reports/charts/*.png` |

## Do NOT Touch

- `styles.py`, `predictor.py`, `persistence.py`, `trainer.py`, `evaluator.py`, `data_loader.py`, `features.py`, `models_view.py`, `mock_data.py`

---

## Step-by-Step Implementation

### Step 1 — RED Phase: Write Tests

**1a.** Add fixture to `tests/ml/conftest.py` (append, don't rewrite):
- `csv_with_predictions_for_charts(tmp_path, synthetic_csv_path_mixed)` — trains a synthetic model, saves it, runs `predict_batch`, returns CSV path with predictions

**1b.** Create `tests/ml/test_charts.py` with 9 tests per spec:
- `test_returns_chart_result_object`
- `test_parity_plot_created`
- `test_delta_histogram_created`
- `test_residuals_plot_created`
- `test_output_dir_created_if_missing`
- `test_raises_on_missing_csv`
- `test_raises_if_no_predictions` (uses `synthetic_csv_path` — no predictions column)
- `test_chart_result_has_n_points`
- `test_chart_result_has_metrics`

**Verify RED:** `python -m pytest tests/ml/test_charts.py -v --tb=short` → 9 ImportError failures

### Step 2 — GREEN Phase: Implement `charts.py`

Create `src/grimperium/ml/charts.py`:

- **`ChartGenerationResult`** — `@dataclass` with fields: `parity_plot: Path`, `delta_histogram: Path`, `residuals_plot: Path`, `n_points: int`, `rmse: float`, `r2: float`
- **`generate_charts(csv_path, output_dir)`** — main entry point:
  - `matplotlib.use("Agg")` before any pyplot import (non-interactive backend)
  - Validates CSV exists, has `H298_predicted` column, has non-NaN rows
  - Computes RMSE and R² from plotted data
  - Calls 3 private plot functions
  - Returns `ChartGenerationResult`

- **`_plot_parity(y_cbs, y_pred, rmse, r2, output_dir)`** — scatter, y=x line, R²/RMSE annotation, cyan `#00D9FF`, dark bg, 8x7, dpi=150
- **`_plot_delta_histogram(delta, output_dir)`** — histogram bins=50, magenta `#FF00FF`, alpha=0.7, zero+mean lines, dark bg, 9x6, dpi=150
- **`_plot_residuals(y_pred, residuals, output_dir)`** — scatter, y=0 and ±5 lines, green `#00FF80`, alpha=0.5, dark bg, 9x6, dpi=150

**Key detail:** Each plot function uses `plt.style.use("dark_background")` and sets `facecolor="#1a1a1a"` on figure and axes.

**Verify GREEN:** `python -m pytest tests/ml/test_charts.py -v --tb=short` → 9/9 passed

### Step 3 — Update `ml/__init__.py`

Add exports:
```python
from grimperium.ml.charts import ChartGenerationResult, generate_charts
from grimperium.ml.persistence import load_model, load_model_metadata, save_model
from grimperium.ml.predictor import predict_batch
```

Update `__all__` to include all new exports (alphabetically sorted).

### Step 4 — Update `results_view.py`

**4a.** Add `_get_charts_dir()` method returning `Path(os.environ.get("GRIMPERIUM_CHARTS_DIR", "reports/charts"))`

**4b.** Enable menu: Remove `disabled=True` and `disabled_reason` from "Visualization Charts" `MenuOption`, use `icon=ICONS.get("results", "📊")`

**4c.** Update `handle_action()`:
- `"charts"` → `self._handle_charts(); return None`
- `"detailed"` → `self.show_in_development("Detailed Metrics"); return None`
- Remove the combined `if action in ["detailed", "charts"]` block

**4d.** Implement `_handle_charts()`:
- Lazy import of `generate_charts` and `ChartGenerationResult`
- Shows "Generating charts..." status
- Calls `generate_charts(csv_path, charts_dir)`
- Displays Panel with n_points, RMSE, R², file paths
- Error handling for `FileNotFoundError`, `ValueError`, generic `Exception`
- Ends with `self.wait_for_enter()`

### Step 5 — Update `.gitignore`

Add `reports/charts/*.png` after existing data/model patterns.

### Step 6 — Verification

1. `python -m pytest tests/ -q --tb=short` — full suite passes
2. `mypy src/grimperium/ml/charts.py --strict`
3. `mypy src/grimperium/cli/views/results_view.py --strict`
4. `ruff check` + `black --check` on changed files
5. Integration test: generate charts from real `data/thermo_pm7.csv` → verify 3 PNGs created

---

## Key Reuse Points

| Existing Code | Location | Reuse |
|--------------|----------|-------|
| `_get_csv_path()` | `results_view.py:53-59` | Reuse for chart CSV path |
| `predict_batch()` | `predictor.py` | Used in test fixture |
| `save_model()` | `persistence.py` | Used in test fixture |
| `train(return_pipeline=True)` | `trainer.py` | Used in test fixture |
| `COLORS`, `ICONS` | `styles.py` | Used in `_handle_charts()` display |
| `Panel` from rich | Already imported in `results_view.py` | Used for chart result display |
| `trained_model_fixture` pattern | `conftest.py:157-183` | Pattern for new fixture |
