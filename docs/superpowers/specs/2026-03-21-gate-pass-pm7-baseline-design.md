# Design: gate_pass + PM7 Baseline Fix

**Date:** 2026-03-21
**Status:** Approved by user

---

## Bug 1 — `gate_pass` Always False

### Root Cause

`compute_all_metrics()` returns `{rmse, mae, r2, mape, max_error}` — no `gate_pass`.
`models_view.py` reads `test_metrics.get("gate_pass", False)`, which always returns `False`.

### Decision: Injection Point

Gate pass is computed in `models_view._handle_train()` / `_handle_retrain()` **before**
calling `save_model()`, not inside `trainer.train()`.

Rationale:
- `trainer.train()` stays a pure data function (raw metrics only); existing tests do not break.
- `gate_pass` IS persisted in the bundle (injected into `test_m` before `save_model()` call).
- Policy (thresholds) lives exclusively in `gate.py`; callers do not need to know thresholds.

### New Module: `src/grimperium/ml/gate.py`

```python
GATE_MAE_MAX  = 3.5   # kcal/mol
GATE_R2_MIN   = 0.97
GATE_RMSE_MAX = 5.0   # kcal/mol

def evaluate_gate(metrics: dict[str, float]) -> bool: ...
```

Gate passes when ALL three criteria are satisfied simultaneously.
`max_error` and `MAPE` are explicitly excluded from the gate (see Rationale below).

### Threshold Rationale

| Criterion | Value | Chemistry justification |
|---|---|---|
| MAE ≤ 3.5 kcal/mol | Strict | 33% margin above current 2.62; a naive mean-delta predictor scores ~5.15 kcal/mol. Anything above 3.5 represents clear regression. |
| R² ≥ 0.97 | Tight | Raises floor from 0.95 to 0.97; current model (0.9968) has 0.027 headroom. |
| RMSE ≤ 5.0 kcal/mol | Complementary | More stable than MAE for detecting outlier amplification; current 3.67 has 1.33 headroom. |
| max_error | Excluded | High values (16 kcal/mol) are driven by organic radicals (`[C]`, `[CH]`) — chemically hard but legitimate members of the training set. Gating on max_error would penalise good models for dataset composition. |
| MAPE | Excluded | H298 heats of formation cross zero; MAPE is numerically undefined for near-zero references (e.g., mol_00001 H298_cbs=0.152 → MAPE=2600%). |

### Changes

- **New:** `src/grimperium/ml/gate.py`
- **Modified:** `src/grimperium/cli/views/models_view.py` — `_handle_train()` and `_handle_retrain()` inject `gate_pass` into `test_m` before building the bundle.
- **No change:** `trainer.py`, `persistence.py`, `metrics.py`, `models_view._handle_test()`.
- **New tests:** `tests/ml/test_gate.py`

---

## Bug 2 — PM7 Baseline Analysis: Relative Error Metrics

### Root Cause

`_handle_pm7_baseline()` computes `re_pct = |pm7 - cbs| / |cbs| * 100`.
ΔHf°298 (heats of formation) cross zero by definition, making relative error
undefined for molecules with small |H298_cbs|. Example: mol_00001 has
H298_cbs=0.152 kcal/mol → RE% = 2616%. The metric is displayed but conveys no
physical information.

### Scale Compatibility Confirmed

For the 1550 molecules with both PM7 and CBS data, both columns are in kcal/mol
on the same reference frame (ΔHf°298, elements in standard state). The bias of
~−5 kcal/mol is PM7's systematic underestimation — physically correct.

MARE and bias ARE valid. Only the relative-error block is wrong.

### Fix

Remove the entire relative-error section. Replace with absolute-error distribution:

```
MARE              X.XX kcal/mol
Bias              X.XX kcal/mol
R²                0.XXXX

Absolute Error Distribution |PM7 − CBS|:
  Median (P50)    X.XX kcal/mol
  P90             X.XX kcal/mol
  P95             X.XX kcal/mol
  < 1 kcal/mol    XX%
  < 2 kcal/mol    XX%
  < 5 kcal/mol    XX%
```

### Refactor for Testability

Extract the stats computation into a module-level pure function:

```python
def _compute_pm7_stats(valid: pd.DataFrame) -> dict[str, float]: ...
```

`_handle_pm7_baseline()` calls this helper and renders the result.
Tests import `_compute_pm7_stats` directly without needing a live CLI view.

### Changes

- **Modified:** `src/grimperium/cli/views/databases_view.py`
  - Extract `_compute_pm7_stats(valid: pd.DataFrame) -> dict[str, float]`
  - `_handle_pm7_baseline()` calls the helper; rendering uses the new key set
- **New tests:** `tests/cli/test_pm7_baseline.py`

---

## TDD Order

1. Write `tests/ml/test_gate.py` (all red)
2. Implement `gate.py` (green)
3. Modify `models_view.py` to inject (no new tests needed — existing display paths already covered)
4. Write `tests/cli/test_pm7_baseline.py` (all red)
5. Extract `_compute_pm7_stats`, update `_handle_pm7_baseline()` (green)
6. Run full quality gate
