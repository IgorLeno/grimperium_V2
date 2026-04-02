# Gate Pass + PM7 Baseline Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `gate_pass` always showing ✗ No; remove meaningless relative-error metrics from PM7 Baseline Analysis.

**Architecture:** Two independent fixes. Bug 1: new `gate.py` module with thresholds + inject `gate_pass` into bundle metrics inside `models_view.py` before saving. Bug 2: extract pure stats helper from `_handle_pm7_baseline()`, swap RE% block for absolute-error distribution. TDD throughout — tests written before production code.

**Tech Stack:** Python 3.12, pytest, numpy, pandas, rich; no new dependencies.

---

## File Map

| Action | Path | Responsibility |
|--------|------|---------------|
| Create | `src/grimperium/ml/gate.py` | Threshold constants + `evaluate_gate(metrics) -> bool` |
| Create | `tests/ml/test_gate.py` | Unit tests for all gate logic |
| Modify | `src/grimperium/cli/views/models_view.py` | Inject `gate_pass` into `test_m` before `save_model()` in `_handle_train` and `_handle_retrain` |
| Modify | `src/grimperium/cli/views/databases_view.py` | Extract `_compute_pm7_stats()`, remove RE% block, add percentile block |
| Create | `tests/cli/test_pm7_baseline.py` | Unit tests for `_compute_pm7_stats` |

No changes to: `metrics.py`, `trainer.py`, `persistence.py`.

---

## Task 1: `gate.py` — Write Failing Tests

**Files:**
- Create: `tests/ml/test_gate.py`

- [ ] **Step 1: Create test file with all cases**

```python
# tests/ml/test_gate.py
"""Tests for grimperium.ml.gate — quality-gate evaluation."""

from __future__ import annotations

import pytest

from grimperium.ml.gate import GATE_MAE_MAX, GATE_R2_MIN, GATE_RMSE_MAX, evaluate_gate

_GOOD = {"mae": 2.5, "r2": 0.98, "rmse": 3.5, "mape": 8.0, "max_error": 15.0}


class TestGateThresholdConstants:
    def test_mae_max(self) -> None:
        assert GATE_MAE_MAX == 3.5

    def test_r2_min(self) -> None:
        assert GATE_R2_MIN == 0.97

    def test_rmse_max(self) -> None:
        assert GATE_RMSE_MAX == 5.0


class TestEvaluateGate:
    def test_passes_when_all_criteria_met(self) -> None:
        assert evaluate_gate(_GOOD) is True

    def test_fails_when_mae_exceeds_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "mae": 3.51}) is False

    def test_passes_when_mae_at_exact_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "mae": 3.5}) is True  # ≤, not <

    def test_fails_when_r2_below_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "r2": 0.969}) is False

    def test_passes_when_r2_at_exact_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "r2": 0.97}) is True  # ≥, not >

    def test_fails_when_rmse_exceeds_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "rmse": 5.01}) is False

    def test_passes_when_rmse_at_exact_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "rmse": 5.0}) is True  # ≤, not <

    def test_high_max_error_does_not_fail_gate(self) -> None:
        """max_error excluded — organic radicals drive it, not model quality."""
        assert evaluate_gate({**_GOOD, "max_error": 999.0}) is True

    def test_high_mape_does_not_fail_gate(self) -> None:
        """MAPE excluded — undefined for H298 near zero."""
        assert evaluate_gate({**_GOOD, "mape": 9999.0}) is True

    def test_all_criteria_fail_returns_false(self) -> None:
        assert evaluate_gate({"mae": 5.0, "r2": 0.90, "rmse": 7.0}) is False

    def test_missing_keys_treated_as_worst_case(self) -> None:
        """Missing mae/r2/rmse default to inf/0.0/inf → gate fails."""
        assert evaluate_gate({}) is False

    def test_current_model_passes(self) -> None:
        """Regression guard: model with MAE=2.62, R²=0.9968, RMSE=3.67 passes."""
        assert evaluate_gate({"mae": 2.62, "r2": 0.9968, "rmse": 3.67}) is True
```

- [ ] **Step 2: Run tests — expect ImportError (module does not exist yet)**

```bash
pytest tests/ml/test_gate.py -v
```

Expected: `ModuleNotFoundError: No module named 'grimperium.ml.gate'`

---

## Task 2: `gate.py` — Implement to Pass Tests

**Files:**
- Create: `src/grimperium/ml/gate.py`

- [ ] **Step 3: Create `gate.py`**

```python
# src/grimperium/ml/gate.py
"""Quality gate for trained DeltaLearner models.

Gate criteria reflect the minimum acceptable performance given 1 537 molecules
with delta mean ≈ 5 kcal/mol, std ≈ 6.45 kcal/mol:

    MAE  ≤ 3.5 kcal/mol  — naive mean-delta predictor scores ~5.15 kcal/mol;
                            3.5 marks a clear regression floor.
    R²   ≥ 0.97           — requires real variance explanation beyond the mean.
    RMSE ≤ 5.0 kcal/mol  — complements MAE; detects outlier amplification.

max_error and MAPE are intentionally excluded:
    - max_error is driven by organic radicals ([C], [CH]) — legitimate dataset members.
    - MAPE is numerically undefined when H298 is near zero (heats of formation cross zero).
"""

from __future__ import annotations

GATE_MAE_MAX: float = 3.5   # kcal/mol
GATE_R2_MIN: float = 0.97
GATE_RMSE_MAX: float = 5.0  # kcal/mol


def evaluate_gate(metrics: dict[str, float]) -> bool:
    """Return True if *metrics* satisfy all quality-gate criteria.

    Args:
        metrics: Dict with keys ``mae``, ``r2``, ``rmse`` (and optionally
                 others that are ignored). Missing keys default to worst-case
                 values (inf / 0.0) so a sparse dict fails the gate.

    Returns:
        True only when MAE ≤ GATE_MAE_MAX, R² ≥ GATE_R2_MIN,
        and RMSE ≤ GATE_RMSE_MAX simultaneously.
    """
    mae = float(metrics.get("mae", float("inf")))
    r2 = float(metrics.get("r2", 0.0))
    rmse = float(metrics.get("rmse", float("inf")))

    return mae <= GATE_MAE_MAX and r2 >= GATE_R2_MIN and rmse <= GATE_RMSE_MAX
```

- [ ] **Step 4: Run tests — expect all green**

```bash
pytest tests/ml/test_gate.py -v
```

Expected: `12 passed`

- [ ] **Step 5: Commit**

```bash
git add src/grimperium/ml/gate.py tests/ml/test_gate.py
git commit -m "feat(ml): add quality gate module with MAE/R2/RMSE thresholds"
```

---

## Task 3: Inject `gate_pass` in `models_view.py`

**Files:**
- Modify: `src/grimperium/cli/views/models_view.py`

`_handle_train()` (line 314) and `_handle_retrain()` (line 409) both unpack
`learner, train_m, test_m, pipeline = result` and then build the bundle.
Add one import and one injection line before each `save_model()` call.

- [ ] **Step 6: Add import at top of models_view.py**

Find the import block at the top of the file (around line 20). Insert after the
`from grimperium.ml.persistence import ...` line:

```python
from grimperium.ml.gate import evaluate_gate
```

- [ ] **Step 7: Inject in `_handle_train()` — line ~316**

Current code (after unpacking result):
```python
            learner, train_m, test_m, pipeline = result

            bundle: dict[str, Any] = {
```

Change to:
```python
            learner, train_m, test_m, pipeline = result
            test_m["gate_pass"] = evaluate_gate(test_m)

            bundle: dict[str, Any] = {
```

- [ ] **Step 8: Inject in `_handle_retrain()` — line ~411, and show gate result**

Same injection pattern — find the unpack in `_handle_retrain`:
```python
            learner, train_m, test_m, pipeline = result

            bundle: dict[str, Any] = {
```

Change to:
```python
            learner, train_m, test_m, pipeline = result
            test_m["gate_pass"] = evaluate_gate(test_m)

            bundle: dict[str, Any] = {
```

After the comparison table is printed (`self.console.print(comp)`), add a gate
result line so retrain shows the same pass/fail feedback as train:

```python
            gate_pass = test_m["gate_pass"]
            gate_icon = ICONS["success"] if gate_pass else ICONS["error"]
            self.console.print(
                f"  Gate Pass: {gate_icon} {'Yes' if gate_pass else 'No'}"
            )
```

- [ ] **Step 9: Verify existing trainer tests still pass (gate_pass NOT in train() return)**

```bash
pytest tests/ml/test_trainer.py -v
```

Expected: all green. `test_metrics_contain_expected_keys` still passes because
`train()` returns `{rmse, mae, r2, mape, max_error}` — gate injection happens
only in the view after `train()` returns.

- [ ] **Step 10: Commit**

```bash
git add src/grimperium/cli/views/models_view.py
git commit -m "fix(cli): inject gate_pass into test_m before saving model bundle"
```

---

## Task 4: PM7 Baseline — Write Failing Tests

**Files:**
- Create: `tests/cli/test_pm7_baseline.py`

`_compute_pm7_stats` will be a module-level function in `databases_view.py`.
Import it directly for testing — no CLI machinery needed.

- [ ] **Step 11: Create test file**

```python
# tests/cli/test_pm7_baseline.py
"""Tests for _compute_pm7_stats — PM7 baseline absolute-error distribution."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grimperium.cli.views.databases_view import _compute_pm7_stats


def _make_df(cbs: list[float], pm7: list[float]) -> pd.DataFrame:
    return pd.DataFrame({"H298_cbs": cbs, "H298_pm7": pm7})


class TestComputePm7StatsAbsoluteMetrics:
    def test_mare_is_mean_absolute_error(self) -> None:
        # Errors: 2, 3, 4 → MARE = 3.0
        df = _make_df([5.0, 5.0, 5.0], [3.0, 2.0, 1.0])
        result = _compute_pm7_stats(df)
        assert result["mare"] == pytest.approx(3.0)

    def test_bias_is_mean_pm7_minus_cbs(self) -> None:
        # pm7 - cbs: -2, -3, -4 → bias = -3.0
        df = _make_df([5.0, 5.0, 5.0], [3.0, 2.0, 1.0])
        result = _compute_pm7_stats(df)
        assert result["bias"] == pytest.approx(-3.0)

    def test_r2_is_coefficient_of_determination(self) -> None:
        # Perfect predictor: pm7 == cbs → r2 = 1.0
        df = _make_df([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
        result = _compute_pm7_stats(df)
        assert result["r2"] == pytest.approx(1.0)

    def test_n_is_row_count(self) -> None:
        df = _make_df([1.0, 2.0, 3.0], [0.5, 1.5, 2.5])
        result = _compute_pm7_stats(df)
        assert result["n"] == 3


class TestComputePm7StatsNoRelativeError:
    def test_mre_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])  # H298_cbs near zero → RE% would be huge
        result = _compute_pm7_stats(df)
        assert "mre_pct" not in result

    def test_mdre_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])
        assert "mdre_pct" not in _compute_pm7_stats(df)

    def test_max_re_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])
        assert "max_re_pct" not in _compute_pm7_stats(df)

    def test_std_re_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])
        assert "std_re_pct" not in _compute_pm7_stats(df)


class TestComputePm7StatsPercentiles:
    def test_p50_is_median_absolute_error(self) -> None:
        # Errors: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        df = _make_df(list(range(10)), [0.0] * 10)
        result = _compute_pm7_stats(df)
        assert result["p50"] == pytest.approx(np.median(range(10)))

    def test_p90_is_90th_percentile_absolute_error(self) -> None:
        df = _make_df(list(range(10)), [0.0] * 10)
        result = _compute_pm7_stats(df)
        assert result["p90"] == pytest.approx(np.percentile(range(10), 90))

    def test_p95_is_95th_percentile_absolute_error(self) -> None:
        df = _make_df(list(range(10)), [0.0] * 10)
        result = _compute_pm7_stats(df)
        assert result["p95"] == pytest.approx(np.percentile(range(10), 95))


class TestComputePm7StatsAbsoluteThresholds:
    def test_pct_lt_1_counts_errors_below_1_kcalmol(self) -> None:
        # Errors: 0.5, 1.5, 2.5, 4.5, 6.5 → only 0.5 < 1 → 20%
        df = _make_df([0.0, 0.0, 0.0, 0.0, 0.0], [-0.5, -1.5, -2.5, -4.5, -6.5])
        result = _compute_pm7_stats(df)
        assert result["pct_lt_1"] == pytest.approx(20.0)

    def test_pct_lt_2_counts_errors_below_2_kcalmol(self) -> None:
        # Errors: 0.5, 1.5, 2.5, 4.5, 6.5 → 0.5 and 1.5 < 2 → 40%
        df = _make_df([0.0, 0.0, 0.0, 0.0, 0.0], [-0.5, -1.5, -2.5, -4.5, -6.5])
        result = _compute_pm7_stats(df)
        assert result["pct_lt_2"] == pytest.approx(40.0)

    def test_pct_lt_5_counts_errors_below_5_kcalmol(self) -> None:
        # Errors: 0.5, 1.5, 2.5, 4.5, 6.5 → first four < 5 → 80%
        df = _make_df([0.0, 0.0, 0.0, 0.0, 0.0], [-0.5, -1.5, -2.5, -4.5, -6.5])
        result = _compute_pm7_stats(df)
        assert result["pct_lt_5"] == pytest.approx(80.0)

    def test_near_zero_h298_cbs_does_not_distort_absolute_metrics(self) -> None:
        """Near-zero H298_cbs that caused RE% blowup is harmless for absolute stats."""
        # mol_00001-style: H298_cbs=0.152, H298_pm7=-3.83 → |error|=3.98
        df = _make_df([0.152, -12.641], [-3.830, -14.660])
        result = _compute_pm7_stats(df)
        assert result["mare"] == pytest.approx((3.982 + 2.019) / 2, abs=0.01)
        assert "mre_pct" not in result
```

- [ ] **Step 12: Run tests — expect ImportError**

```bash
pytest tests/cli/test_pm7_baseline.py -v
```

Expected: `ImportError: cannot import name '_compute_pm7_stats'`

---

## Task 5: Extract `_compute_pm7_stats` + Fix Display

**Files:**
- Modify: `src/grimperium/cli/views/databases_view.py`

The function `_handle_pm7_baseline` currently occupies lines ~333-420.
Extract the calculation into a module-level helper, replace the RE% table rows
with percentile rows.

- [ ] **Step 13: Add `_compute_pm7_stats` above `DatabasesView` class**

Insert this function before the `class DatabasesView` line (after the imports block).
`mae` and `r2_score` are already imported from `grimperium.core.metrics` (line 33).

```python
def _compute_pm7_stats(valid: pd.DataFrame) -> dict[str, float]:
    """Compute PM7 vs CBS absolute-error statistics from a validated DataFrame.

    Args:
        valid: DataFrame with non-null ``H298_cbs`` and ``H298_pm7`` columns.

    Returns:
        Dict with keys: mare, bias, r2, n, p50, p90, p95,
        pct_lt_1, pct_lt_2, pct_lt_5.

    Note:
        Relative-error metrics (MRE%, MdRE%) are intentionally absent.
        H298 heats of formation cross zero, making percentage errors
        numerically undefined for molecules with small |H298_cbs|.
    """
    h298_cbs = valid["H298_cbs"].to_numpy(dtype=float)
    h298_pm7 = valid["H298_pm7"].to_numpy(dtype=float)

    abs_err = np.abs(h298_pm7 - h298_cbs)
    n = len(abs_err)

    return {
        "mare": float(mae(h298_cbs, h298_pm7)),
        "bias": float(np.mean(h298_pm7 - h298_cbs)),
        "r2": float(r2_score(h298_cbs, h298_pm7)),
        "n": n,
        "p50": float(np.median(abs_err)),
        "p90": float(np.percentile(abs_err, 90)),
        "p95": float(np.percentile(abs_err, 95)),
        "pct_lt_1": float(np.sum(abs_err < 1.0) / n * 100),
        "pct_lt_2": float(np.sum(abs_err < 2.0) / n * 100),
        "pct_lt_5": float(np.sum(abs_err < 5.0) / n * 100),
    }
```

- [ ] **Step 14: Replace the calculation block inside `_handle_pm7_baseline()`**

Find the block starting at line ~358 (`h298_cbs = valid["H298_cbs"].to_numpy()`).
It ends just before `# Display table` (~line 385). Replace that entire block with:

```python
        stats = _compute_pm7_stats(valid)
```

- [ ] **Step 15: Replace the RE% table rows with absolute-error distribution rows**

In the same method, find the `# Display table` section (the `Table(...)` and
`table.add_row(...)` calls). The current table rows are:

```python
        table.add_row("MARE (kcal/mol)", f"{mare:.4f}")
        table.add_row("Bias (kcal/mol)", f"{bias:.4f}")
        table.add_row("R\u00b2", f"{r2:.4f}")
        table.add_row("", "")
        table.add_row("MRE%", f"{mre_pct:.2f}%")
        table.add_row("MdRE%", f"{mdre_pct:.2f}%")
        table.add_row("Std RE% (ddof=0)", f"{std_re_pct:.4f}%")
        table.add_row("Max RE%", f"{max_re_pct:.2f}%")
        table.add_row("", "")
        table.add_row("RE% < 1%", f"{pct_lt_1:.1f}%")
        table.add_row("RE% < 5%", f"{pct_lt_5:.1f}%")
        table.add_row("RE% < 10%", f"{pct_lt_10:.1f}%")
        table.add_row("", "")
```

Replace with:

```python
        table.add_row("MARE (kcal/mol)", f"{stats['mare']:.4f}")
        table.add_row("Bias (kcal/mol)", f"{stats['bias']:.4f}")
        table.add_row("R\u00b2", f"{stats['r2']:.4f}")
        table.add_row("", "")
        table.add_row("Absolute Error Distribution", "")
        table.add_row("  Median |error| (P50)", f"{stats['p50']:.3f} kcal/mol")
        table.add_row("  P90", f"{stats['p90']:.3f} kcal/mol")
        table.add_row("  P95", f"{stats['p95']:.3f} kcal/mol")
        table.add_row("", "")
        table.add_row("  |error| < 1 kcal/mol", f"{stats['pct_lt_1']:.1f}%")
        table.add_row("  |error| < 2 kcal/mol", f"{stats['pct_lt_2']:.1f}%")
        table.add_row("  |error| < 5 kcal/mol", f"{stats['pct_lt_5']:.1f}%")
        table.add_row("", "")
```

- [ ] **Step 16: Remove entire relative-error block and fix interpretation panel**

After Step 14 (`stats = _compute_pm7_stats(valid)` replaced the CBS/PM7
array construction block), there is still a large dead block in the method body.
Remove it entirely — from the `nonzero_mask` assignment down through
the closing `else: mre_pct = mdre_pct = ... = 0.0` branch (the full `if/else`
on whether `len(cbs_nz) > 0`). This whole section no longer exists in the
new code. Leaving any fragment of it will cause a `NameError` or `SyntaxError`
at runtime.

Old block to delete (lines ~367–383):
```python
        # Relative metrics (exclude H298_cbs == 0)
        nonzero_mask = valid["H298_cbs"] != 0.0
        cbs_nz = valid.loc[nonzero_mask, "H298_cbs"].to_numpy()
        pm7_nz = valid.loc[nonzero_mask, "H298_pm7"].to_numpy()

        if len(cbs_nz) > 0:
            re_pct = np.abs(pm7_nz - cbs_nz) / np.abs(cbs_nz) * 100
            mre_pct = float(np.mean(re_pct))
            mdre_pct = float(np.median(re_pct))
            std_re_pct = float(np.std(re_pct, ddof=0))
            max_re_pct = float(np.max(re_pct))
            n_re = len(re_pct)
            pct_lt_1 = float(np.sum(re_pct < 1.0) / n_re * 100)
            pct_lt_5 = float(np.sum(re_pct < 5.0) / n_re * 100)
            pct_lt_10 = float(np.sum(re_pct < 10.0) / n_re * 100)
        else:
            mre_pct = mdre_pct = std_re_pct = max_re_pct = 0.0
            pct_lt_1 = pct_lt_5 = pct_lt_10 = 0.0
```

Delete this block completely. After Step 14, the method body should jump
directly from `stats = _compute_pm7_stats(valid)` to `# Display table`.

**Also update the closing table row** (line ~410) that uses `len(valid)`:
```python
# Old:
table.add_row("Molecules analyzed", f"{len(valid):,}")
# New:
table.add_row("Molecules analyzed", f"{stats['n']:,}")
```

**Also update the interpretation panel** (lines ~417-424).
The panel references `mare`, `bias`, and `mre_pct` directly — `mre_pct` will
cause a `NameError` after the refactor. Replace the entire `interpretation`
string with:

```python
        interpretation = (
            f"[bold]PM7 Baseline Error Context:[/bold]\n\n"
            f"\u2022 PM7 mean absolute error: {stats['mare']:.4f} kcal/mol vs CBS reference\n"
            f"\u2022 Systematic bias: {stats['bias']:+.4f} kcal/mol "
            f"({'PM7 overestimates' if stats['bias'] > 0 else 'PM7 underestimates' if stats['bias'] < 0 else 'unbiased'})\n"
            f"\u2022 Median absolute error (P50): {stats['p50']:.3f} kcal/mol\n"
            f"\u2022 This is the error that delta-learning needs to correct"
        )
```

- [ ] **Step 17: Run PM7 baseline tests**

```bash
pytest tests/cli/test_pm7_baseline.py -v
```

Expected: all green.

- [ ] **Step 18: Run full trainer + gate tests to check no regressions**

```bash
pytest tests/ml/ -v
```

Expected: all green including `test_metrics_contain_expected_keys`.

- [ ] **Step 19: Commit**

```bash
git add src/grimperium/cli/views/databases_view.py tests/cli/test_pm7_baseline.py
git commit -m "fix(cli): replace PM7 baseline RE% with absolute error distribution"
```

---

## Task 6: Quality Gate

- [ ] **Step 20: Run linter and type checker**

```bash
black src/grimperium/ml/gate.py tests/ml/test_gate.py \
      src/grimperium/cli/views/databases_view.py \
      src/grimperium/cli/views/models_view.py \
      tests/cli/test_pm7_baseline.py
ruff check src/ tests/
mypy src/ --strict
```

Expected: no errors.

- [ ] **Step 21: Run targeted test suite**

```bash
pytest tests/ml/ tests/cli/ -v --cov=src/grimperium/ml/gate \
       --cov=src/grimperium/cli/views/databases_view \
       --cov=src/grimperium/cli/views/models_view
```

Expected: all green.

- [ ] **Step 22: Final commit (if any formatting fixes were needed)**

```bash
git add -p
git commit -m "chore: apply black/ruff formatting to new files"
```

---

## Verification Checklist

After all tasks complete:

- [ ] `gate_pass = True` when training with `data/thermo_pm7.csv` (MAE=2.62 < 3.5, R²=0.9968 > 0.97, RMSE=3.67 < 5.0)
- [ ] Existing bundle `models/delta_learner_v1.joblib` still loads correctly (no `gate_pass` key → `_handle_test()` shows `False` as before; next retrain will show `True`)
- [ ] PM7 Baseline Analysis in CLI shows no RE%, shows P50/P90/P95 rows
- [ ] `tests/ml/test_trainer.py::test_metrics_contain_expected_keys` still passes (unchanged)
