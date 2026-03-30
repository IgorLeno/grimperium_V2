# CBS SUSPECT Filter — PM7 Baseline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Exclude `cbs_quality_flag == "SUSPECT"` rows from the PM7 baseline analysis and document the known data issue.

**Architecture:** Add a module-level helper `_filter_suspect_rows()` in `databases_view.py` that encapsulates the filtering and count logic, then call it from `_handle_pm7_baseline()`. This keeps the helper unit-testable without mocking the full view. `trainer.py` is already safe — `data_loader.py:51` already filters `cbs_quality_flag == "OK"` before any training data reaches the model.

**Tech Stack:** pandas, Rich (console output), pytest

---

## Pre-check: Confirm trainer.py is already safe

Before any code changes, verify the existing guard in `data_loader.py`:

```python
# src/grimperium/ml/data_loader.py, line 51
df = df[(df["status"] == "OK") & (df["cbs_quality_flag"] == "OK")].copy()
```

This already blocks SUSPECT rows from training. No changes needed in `trainer.py` or `data_loader.py`.

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `src/grimperium/cli/views/databases_view.py` | **Modify** | Add `_filter_suspect_rows()` helper; call it in `_handle_pm7_baseline()` |
| `tests/cli/test_pm7_baseline.py` | **Modify** | Add 3 test scenarios for SUSPECT filtering |
| `docs/known_issues.md` | **Create** | Document the CBS SUSPECT data anomaly |

---

## Task 1 — Write failing tests for SUSPECT filtering

**Files:**
- Modify: `tests/cli/test_pm7_baseline.py`

These tests target `_filter_suspect_rows()`, which doesn't exist yet — they will fail with `ImportError`.

- [ ] **Step 1: Add the import and 3 test cases to `test_pm7_baseline.py`**

Add the import at the **top of the file** alongside the existing imports (around line 10):

```python
from grimperium.cli.views.databases_view import _compute_pm7_stats, _filter_suspect_rows
```

(Replace the existing `_compute_pm7_stats`-only import line.)

Then add the helper and test class at the **bottom** of the file:

```python
def _make_df_with_flag(
    cbs: list[float],
    pm7: list[float],
    flags: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        {"H298_cbs": cbs, "H298_pm7": pm7, "cbs_quality_flag": flags}
    )


class TestFilterSuspectRows:
    def test_suspect_rows_are_excluded_and_count_returned(self) -> None:
        df = _make_df_with_flag(
            cbs=[1.0, 2.0, -17307.0],
            pm7=[1.1, 2.1, -33.2],
            flags=["OK", "OK", "SUSPECT"],
        )
        filtered, count = _filter_suspect_rows(df)
        assert count == 1
        assert len(filtered) == 2
        assert (filtered["cbs_quality_flag"] == "SUSPECT").sum() == 0

    def test_no_suspect_rows_returns_full_df_and_zero_count(self) -> None:
        df = _make_df_with_flag(
            cbs=[1.0, 2.0],
            pm7=[1.1, 2.1],
            flags=["OK", "OK"],
        )
        filtered, count = _filter_suspect_rows(df)
        assert count == 0
        assert len(filtered) == 2

    def test_missing_cbs_quality_flag_column_returns_df_unchanged_and_zero_count(
        self,
    ) -> None:
        df = pd.DataFrame({"H298_cbs": [1.0, 2.0], "H298_pm7": [1.1, 2.1]})
        filtered, count = _filter_suspect_rows(df)
        assert count == 0
        assert len(filtered) == 2
```

- [ ] **Step 2: Run tests to confirm they fail with ImportError**

```bash
pytest tests/cli/test_pm7_baseline.py::TestFilterSuspectRows -v
```

Expected: `ImportError: cannot import name '_filter_suspect_rows'`

---

## Task 2 — Implement `_filter_suspect_rows` and wire into `_handle_pm7_baseline`

**Files:**
- Modify: `src/grimperium/cli/views/databases_view.py`

- [ ] **Step 3: Add `_filter_suspect_rows` as a module-level function**

Insert immediately after the `_compute_pm7_stats` function (around line 74):

```python
def _filter_suspect_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Remove rows flagged as SUSPECT CBS reference values.

    Args:
        df: DataFrame that may contain a ``cbs_quality_flag`` column.

    Returns:
        Tuple of (filtered_df, suspect_count). If the column is absent,
        returns (df, 0) unchanged — backward-compatible with datasets that
        predate the flag.
    """
    if "cbs_quality_flag" not in df.columns:
        return df, 0
    suspect_mask = df["cbs_quality_flag"] == "SUSPECT"
    suspect_count: int = int(suspect_mask.sum())
    return df[~suspect_mask].copy(), suspect_count
```

- [ ] **Step 4: Run the 3 new tests — they should now pass**

```bash
pytest tests/cli/test_pm7_baseline.py::TestFilterSuspectRows -v
```

Expected: 3 PASSED

- [ ] **Step 5: Wire `_filter_suspect_rows` into `_handle_pm7_baseline`**

In `_handle_pm7_baseline` (around line 387), replace:

```python
        valid = df.dropna(subset=["H298_pm7", "H298_cbs"])
        if len(valid) == 0:
            self.show_error("No valid rows with both PM7 and CBS values.")
            self.wait_for_enter()
            return

        stats = _compute_pm7_stats(valid)
```

with:

```python
        valid = df.dropna(subset=["H298_pm7", "H298_cbs"])
        valid, suspect_count = _filter_suspect_rows(valid)
        if len(valid) == 0:
            self.show_error("No valid rows with both PM7 and CBS values.")
            self.wait_for_enter()
            return

        stats = _compute_pm7_stats(valid)
```

Replace the `self.console.print()` blank-line call that follows `self.console.print(table)` (line 424) with:

```python
        if suspect_count > 0:
            self.console.print(
                f"[yellow]⚠ {suspect_count} molecule(s) with "
                "cbs_quality_flag=SUSPECT excluded from this analysis. "
                "These values originate from the CBS source dataset and are "
                "flagged as unreliable. Retained in CSV for traceability.[/yellow]"
            )
        self.console.print()
```

Note: if ruff's `ISC001` rule is enabled in the project, the implicit string concatenation in the f-string may be flagged. In that case, combine all four string literals into one by removing the line breaks between them.

- [ ] **Step 6: Run the full test suite for the CLI tests**

```bash
pytest tests/cli/ -v
```

Expected: all previously passing tests still pass; 3 new tests pass.

- [ ] **Step 7: Run quality gates**

```bash
black src/grimperium/cli/views/databases_view.py tests/cli/test_pm7_baseline.py
ruff check src/grimperium/cli/views/databases_view.py tests/cli/test_pm7_baseline.py
mypy src/grimperium/cli/views/databases_view.py --strict
```

Expected: no errors.

- [ ] **Step 8: Commit**

```bash
git add src/grimperium/cli/views/databases_view.py tests/cli/test_pm7_baseline.py
git commit -m "fix(cli): exclude cbs_quality_flag=SUSPECT from PM7 baseline analysis

13 molecules in thermo_pm7.csv have H298_cbs values in the range of -17k
to -145k kcal/mol (possible unit conversion error upstream). These were
inflating MARE from 6 kcal/mol to 757 kcal/mol and flipping R² to -0.007.

Adds _filter_suspect_rows() helper (tested) and wires it into
_handle_pm7_baseline(). trainer.py was already safe via data_loader.py:51."
```

---

## Task 3 — Create `docs/known_issues.md`

**Files:**
- Create: `docs/known_issues.md`

- [ ] **Step 9: Create the file**

```markdown
# Known Issues — Grimperium

## CBS_SUSPECT: Molecules with anomalous H298_cbs in source dataset

**Status:** Known, documented, NOT corrected (upstream data)
**Discovered:** March 2026 (Batch 12)
**How found:** PM7 Baseline — MARE = 757 kcal/mol vs P50 = 5.14 kcal/mol
(mean >> median by 147x is a classic outlier signal)

### Description

13 molecules in `thermo_pm7.csv` carry `H298_cbs` values in the range of
−17,000 to −145,000 kcal/mol. A CHON molecule with 7–11 heavy atoms has a
physically plausible H298_cbs range of roughly −300 to +50 kcal/mol.

Most probable cause: CBS total electronic energy in Hartrees stored without
conversion to kcal/mol (e.g., −232 Ha × 627.5 kcal/mol/Ha ≈ −145,620 kcal/mol).
Origin is upstream of this repository and is not traceable from the current import
history.

### Quantified Impact (without filter)

| Metric | Corrupted | Clean (13 rows removed) |
|---|---|---|
| MARE | 757 kcal/mol | 6.22 kcal/mol |
| Bias | +746 kcal/mol | −5.00 kcal/mol |
| R² | −0.0068 | 0.9845 |

### Treatment

`cbs_quality_flag = "SUSPECT"` marks these 13 rows in the CSV.

| File | Symbol | Filter |
|---|---|---|
| `databases_view.py` | `_handle_pm7_baseline` | `_filter_suspect_rows()` |
| `database_analyzer.py` | `_target_delta_stats` | `cbs_quality_flag == "OK"` |
| `database_analyzer.py` | `_top_delta_outliers` | `cbs_quality_flag == "OK"` |
| `ml/data_loader.py` | `load_ml_data` | `cbs_quality_flag == "OK"` |

Rows are **retained in the CSV for traceability** but excluded from all
analysis and training.

### Reproduce

```python
import pandas as pd
df = pd.read_csv("data/thermo_pm7.csv")
suspects = df[df["cbs_quality_flag"] == "SUSPECT"][
    ["mol_id", "smiles", "nheavy", "H298_cbs"]
]
print(suspects)  # 13 rows
```
```

- [ ] **Step 10: Verify file was created, commit**

```bash
git add docs/known_issues.md
git commit -m "docs: add known_issues.md documenting CBS_SUSPECT data anomaly"
```

---

## Acceptance Criteria

- [ ] `_handle_pm7_baseline()` shows MARE ≈ 6 kcal/mol and R² ≈ 0.98 on current dataset
- [ ] Yellow warning with SUSPECT count appears in the output
- [ ] `trainer.py` / `data_loader.py` confirmed already guards against SUSPECT (no changes needed)
- [ ] 3 new tests in `TestFilterSuspectRows` all pass
- [ ] `black`, `ruff`, `mypy --strict` clean
- [ ] `docs/known_issues.md` committed
