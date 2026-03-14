# Plan: DatabaseAnalyzer Bug Fix + cbs_quality_flag

## Context

After running the new DatabaseAnalyzer on the PM7 database (1009 OK molecules), two issues surfaced:

1. **False positive pipeline failures:** The analyzer's `_pipeline_health()` compares `crest_status != "ok"` and `mopac_status != "ok"` (lowercase), but the actual CSV values are `"SUCCESS"` for crest and `"OK"` (uppercase) for mopac. This produces 1021 false "failed" counts.

2. **6 molecules with suspect CBS reference values:** Their `H298_cbs` values produce implausibly large deltas (-17k to -145k kcal/mol), likely unit errors in the reference data. These need to be flagged for exclusion from ML training without deleting the rows.

---

## BATCH 1: Fix `_pipeline_health()` status comparisons

**File:** `src/grimperium/crest_pm7/database_analyzer.py` (lines 469-494)

### Current (broken)

```python
# Line 484-486: compares against "ok" — misses "SUCCESS"
result["crest_failed"] = int(
    ((df["crest_status"] != "ok") & has_status).sum()
)
# Line 490-492: compares against "ok" — misses "OK" (uppercase)
result["mopac_failed"] = int(
    ((df["mopac_status"] != "ok") & has_status).sum()
)
```

### Fix

Replace exact string comparisons with case-insensitive checks that recognize actual pipeline values:

```python
# crest_status: healthy values are "ok" or "success" (case-insensitive)
if "crest_status" in df.columns:
    has_status = df["crest_status"].notna() & (df["crest_status"] != "")
    healthy = df["crest_status"].str.lower().isin({"ok", "success"})
    result["crest_failed"] = int((~healthy & has_status).sum())

# mopac_status: healthy value is "ok" (case-insensitive)
if "mopac_status" in df.columns:
    has_status = df["mopac_status"].notna() & (df["mopac_status"] != "")
    healthy = df["mopac_status"].str.lower() == "ok"
    result["mopac_failed"] = int((~healthy & has_status).sum())
```

### Additional: scope to OK rows only

The user's plan specifies that pipeline health should operate on `status == "OK"` rows only. Currently `_pipeline_health()` receives the full `df`. Change the call site in `analyze()` (line ~203):

```python
# Before:
pipeline = self._pipeline_health(df)
# After:
pipeline = self._pipeline_health(ok_df)
```

This prevents counting PENDING/RUNNING rows (which naturally have empty crest/mopac_status) in the denominator.

---

## BATCH 2: Add `cbs_quality_flag` column

### 2a. Schema: Add field to Pydantic model

**File:** `src/grimperium/crest_pm7/batch/models.py` (after line 209)

```python
cbs_quality_flag: str = Field(
    default="", description='CBS reference quality: "OK", "SUSPECT", or ""'
)
```

### 2b. Schema: Add to CSV column groups

**File:** `src/grimperium/crest_pm7/batch/csv_manager.py` (line 60-68)

Add `"cbs_quality_flag"` to `MOLECULAR_PROPERTIES_COLUMNS` after `"abs_diff_%"`.

### 2c. Auto-detection in `csv_enhancements.py`

**File:** `src/grimperium/crest_pm7/csv_enhancements.py`

**Verified internal API (from `csv_manager.py`):**
- `csv_manager._ensure_loaded() -> pd.DataFrame` (line 343) — returns the DataFrame, loading if needed
- `csv_manager._get_row_index(mol_id: str) -> int` (line 393) — returns integer index, raises `KeyError` if not found
- `csv_manager._update_extra_fields(mol_id, field_updates)` (line 858) — writes `{col: val}` dict to CSV, skips columns not in `df.columns`

**Step 1:** Add a static method to `CSVManagerExtensions`:

```python
@staticmethod
def compute_cbs_quality_flag(
    h298_cbs: float | None,
    h298_pm7: float | None,
    nheavy: int | None,
) -> str:
    """Flag CBS reference if delta is implausibly large.

    Threshold: |delta| > 500 * nheavy suggests unit error.
    Returns "OK" or "SUSPECT".
    """
    if h298_cbs is None or h298_pm7 is None or not nheavy:
        return "OK"
    delta = abs(h298_cbs - h298_pm7)
    if delta > 500 * nheavy:
        return "SUSPECT"
    return "OK"
```

**Step 2:** In `update_molecule_with_mopac_results()`, **after** building the `updates` dict (after line 288) and **before** the `_update_extra_fields` call (line 291), insert:

```python
# CBS quality flag — preserve manual SUSPECT, auto-detect on first write
try:
    df = csv_manager._ensure_loaded()
    idx = csv_manager._get_row_index(mol_id)
    current_flag = ""
    if "cbs_quality_flag" in df.columns:
        current_flag = str(df.at[idx, "cbs_quality_flag"] or "")
    if current_flag != "SUSPECT":
        nheavy_val = None
        if "nheavy" in df.columns:
            raw = df.at[idx, "nheavy"]
            if pd.notna(raw):
                nheavy_val = int(raw)
        updates["cbs_quality_flag"] = CSVManagerExtensions.compute_cbs_quality_flag(
            h298_cbs, h298_pm7, nheavy_val
        )
except (KeyError, ValueError, TypeError):
    pass  # mol_id not found or nheavy unreadable — skip flag
```

**Note:** `_update_extra_fields` (line 291) will only write `cbs_quality_flag` if the column exists in `df.columns`. Since we're adding it to `MOLECULAR_PROPERTIES_COLUMNS` (BATCH 2b), it will be present in new CSVs. For existing CSVs, the migration script (BATCH 2d) adds the column.

### 2d. Migration script (inline execution)

Run once to add the column to the existing CSV:

```python
import pandas as pd

CSV_PATH = "data/thermo_pm7.csv"
SUSPECT_IDS = ["mol_01397", "mol_00809", "mol_01167", "mol_00862", "mol_01152", "mol_00439"]

df = pd.read_csv(CSV_PATH, dtype=str)
df["cbs_quality_flag"] = ""
df.loc[df["status"] == "OK", "cbs_quality_flag"] = "OK"
df.loc[df["mol_id"].isin(SUSPECT_IDS), "cbs_quality_flag"] = "SUSPECT"
df.to_csv(CSV_PATH, index=False)
print(f"SUSPECT count: {(df['cbs_quality_flag'] == 'SUSPECT').sum()}")
print(f"OK count: {(df['cbs_quality_flag'] == 'OK').sum()}")
```

Expected output: `SUSPECT count: 6, OK count: 1003`

---

## BATCH 3: Report `cbs_quality_flag` in DatabaseAnalyzer

### 3a. Add fields to `AnalysisReport`

**File:** `src/grimperium/crest_pm7/database_analyzer.py`

After `ok_count` (line ~30), add:

```python
cbs_suspect_count: int     # cbs_quality_flag == "SUSPECT"
cbs_ok_count: int          # cbs_quality_flag == "OK" (available for training)
```

### 3b. Compute in `analyze()`

After the `ok_count` computation block, add:

```python
cbs_suspect = 0
cbs_ok = 0
if "cbs_quality_flag" in ok_df.columns:
    cbs_suspect = int((ok_df["cbs_quality_flag"] == "SUSPECT").sum())
    cbs_ok = int((ok_df["cbs_quality_flag"] == "OK").sum())
if cbs_suspect > 0:
    alerts.append(
        f"⚠️ {cbs_suspect} molecules with SUSPECT CBS reference "
        "(excluded from training)"
    )
```

Pass `cbs_suspect_count=cbs_suspect, cbs_ok_count=cbs_ok` to the `AnalysisReport` constructor.

### 3c. Render in Status Overview panel

**File:** `src/grimperium/cli/views/databases_view.py`

In `render_analysis_report()`, in the Status Overview panel (after orphan_running block), add:

```python
if report.cbs_ok_count > 0 or report.cbs_suspect_count > 0:
    total_flagged = report.cbs_ok_count + report.cbs_suspect_count
    ok_pct = report.cbs_ok_count / total_flagged * 100 if total_flagged else 0
    status_lines.append("")
    status_lines.append("[bold]CBS Quality:[/bold]")
    status_lines.append(
        f"  Available for training: {report.cbs_ok_count:,} ({ok_pct:.1f}%)"
    )
    if report.cbs_suspect_count > 0:
        status_lines.append(
            f"  [yellow]SUSPECT (excluded): {report.cbs_suspect_count}[/yellow]"
        )
```

---

## BATCH 4: Tests

### Update `tests/cli/test_database_analyzer.py`

Add 3 tests:

1. **`test_crest_status_success_not_failed`** — DataFrame with `status="OK"`, `crest_status="SUCCESS"` → `crest_failed_count == 0`
2. **`test_mopac_status_ok_uppercase_not_failed`** — DataFrame with `status="OK"`, `mopac_status="OK"` (uppercase) → `mopac_failed_count == 0`
3. **`test_cbs_suspect_flag`** — DataFrame with 10 OK rows, 2 with `cbs_quality_flag="SUSPECT"` → `cbs_suspect_count == 2`, `cbs_ok_count == 8`, alert contains "SUSPECT"

---

## Files Summary

| File | Change |
|------|--------|
| `src/grimperium/crest_pm7/database_analyzer.py` | Fix `_pipeline_health()` status comparison; add `cbs_suspect_count`/`cbs_ok_count` to report; scope pipeline health to OK rows |
| `src/grimperium/crest_pm7/csv_enhancements.py` | Add `compute_cbs_quality_flag()` + integrate into `update_molecule_with_mopac_results()` |
| `src/grimperium/crest_pm7/batch/models.py` | Add `cbs_quality_flag` field to `BatchRowCSV` |
| `src/grimperium/crest_pm7/batch/csv_manager.py` | Add `"cbs_quality_flag"` to `MOLECULAR_PROPERTIES_COLUMNS` |
| `src/grimperium/cli/views/databases_view.py` | Add CBS Quality sub-section to Status Overview panel |
| `tests/cli/test_database_analyzer.py` | Add 3 new tests |
| `data/thermo_pm7.csv` | Migration: add column, flag 6 SUSPECT molecules |

### Files NOT modified
- `src/grimperium/crest_pm7/batch/*.py` (except `csv_manager.py` and `models.py`)
- `src/grimperium/cli/database_registry.py`
- `src/grimperium/cli/mock_data.py`

---

## Execution Order

```
BATCH 1 (analyzer fix)  →  BATCH 2 (cbs_quality_flag)  →  BATCH 3 (report)  →  BATCH 4 (tests)
```

All sequential — each batch depends on the previous.

---

## Verification

1. Run migration script → output: `SUSPECT count: 6, OK count: 1003`
2. `python -m pytest tests/cli/test_database_analyzer.py -v` — all tests pass (13 existing + 3 new)
3. `python -m pytest tests/ -v` — zero regressions
4. `ruff check src/grimperium/crest_pm7/database_analyzer.py` — clean
5. `ruff check src/grimperium/crest_pm7/csv_enhancements.py` — clean
6. `ruff check src/grimperium/crest_pm7/batch/csv_manager.py` — clean
7. Manual: run Analyze Database on PM7 → `crest_failed_count = 0`, `mopac_failed_count = 0`
8. Manual: Status Overview shows CBS Quality with 1003 available, 6 SUSPECT
