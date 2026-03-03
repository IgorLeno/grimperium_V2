# Plan: Atomic CSV Save + Backup + Recovery

## Context

During PM7 batch execution, Ctrl-C can truncate `thermo_pm7.csv` to 0 bytes because `save_csv()` writes directly via `df.to_csv(self.csv_path)`, which truncates before writing. This CSV has 58+ columns of scientific data (RDKit descriptors, CREST conformers, MOPAC energies) that can represent hours/days of computation. The fix implements atomic write-then-rename, automatic backup, and recovery on load.

**Existing pattern to reuse**: `src/grimperium/crest_pm7/batch/detail_manager.py` already implements atomic writes for JSON using `tempfile.mkstemp` + `os.fsync` + `os.rename`. We'll follow the same approach.

---

## Files to Modify/Create

| File | Action |
|------|--------|
| `src/grimperium/core/csv_utils.py` | **CREATE** - Shared atomic CSV write utility |
| `src/grimperium/crest_pm7/batch/csv_manager.py` | **MODIFY** - Use atomic save + add recovery to load |
| `src/grimperium/cli/dataset_manager.py` | **MODIFY** - Use atomic save in `create_working_csv()` |
| `tests/test_csv_atomic_save.py` | **CREATE** - 6+ tests for atomic save/recovery |

---

## Step 1: Create `src/grimperium/core/csv_utils.py`

Shared utility function reusable by both `BatchCSVManager` and `DatasetManager`.

```python
def atomic_to_csv(path: Path, df: pd.DataFrame) -> None:
```

Algorithm:
1. Validate inputs (path not None, df not None)
2. Ensure parent directory exists
3. Write to temp file (same directory, `.tmp` suffix) using `tempfile.NamedTemporaryFile`
4. `df.to_csv(f, index=False)` + `f.flush()` + `os.fsync(f.fileno())`
5. If `path` exists and has size > 0: `shutil.copy2(path, path.with_suffix(path.suffix + '.bak'))`
6. `os.replace(tmp_path, path)` (atomic on same filesystem)
7. Cleanup temp on any exception

**Why temp file in same directory**: `os.replace` is atomic only on the same filesystem. Same pattern used in `detail_manager.py`.

---

## Step 2: Modify `BatchCSVManager.save_csv()`

Replace the body with a call to `atomic_to_csv`:

```python
def save_csv(self) -> None:
    if self.csv_path is None:
        raise RuntimeError("csv_path is None - cannot save CSV")
    if self.df is None:
        raise RuntimeError("No DataFrame loaded - call load_csv() first")
    atomic_to_csv(self.csv_path, self.df)
    LOG.info("Atomically saved %d molecules to %s", len(self.df), self.csv_path)
```

**Public signature unchanged** - all 11 callers continue working without modification.

---

## Step 3: Modify `BatchCSVManager.load_csv()` with recovery

Add three pieces of logic **before** the existing `pd.read_csv` call:

### 3a. Clean orphan `.tmp` files
```python
for tmp in self.csv_path.parent.glob("*.csv.tmp"):
    tmp.unlink()
    LOG.warning("Removing orphan temp file: %s", tmp)
```

### 3b. Detect empty/corrupted CSV and fall back to `.bak`
```python
backup_path = self.csv_path.with_suffix(self.csv_path.suffix + ".bak")

try:
    if self.csv_path.stat().st_size == 0:
        raise pd.errors.EmptyDataError("CSV file is 0 bytes")
    self.df = pd.read_csv(self.csv_path, dtype={...})
except (pd.errors.EmptyDataError, pd.errors.ParserError):
    # Try backup
    if backup_path.exists() and backup_path.stat().st_size > 0:
        LOG.error("CSV corrupted/empty, restoring from backup: %s", backup_path)
        shutil.copy2(backup_path, self.csv_path)
        self.df = pd.read_csv(self.csv_path, dtype={...})
        LOG.warning("Recovered %d molecules from backup", len(self.df))
    else:
        raise CSVCorruptedError(
            f"CSV and backup both corrupted/missing. "
            f"Restore from source: {self.csv_path}"
        )
```

### 3c. Add `CSVCorruptedError` exception class
Simple custom exception defined at module level in `csv_manager.py`.

**Existing validation logic** (column check, status normalization, duplicate mol_id check) remains unchanged after the read.

---

## Step 4: Modify `DatasetManager.create_working_csv()`

Replace `df.to_csv(self.working_csv, index=False)` with:
```python
from grimperium.core.csv_utils import atomic_to_csv
atomic_to_csv(self.working_csv, df)
```

Minimal change - the `mkdir` call above can stay since `atomic_to_csv` also ensures parent exists.

---

## Step 5: Create `tests/test_csv_atomic_save.py`

6 tests covering:

1. **`test_atomic_save_survives_interruption`** - Monkeypatch `DataFrame.to_csv` to raise during temp write. Assert original CSV intact.
2. **`test_atomic_save_preserves_all_columns`** - Full 58-column CSV survives save + reload with `.bak` created.
3. **`test_load_csv_recovers_from_backup`** - 0-byte CSV + valid `.bak` -> auto-recovery.
4. **`test_orphan_tmp_files_cleaned`** - Stale `.tmp` files removed on `load_csv()`.
5. **`test_both_csv_and_backup_corrupted`** - Raises `CSVCorruptedError` with actionable message.
6. **`test_calculated_fields_persist_across_cycles`** - 10 save/load cycles preserve float precision.

---

## Implementation Order

1. Create `csv_utils.py` with `atomic_to_csv()`
2. Add `CSVCorruptedError` to `csv_manager.py`
3. Modify `save_csv()` in `csv_manager.py`
4. Modify `load_csv()` in `csv_manager.py`
5. Modify `create_working_csv()` in `dataset_manager.py`
6. Create tests
7. Run quality gates

---

## Verification

```bash
# Tests pass
pytest tests/test_csv_atomic_save.py tests/test_csv_manager_retry.py -v

# Existing tests unbroken
pytest tests/ -v

# Type checking
mypy --strict src/grimperium/core/csv_utils.py src/grimperium/crest_pm7/batch/csv_manager.py

# Linting + formatting
ruff check src/grimperium/core/csv_utils.py src/grimperium/crest_pm7/batch/csv_manager.py src/grimperium/cli/dataset_manager.py
black --check src/grimperium/
```
