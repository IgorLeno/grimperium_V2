# Implementation Plan: Fix CLI Bug + Documentation Sync

**Created:** 2026-01-19
**Status:** Ready for Approval
**Estimated Time:** 30-45 minutes
**Priority:** HIGH (CLI bug is blocking PM7 calculations)

---

## 🎯 Objectives

1. **Fix Critical CLI Bug**: Remove hardcoded test paths preventing CREST PM7 batch calculations
2. **Sync Documentation**: Update docs to reflect current codebase state (dataset names, CLI structure)

---

## 📋 Problem Summary

### CLI Bug (CRITICAL)
The user encounters this error when trying to calculate PM7 values:
```
✗ Batch CSV not found: data/test_batch_final.csv
Create a batch CSV with columns: mol_id, smiles, nheavy, status
```

**Root Cause:** `databases_view.py:580` uses hardcoded test path instead of production CSV.

### Documentation Drift (HIGH)
Several documentation files reference:
- Old dataset names (`thermo_cbs_clean.csv`, `thermo_batch_final.csv`)
- Missing new CLI view (`batch_view.py`)
- Stale line number references

---

## 🔧 Implementation Plan

### PART 1: Fix CLI Bug (Critical)

#### File: `src/grimperium/cli/views/databases_view.py`

**Change 1: Line 580 - Fix batch CSV path**
```python
# BEFORE (line 580):
csv_path = Path("data/test_batch_final.csv")

# AFTER:
csv_path = DATA_DIR / "thermo_pm7.csv"
```

**Change 2: Line 581 - Fix conformer details directory**
```python
# BEFORE (line 581):
detail_dir = Path("data/molecules_pm7/conformer_details")

# AFTER:
detail_dir = DATA_DIR / "molecules_pm7" / "conformer_details"
```

**Change 3: Line 364 - Fix working CSV reference in _refresh_database**
```python
# BEFORE (line 364):
manager = DatasetManager(
    source_csv=DATA_DIR / "thermo_cbs_chon.csv",
    working_csv=DATA_DIR / "thermo_batch_final.csv",  # ❌ Old name
)

# AFTER:
manager = DatasetManager(
    source_csv=DATA_DIR / "thermo_cbs_chon.csv",
    working_csv=DATA_DIR / "thermo_pm7.csv",  # ✅ Current name
)
```

**Imports Check:**
Verify `DATA_DIR` is imported at top of file (should already be present from line 20).

---

#### File: `src/grimperium/cli/views/batch_view.py`

**Change 4: Line 41 - Fix batch tracking CSV path**
```python
# BEFORE (line 41):
DEFAULT_CSV_PATH = Path("data/batch_tracking.csv")

# AFTER:
DEFAULT_CSV_PATH = DATA_DIR / "batch_tracking.csv"
```

**Required:** Add import at top of file:
```python
from grimperium.cli.constants import DATA_DIR
```

---

### PART 2: Documentation Sync

#### File: `docs/CREST_INTEGRATION.md`

**Change 5: Lines 44-46 - Update dataset references**
```markdown
# BEFORE:
Database Files:
- `thermo_cbs_clean.csv` - Primary (CBS calculations)
- Backup: `thermo_original.csv` (deprecated CBS Original)

# AFTER:
Database Files:
- `thermo_cbs_chon.csv` - Primary (29,568 CHON-only molecules, CBS-QB3 level)
- `thermo_pm7.csv` - PM7 optimization results for batch processing
```

---

#### File: `docs/CLAUDE.md`

**Change 6: Lines 293-321 - Add BatchView to CLI structure**

Add after line 301 (after about_view.py):
```markdown
   ├─ batch_view.py            (Batch processing management)
```

Update the table to include:
```markdown
| batch_view.py           ← Batch execution tracking
```

---

#### File: `README.md`

**Verification:** Check lines 70-90 for any references to old dataset names. Based on exploration, README.md is already accurate (references `thermo_cbs_chon.csv` and `thermo_pm7.csv` correctly).

**Action:** No changes needed - documentation is current.

---

### PART 3: Update Planning Document (Optional)

#### File: `docs/plans/magical-napping-moonbeam.md`

**Change 7: Add deprecation notice**

Add at top of file:
```markdown
> **STATUS:** Historical planning document for BATCH 12.
> **Note:** Line numbers may be stale. Refer to current code for actual locations.
```

**Reason:** This prevents future confusion about outdated line references.

---

## 📁 Critical Files to Modify

### Code Changes (4 files)
1. ✅ `src/grimperium/cli/views/databases_view.py` (3 changes: lines 364, 580, 581)
2. ✅ `src/grimperium/cli/views/batch_view.py` (1 change: line 41 + import)

### Documentation Changes (2 files)
3. ✅ `docs/CREST_INTEGRATION.md` (1 change: lines 44-46)
4. ✅ `docs/CLAUDE.md` (1 change: add batch_view.py to structure)
5. 🔍 `docs/plans/magical-napping-moonbeam.md` (optional: add deprecation notice)

---

## ✅ Verification Steps

### 1. Code Quality Checks
```bash
# Type hints (should pass)
mypy src/grimperium/cli/views/databases_view.py --strict

# Linting (should pass)
ruff check src/grimperium/cli/views/databases_view.py
ruff check src/grimperium/cli/views/batch_view.py

# Formatting (should pass)
black --check src/grimperium/cli/views/
```

### 2. File Existence Verification
```bash
# Verify production files exist
ls -lh data/thermo_pm7.csv          # Should exist (PM7 results)
ls -lh data/thermo_cbs_chon.csv     # Should exist (29,568 molecules)

# Verify test file does NOT exist
ls data/test_batch_final.csv        # Should NOT exist (confirms bug)
```

### 3. CLI Functional Test (Manual)

**Test Scenario:** Reproduce the user's error, then verify fix
```bash
# Start CLI
python -m grimperium.cli

# Navigate: DATABASES → Calculate PM7 Values
# Enter: 3 molecules, 20 min CREST timeout, 10 min MOPAC timeout
# Confirm: yes

# EXPECTED BEFORE FIX:
# ✗ Batch CSV not found: data/test_batch_final.csv

# EXPECTED AFTER FIX:
# ✓ Batch processor initializes successfully
# ✓ Loads thermo_pm7.csv
# ✓ Proceeds to molecule selection
```

### 4. Documentation Verification
```bash
# Search for old dataset names in docs (should return 0 matches in active sections)
grep -r "thermo_cbs_clean" docs/*.md
grep -r "thermo_batch_final" docs/*.md
grep -r "test_batch_final" docs/*.md

# After fix, only historical references in DATASETS.md (deprecated section) should remain
```

### 5. Unit Test Coverage
```bash
# Run existing tests (should still pass)
pytest tests/cli/test_views/test_databases_view.py -v

# Check if tests need updates for new paths
# (Tests may need to mock DATA_DIR / "thermo_pm7.csv" instead of old paths)
```

---

## 🚨 Potential Issues & Mitigations

### Issue 1: thermo_pm7.csv Schema Mismatch
**Problem:** `thermo_pm7.csv` currently has 26 columns, but `BatchCSVManager` expects 45 columns (full schema).

**Mitigation:** `BatchCSVManager.load_csv()` validates only required columns (mol_id, smiles, nheavy, status). Optional columns use defaults. **No action needed** - code already handles this gracefully.

**Evidence:** csv_manager.py:150-207 shows optional column handling.

---

### Issue 2: Conformer Details Directory Missing
**Problem:** `data/molecules_pm7/conformer_details/` directory may not exist.

**Mitigation:** Line 599 already creates directory if missing:
```python
detail_dir.mkdir(parents=True, exist_ok=True)
```
**No action needed** - code already handles this.

---

### Issue 3: Test File References in Tests
**Problem:** Unit tests in `tests/cli/test_views/` may reference old paths.

**Mitigation:**
1. Run tests after code changes
2. Update test fixtures if they fail
3. Use `tmp_path` fixtures for test CSV files (best practice)

**Action:** Verify tests pass; update if needed (should be minimal).

---

## 📊 Impact Analysis

### User Impact
- ✅ **IMMEDIATE FIX:** PM7 calculation workflow becomes functional
- ✅ **DOCUMENTATION:** Users see accurate file references in guides
- ✅ **CONSISTENCY:** Code and docs aligned, reducing confusion

### Code Impact
- ✅ **4 line changes** in production code (databases_view.py, batch_view.py)
- ✅ **Low risk:** Using existing constants (`DATA_DIR`) with proven paths
- ✅ **No API changes:** Internal path fixes only, no public interface changes

### Test Impact
- ⚠️ **Potential:** Some tests may need fixture updates
- ✅ **Low risk:** Most tests use mocks, not actual file paths

---

## 🎯 Success Criteria

**Code:**
1. ✅ All hardcoded `Path("data/...")` replaced with `DATA_DIR / ...`
2. ✅ mypy, ruff, black pass with 0 errors
3. ✅ Existing tests pass (or updated tests pass)

**CLI Functionality:**
4. ✅ User can navigate to DATABASES → Calculate PM7 Values
5. ✅ No "Batch CSV not found" error appears
6. ✅ Batch processor initializes successfully with real `thermo_pm7.csv`

**Documentation:**
7. ✅ No references to `thermo_cbs_clean.csv`, `thermo_batch_final.csv`, or `test_batch_final.csv` in active documentation sections
8. ✅ `batch_view.py` listed in CLAUDE.md CLI structure
9. ✅ CREST_INTEGRATION.md references current dataset files

---

## 📝 Implementation Sequence

### Phase 1: Code Fixes (20 min)
1. Fix `databases_view.py` lines 364, 580, 581
2. Fix `batch_view.py` line 41 + add import
3. Run `mypy`, `ruff`, `black` checks
4. Run existing tests

### Phase 2: Documentation Updates (10 min)
5. Update `docs/CREST_INTEGRATION.md` lines 44-46
6. Update `docs/CLAUDE.md` CLI structure section
7. (Optional) Add deprecation notice to `docs/plans/magical-napping-moonbeam.md`

### Phase 3: Verification (10 min)
8. Verify files exist (`ls data/*.csv`)
9. Search for old dataset names in docs (`grep -r`)
10. Manual CLI test (DATABASES → Calculate PM7 Values)

### Phase 4: Commit & Document (5 min)
11. Commit code changes: `fix(cli): replace hardcoded test paths with production CSV paths`
12. Commit docs changes: `docs: update dataset references to current file names`
13. Update CHANGELOG.md if needed

---

## 🔐 Quality Gates

Before considering this task complete:
- [ ] All 4 code files modified successfully
- [ ] mypy --strict passes (0 errors)
- [ ] ruff check passes (0 errors)
- [ ] black --check passes (0 changes needed)
- [ ] pytest tests/cli/test_views/ passes
- [ ] Manual CLI test confirms PM7 calculation workflow works
- [ ] grep confirms no stale dataset names in active docs
- [ ] Code review: all hardcoded paths replaced with constants

---

## 📚 References

**Exploration Agents:**
- Agent ab32d21: CLI bug scope analysis
- Agent a0acbf7: Documentation mismatch analysis

**Key Files:**
- constants.py: `DATA_DIR` definition (line 115)
- csv_manager.py: BatchCSVManager schema (lines 127-148)
- dataset_manager.py: DatasetManager initialization (lines 118-144)

**Datasets:**
- `data/thermo_cbs_chon.csv`: 29,568 CHON molecules (source)
- `data/thermo_pm7.csv`: PM7 optimization results (working CSV)

---

**END OF PLAN**
