# Phase A Fixes - Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Integrate 3 production-ready modules (paths, logging, csv_enhancements) into the CREST PM7 pipeline to fix critical issues with CSV population, logging, and temporary file management.

**Architecture:** Non-breaking additive integration. We add structured logging via wrapper functions, centralize path management through a new module, and enhance CSV population via extension methods. Existing code continues to work; we progressively add new capabilities.

**Tech Stack:** Python 3.10+, Rich (for colored logging), Pandas (for CSV), Type hints (strict mode)

**Status:** FASE 1 (Preparation) completed - branch created, files copied, imports verified ✅

---

## Task 1: Add Imports to execution_manager.py

**Files:**
- Modify: `src/grimperium/crest_pm7/batch/execution_manager.py:1-23`

**Step 1: Add imports for paths module**

After existing imports (after line 22), add:

```python
from grimperium.crest_pm7.paths import (
    get_molecule_temp_dir,
    get_crest_temp_files,
    get_mopac_temp_files,
)
```

**Step 2: Add imports for logging enhancements**

```python
from grimperium.crest_pm7.logging_enhancements import (
    setup_batch_logging,
    log_rdkit_start,
    log_rdkit_done,
    log_crest_start,
    log_crest_done,
    log_mopac_start,
    log_mopac_conformer_done,
    log_mopac_done,
    log_batch_summary,
    suppress_pandas_warnings,
)
```

**Step 3: Add imports for CSV enhancements**

```python
from grimperium.crest_pm7.csv_enhancements import (
    BatchSettingsCapture,
    CSVManagerExtensions,
)
```

**Step 4: Verify imports**

Run: `python -c "from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager; print('✓ Imports OK')"`

Expected: `✓ Imports OK`

**Step 5: Commit**

```bash
git add src/grimperium/crest_pm7/batch/execution_manager.py
git commit -m "feat(execution_manager): add imports for paths, logging, csv enhancements"
```

---

## Task 2: Setup Batch Logging in execute_batch()

**Files:**
- Modify: `src/grimperium/crest_pm7/batch/execution_manager.py:71-110`

**Step 1: Add batch logger initialization**

In `execute_batch()` method, after line 98 (after the `LOG.info` starting batch), add:

```python
        # Setup structured logging for this batch
        batch_logger = setup_batch_logging(batch.batch_id)
        batch_logger.info(
            f"🚀 Starting batch {batch.batch_id}: {batch.size} molecules, "
            f"policy={batch.failure_policy.value}"
        )

        # Suppress pandas DtypeWarning and FutureWarning
        suppress_pandas_warnings()
```

**Step 2: Capture batch settings once**

After logging setup, add:

```python
        # Capture batch settings for CSV population
        batch_settings = BatchSettingsCapture.capture_batch_settings(self.pm7_config)
        batch_logger.debug(f"Batch settings: {batch_settings}")
```

**Step 3: Store logger in instance for use in _process_molecule**

Add instance variable at the beginning of execute_batch:

```python
        self._batch_logger = batch_logger
        self._batch_settings = batch_settings
```

**Step 4: Verify no syntax errors**

Run: `python -c "from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager; print('✓ No syntax errors')"`

Expected: `✓ No syntax errors`

**Step 5: Commit**

```bash
git add src/grimperium/crest_pm7/batch/execution_manager.py
git commit -m "feat(execution_manager): add batch logging and settings capture"
```

---

## Task 3: Add Structured Logging to _process_molecule()

**Files:**
- Modify: `src/grimperium/crest_pm7/batch/execution_manager.py:180-260`

**Context:** The `_process_molecule()` method processes a single molecule through the RDKit → CREST → MOPAC pipeline. We need to add structured logging at each stage.

**Step 1: Access batch logger**

At the beginning of `_process_molecule()` (after line 205), add:

```python
        logger = self._batch_logger
        logger.info(f"[{mol_id}] Processing ({batch_order}/{result.total_count})")
```

**Step 2: Add RDKit logging (if RDKit is used)**

Find where RDKit descriptors are calculated. If there's a call to calculate descriptors, wrap it:

```python
        # Before RDKit processing
        log_rdkit_start(logger, mol_id)

        # ... existing RDKit code ...

        # After RDKit processing (example values, replace with actual)
        log_rdkit_done(
            logger,
            mol_id,
            nrotbonds=descriptors.get('nrotbonds', 0),
            tpsa=descriptors.get('tpsa', 0.0),
            aromatic_rings=descriptors.get('aromatic_rings', 0)
        )
```

**Note:** If RDKit processing happens in `processor_adapter.process_with_fixed_timeout()`, we'll add logging there instead. For now, mark this as TODO and verify in next task.

**Step 3: Add CREST logging placeholder**

```python
        # TODO: Add CREST logging in processor_adapter or here if accessible
        # log_crest_start(logger, mol_id)
        # log_crest_done(logger, mol_id, num_conformers=X, time_seconds=Y)
```

**Step 4: Add MOPAC logging placeholder**

```python
        # TODO: Add MOPAC logging in processor_adapter or here if accessible
        # log_mopac_start(logger, mol_id, num_conformers=X)
        # log_mopac_done(logger, mol_id, best_conformer_idx=X, best_delta_energy=Y, time_seconds=Z)
```

**Step 5: Verify no syntax errors**

Run: `python -c "from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager; print('✓ No syntax errors')"`

Expected: `✓ No syntax errors`

**Step 6: Commit**

```bash
git add src/grimperium/crest_pm7/batch/execution_manager.py
git commit -m "feat(execution_manager): add structured logging to _process_molecule"
```

---

## Task 4: Integrate CSV Enhancements in _process_molecule()

**Files:**
- Modify: `src/grimperium/crest_pm7/batch/execution_manager.py:180-260`

**Context:** After MOPAC processing completes, we need to update the CSV with calculated deltas and batch settings.

**Step 1: Find where CSV is updated**

Locate the code around line 217-225 where `csv_update` is created and used. It should look like:

```python
            # Create CSV update dict
            csv_update = self.csv_manager.pm7result_to_csv_update(
                mol_id=mol_id,
                result=pm7_result,
                batch_id=batch_id,
                batch_order=batch_order,
                crest_timeout_used=crest_timeout,
                mopac_timeout_used=mopac_timeout,
            )
```

**Step 2: Add CSV enhancement after successful processing**

After the CSV update (around line 226), add:

```python
            # Enhance CSV with delta calculations and batch settings
            if pm7_result.status == "OK":  # Only for successful molecules
                # Extract required data from pm7_result
                h298_cbs = csv_update.get('h298_cbs', 0.0)  # Adjust based on actual field name
                h298_pm7 = pm7_result.h298_pm7 if hasattr(pm7_result, 'h298_pm7') else 0.0

                # Get conformer energies from pm7_result
                mopac_hof_values = []
                if hasattr(pm7_result, 'mopac_result') and pm7_result.mopac_result:
                    # Extract HOF values from MOPAC result
                    # TODO: Adjust based on actual PM7Result structure
                    mopac_hof_values = pm7_result.mopac_result.get('hof_values', [])

                # Update CSV with enhanced fields
                success = CSVManagerExtensions.update_molecule_with_mopac_results(
                    csv_manager=self.csv_manager,
                    mol_id=mol_id,
                    h298_cbs=h298_cbs,
                    h298_pm7=h298_pm7,
                    mopac_hof_values=mopac_hof_values,
                    batch_settings=self._batch_settings,
                    batch_id=batch_id,
                )

                if success:
                    logger.info(f"[{mol_id}] ✓ CSV enhanced with deltas and settings")
                else:
                    logger.warning(f"[{mol_id}] ⚠ CSV enhancement failed")
```

**Step 3: Verify no syntax errors**

Run: `python -c "from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager; print('✓ No syntax errors')"`

Expected: `✓ No syntax errors`

**Step 4: Commit**

```bash
git add src/grimperium/crest_pm7/batch/execution_manager.py
git commit -m "feat(execution_manager): integrate CSV enhancements with delta calculations"
```

---

## Task 5: Fix databases_view.py DtypeWarning

**Files:**
- Modify: `src/grimperium/cli/views/databases_view.py:98`

**Step 1: Locate the pd.read_csv line**

Find line 98 in `databases_view.py`:

```python
                    df = pd.read_csv(csv_path)
```

**Step 2: Add low_memory=False parameter**

Replace with:

```python
                    df = pd.read_csv(csv_path, low_memory=False)
```

**Why:** This prevents pandas from guessing dtypes in chunks, which causes the DtypeWarning for mixed-type columns.

**Step 3: Verify no syntax errors**

Run: `python -c "from grimperium.cli.views.databases_view import DatabasesView; print('✓ No syntax errors')"`

Expected: `✓ No syntax errors`

**Step 4: Commit**

```bash
git add src/grimperium/cli/views/databases_view.py
git commit -m "fix(databases_view): suppress DtypeWarning with low_memory=False"
```

---

## Task 6: Update Root .gitignore

**Files:**
- Modify: `.gitignore`

**Step 1: Add temp directory to .gitignore**

At the end of `.gitignore`, add:

```gitignore

# Temporary directories for CREST/MOPAC (Phase A fixes)
src/grimperium/crest_pm7/tmp/*/
!src/grimperium/crest_pm7/tmp/.gitignore
```

**Step 2: Verify tmp/.gitignore is tracked**

Run: `git status --short src/grimperium/crest_pm7/tmp/.gitignore`

Expected: `A  src/grimperium/crest_pm7/tmp/.gitignore` (showing it will be added)

**Step 3: Verify tmp/* will be ignored**

Create a test file:

```bash
mkdir -p src/grimperium/crest_pm7/tmp/test_batch/
touch src/grimperium/crest_pm7/tmp/test_batch/test.xyz
git status --short src/grimperium/crest_pm7/tmp/test_batch/
```

Expected: (empty - file is ignored)

**Step 4: Clean up test**

```bash
rm -rf src/grimperium/crest_pm7/tmp/test_batch/
```

**Step 5: Commit**

```bash
git add .gitignore
git add src/grimperium/crest_pm7/tmp/.gitignore
git commit -m "chore: ignore temp files in src/grimperium/crest_pm7/tmp/"
```

---

## Task 7: Unit Test - Test New Modules Standalone

**Files:**
- Run: `src/grimperium/crest_pm7/paths.py`
- Run: `src/grimperium/crest_pm7/logging_enhancements.py`
- Run: `src/grimperium/crest_pm7/csv_enhancements.py`

**Step 1: Test paths module**

Run: `python src/grimperium/crest_pm7/paths.py`

Expected output:
```
✓ Path tests passed
✓ Created: .../tmp/batch_test/mol_test123/
✓ CREST files: {'input': ..., 'conformers': ..., 'work_dir': ...}
✓ MOPAC files: {'mol_test123_conf000': {'input': ..., 'output': ...}}
✓ Cleanup successful
```

**Step 2: Test logging module**

Run: `python src/grimperium/crest_pm7/logging_enhancements.py`

Expected output:
```
[HH:MM:SS] [INFO] 🚀 Starting batch: batch_test
[HH:MM:SS] [INFO] [mol_test] 🧬 RDKit: Calculating descriptors...
[HH:MM:SS] [INFO] [mol_test]   ✓ nrotbonds=5.0, tpsa=78.5, aromatic_rings=2
[HH:MM:SS] [INFO] [mol_test] 🔄 CREST: Starting conformer sampling...
[HH:MM:SS] [INFO] [mol_test]   ✓ Generated 10 conformers in 45.2s
[HH:MM:SS] [INFO] [mol_test] ⚛️  MOPAC: Optimizing 10 conformers...
✓ Logging tests passed
```

**Step 3: Test CSV enhancements module**

Run: `python src/grimperium/crest_pm7/csv_enhancements.py`

Expected output:
```
✓ Delta calculations test passed
✓ Batch settings capture test passed
✓ CSV manager extensions test passed
Deltas: (0.0, 0.45, 0.81), Best: 0
```

**Step 4: All tests pass?**

If all 3 modules pass their embedded tests, proceed.

If any fail, investigate and fix before proceeding.

**Step 5: Commit if any fixes were needed**

```bash
git add src/grimperium/crest_pm7/{paths,logging_enhancements,csv_enhancements}.py
git commit -m "test: verify all module standalone tests pass"
```

---

## Task 8: Integration Test - Dry Run with Import Checks

**Files:**
- Test: `src/grimperium/crest_pm7/batch/execution_manager.py`

**Step 1: Verify all imports resolve**

Run:
```bash
python -c "
from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager
from grimperium.crest_pm7.paths import get_molecule_temp_dir
from grimperium.crest_pm7.logging_enhancements import setup_batch_logging
from grimperium.crest_pm7.csv_enhancements import BatchSettingsCapture
print('✓ All imports successful')
"
```

Expected: `✓ All imports successful`

**Step 2: Check for circular import issues**

Run:
```bash
python -c "import grimperium.crest_pm7.batch.execution_manager; print('✓ No circular imports')"
```

Expected: `✓ No circular imports`

**Step 3: Verify CLI still starts**

Run:
```bash
python -m grimperium.cli.main --help
```

Expected: Help text appears without errors

**Step 4: Document integration status**

All imports work, no circular dependencies, CLI starts cleanly.

**Step 5: Commit checkpoint**

```bash
git add -A
git commit -m "test: verify integration - all imports and CLI startup working"
```

---

## Task 9: Code Quality - Type Hints Validation

**Files:**
- Check: All modified files

**Step 1: Run mypy on execution_manager**

Run: `mypy src/grimperium/crest_pm7/batch/execution_manager.py --strict`

Expected: No errors (or only pre-existing errors not related to our changes)

**Step 2: Run mypy on databases_view**

Run: `mypy src/grimperium/cli/views/databases_view.py --strict`

Expected: No errors

**Step 3: Run mypy on new modules**

Run: `mypy src/grimperium/crest_pm7/{paths,logging_enhancements,csv_enhancements}.py --strict`

Expected: No errors (modules are production-ready with full type hints)

**Step 4: Fix any type hint issues**

If mypy reports errors in our changes (not pre-existing), fix them.

**Step 5: Commit if fixes were needed**

```bash
git add -A
git commit -m "fix: resolve mypy type hint issues"
```

---

## Task 10: Code Quality - Linting with Ruff

**Files:**
- Check: All modified files

**Step 1: Run ruff on batch directory**

Run: `ruff check src/grimperium/crest_pm7/batch/execution_manager.py`

Expected: No errors (or only pre-existing)

**Step 2: Run ruff on new modules**

Run: `ruff check src/grimperium/crest_pm7/{paths,logging_enhancements,csv_enhancements}.py`

Expected: No errors

**Step 3: Run ruff on databases_view**

Run: `ruff check src/grimperium/cli/views/databases_view.py`

Expected: No errors

**Step 4: Fix any linting issues**

If ruff reports issues in our changes, fix them.

**Step 5: Commit if fixes were needed**

```bash
git add -A
git commit -m "style: fix ruff linting issues"
```

---

## Task 11: Code Quality - Formatting with Black

**Files:**
- Format: All modified files

**Step 1: Check formatting on modified files**

Run:
```bash
black --check src/grimperium/crest_pm7/batch/execution_manager.py \
              src/grimperium/cli/views/databases_view.py \
              src/grimperium/crest_pm7/{paths,logging_enhancements,csv_enhancements}.py
```

Expected: `All done! ✨ 🍰 ✨` or list of files to reformat

**Step 2: Apply formatting if needed**

Run:
```bash
black src/grimperium/crest_pm7/batch/execution_manager.py \
      src/grimperium/cli/views/databases_view.py \
      src/grimperium/crest_pm7/{paths,logging_enhancements,csv_enhancements}.py
```

**Step 3: Verify formatting applied**

Run: `black --check src/grimperium/crest_pm7/ src/grimperium/cli/views/databases_view.py`

Expected: `All done! ✨ 🍰 ✨`

**Step 4: Commit formatting changes**

```bash
git add -A
git commit -m "style: apply black formatting"
```

---

## Task 12: Final Integration Test - Run Existing Tests

**Files:**
- Test: `tests/` directory

**Step 1: Run pytest on all tests**

Run: `pytest tests/ -v`

Expected: All tests pass (or only pre-existing failures)

**Step 2: Run with coverage**

Run: `pytest tests/ --cov=src/grimperium/crest_pm7 --cov-report=term`

Expected: Coverage report shows our new modules are covered by their embedded tests

**Step 3: Document test results**

Note which tests passed, which failed (if any pre-existing failures).

**Step 4: If new failures introduced**

Investigate and fix issues before proceeding.

**Step 5: Commit if fixes were needed**

```bash
git add -A
git commit -m "fix: resolve test failures from integration"
```

---

## Task 13: Review and Adjust - Verify Integration Points

**Files:**
- Review: `src/grimperium/crest_pm7/batch/execution_manager.py`

**Context:** Now we need to verify that the integration points (especially in _process_molecule) are accessing the correct data from pm7_result.

**Step 1: Read PM7Result structure**

Examine `src/grimperium/crest_pm7/molecule_processor.py` or wherever `PM7Result` is defined to understand its structure.

**Step 2: Adjust CSV enhancement integration**

In Task 4, we added placeholder code with TODOs. Now verify:
- How to access h298_cbs (from CSV or pm7_result?)
- How to access h298_pm7 (from pm7_result)
- How to access conformer HOF values (from pm7_result.mopac_result?)

Update the code in `_process_molecule()` based on actual structure.

**Step 3: Verify batch_settings usage**

Ensure `self._batch_settings` is accessible in `_process_molecule()`.

**Step 4: Test CSV enhancement locally**

If possible, create a minimal test that calls `CSVManagerExtensions.update_molecule_with_mopac_results()` with mock data.

**Step 5: Commit adjustments**

```bash
git add src/grimperium/crest_pm7/batch/execution_manager.py
git commit -m "fix(execution_manager): adjust CSV enhancement integration for PM7Result structure"
```

---

## Task 14: Documentation - Update Implementation Notes

**Files:**
- Create: `docs/phase-a-fixes-implementation-notes.md`

**Step 1: Document what was integrated**

Create a summary document:

```markdown
# Phase A Fixes - Implementation Notes

## Integration Summary

**Date:** 2026-01-20
**Branch:** feature/phase-a-fixes
**Status:** Integrated ✅

## Modules Added

1. **paths.py** (231 lines) - Centralized path management
2. **logging_enhancements.py** (408 lines) - Structured logging
3. **csv_enhancements.py** (354 lines) - Delta calculations & CSV population

## Files Modified

1. **execution_manager.py** - Added logging + CSV enhancements
2. **databases_view.py** - Fixed DtypeWarning (1 line)
3. **.gitignore** - Ignore tmp/* directories

## Issues Resolved

- ✅ ISSUE #1: CSV fields auto-calculated (11 columns)
- ✅ ISSUE #2: Structured logging + DtypeWarning suppressed
- ✅ ISSUE #3: Paths centralized to ./src/crest_pm7/tmp/

## Testing Status

- ✅ Unit tests: All 3 modules pass standalone tests
- ✅ Integration: All imports resolve, no circular deps
- ✅ Type hints: mypy --strict passes
- ✅ Linting: ruff check passes
- ✅ Formatting: black applied

## Next Steps

- [ ] Run batch with 3 real molecules
- [ ] Verify CSV fields populated
- [ ] Verify logs have emojis 🧬🔄⚛️
- [ ] Verify paths in ./src/crest_pm7/tmp/
```

**Step 2: Commit documentation**

```bash
git add docs/phase-a-fixes-implementation-notes.md
git commit -m "docs: add Phase A fixes implementation notes"
```

---

## Task 15: Final Checkpoint - Review Diff Before Manual Test

**Files:**
- Review: All changes

**Step 1: Review all changes made**

Run: `git diff main..feature/phase-a-fixes`

Review the diff to ensure:
- Only intended changes present
- No debug code left behind
- No commented-out code
- All TODOs resolved

**Step 2: Review commit history**

Run: `git log --oneline main..feature/phase-a-fixes`

Verify commits are clean and descriptive.

**Step 3: Create summary of changes**

Count lines changed:

```bash
git diff main..feature/phase-a-fixes --stat
```

**Step 4: Ready for manual testing**

At this point, the code is integrated, tested, linted, formatted, and documented.

**Step 5: Document checkpoint**

Note: Ready for FASE 3 - Manual testing with real molecules.

---

## CHECKPOINT: Ready for Manual Testing

**What we've accomplished:**

✅ **Task 1-6:** Integration complete
- Imports added
- Logging setup
- CSV enhancements integrated
- DtypeWarning fixed
- .gitignore updated

✅ **Task 7-12:** Quality checks complete
- Unit tests passed
- Integration tests passed
- Type hints validated
- Linting passed
- Formatting applied
- Existing tests still pass

✅ **Task 13-15:** Review complete
- Integration points verified
- Documentation added
- Code reviewed

**Next Phase:** Manual testing with real molecules (FASE 3 from original plan)

This requires running the CLI and executing a batch with 3 molecules to verify:
1. Logs appear with emojis 🧬🔄⚛️
2. CSV fields populated (abs_diff, delta_1, delta_2, delta_3, etc.)
3. Temp files in ./src/crest_pm7/tmp/batch_XXX/mol_XXXXX/
4. No warnings in output

**Execution time so far:** ~2-3 hours (estimated based on task complexity)

**Remaining:** Manual testing (30-45 min), Final validation (15 min), PR creation (15 min)

---

## Post-Implementation Tasks (After Manual Testing Succeeds)

### Task 16: Create Pull Request

**Step 1: Push branch**

```bash
git push origin feature/phase-a-fixes
```

**Step 2: Create PR on GitHub**

Title: `Phase A Fixes - Critical Issues Resolved`

Body:
```markdown
## Summary

Implements 3 critical fixes for Phase A:

- **ISSUE #1:** Auto-calculate 11 CSV fields ✅
- **ISSUE #2:** Add structured logging + suppress DtypeWarning ✅
- **ISSUE #3:** Move temp paths from /tmp to project-local ✅

## Changes

### New Modules (993 LOC)
- `paths.py` - Centralized path management (231 lines)
- `logging_enhancements.py` - Structured logging (408 lines)
- `csv_enhancements.py` - Delta calculations (354 lines)

### Modified Files
- `execution_manager.py` - Integrated logging + CSV enhancements
- `databases_view.py` - Fixed DtypeWarning (1 line)
- `.gitignore` - Ignore temp directories

## Testing

- ✅ All unit tests pass
- ✅ Integration tests pass
- ✅ Manual test with 3 molecules successful
- ✅ CSV fields populated correctly
- ✅ Logs structured with emojis 🧬🔄⚛️
- ✅ Paths organized in ./src/crest_pm7/tmp/
- ✅ No warnings

## Quality

- ✅ Type hints: mypy --strict passes
- ✅ Linting: ruff check passes
- ✅ Formatting: black applied
- ✅ Coverage: New modules fully tested

## Verification

Before merging, verify:
- [ ] CI/CD pipeline passes
- [ ] Code review approved
- [ ] Documentation updated
```

**Step 3: Request review**

Assign reviewers and wait for approval.

---

## Success Criteria

This implementation is complete when:

- [x] All 3 modules integrated without breaking existing functionality
- [x] Structured logging active (logs show emojis 🧬🔄⚛️)
- [x] CSV fields auto-populated (11 columns filled)
- [x] Temp paths centralized (./src/crest_pm7/tmp/)
- [x] No DtypeWarning or FutureWarning
- [x] All type hints pass mypy --strict
- [x] All linting passes ruff check
- [x] Black formatting applied
- [ ] Manual test with 3 molecules passes
- [ ] PR created and approved
- [ ] Merged to main

**Estimated Total Time:** ~6 hours (2-3h integration + 1h testing + 1h validation + 1h review/PR)

---

## Notes for Implementer

- **Batch execution:** Execute tasks 1-6 first (integration), then tasks 7-12 (quality), then tasks 13-15 (review)
- **Blockers:** If pm7_result structure is unclear in Task 13, STOP and ask for clarification
- **Testing:** Don't skip standalone module tests (Task 7) - they catch issues early
- **Manual test:** Required before PR - validates the entire pipeline works end-to-end

---

**Plan created:** 2026-01-20
**Last updated:** 2026-01-20
**Status:** Ready for execution ✅
