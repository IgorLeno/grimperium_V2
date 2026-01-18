# BATCH 12: CLI Critical Fixes - Implementation Plan

**Project:** Grimperium Delta-Learning Framework
**Phase:** C - CLI Interactive Application
**Bugs:** 11 (3 CRITICAL, 4 HIGH, 4 MEDIUM)
**Coverage Goal:** 242 → 248+ tests, 82% → 85%+ coverage
**Estimated Time:** 90-120 minutes

---

## 🎯 Success Criteria

All fixes must pass 6 quality aspects before commit:
1. **Type Hints:** 100% (`mypy --strict src/`)
2. **Tests:** ≥85% coverage with all edge cases
3. **Linting:** Clean (`ruff check src/`)
4. **Formatting:** Aligned (`black --check src/`)
5. **Performance:** Optimized (no redundant operations)
6. **Correctness:** Validated (domain logic + edge cases)

---

## 📋 Verified Bug Locations

### Critical Blockers
- **Bug #1:** CBS Original reference → `src/grimperium/cli/mock_data.py:79-89`
- **Bug #2:** CREST molecule count → `src/grimperium/cli/views/databases_view.py:48-90`
- **Bug #3:** Calc flow validation → `src/grimperium/cli/views/calc_view.py:146-180`

### High Priority
- **Bug #4:** Update button → `src/grimperium/cli/views/databases_view.py:226-455`
- **Bug #5:** Duplicate headers → `src/grimperium/cli/views/settings_view.py:143-165`
- **Bug #6:** No back button → `src/grimperium/cli/views/settings_view.py:67-104`

### Medium Priority
- **Bug #7:** MOPAC alignment → `src/grimperium/cli/settings_manager.py:557-573`
- **Bug #8:** xTB hierarchy → `src/grimperium/cli/views/settings_view.py:87-91`
- **Bug #9:** CREST toggles → `src/grimperium/cli/settings_manager.py:460-476`
- **Bug #10:** MOPAC clarity → `src/grimperium/cli/settings_manager.py:575-588`
- **Bug #11:** Visual feedback → Multiple files

---

## 🔧 Implementation Approach

### Phase 1: Critical Fixes (30 min)

#### Bug #1: Remove CBS Original Reference
**File:** `src/grimperium/cli/mock_data.py`
**Change:** Remove lines 79-89 (CBS Original database entry)
**Impact:** Database menu will no longer show deleted dataset
**Tests:**
- Verify "CBS Original" not in menu options
- Ensure remaining databases still load correctly

#### Bug #2: Dynamic Molecule Count
**File:** `src/grimperium/cli/views/databases_view.py`
**Change:** Replace hardcoded `molecules=0` with dynamic count from JSON
- Load `n_molecules` from `phase_a_results.json`
- Fallback to actual CSV row count if JSON missing
- Handle empty/missing file gracefully

**Tests:**
- Empty JSON → should count CSV rows
- Valid JSON → should use `n_molecules` field
- Missing file → should show 0 with clear message

#### Bug #3: Calculation Flow Validation
**File:** `src/grimperium/cli/views/calc_view.py`
**Change:** Add multi-stage validation before calculation
- Stage 1: Dataset selection with existence check
- Stage 2: Molecule loading with error handling
- Stage 3: SMILES validation using RDKit
- Stage 4: User confirmation before calculation
- Stage 5: Execute with progress feedback

**New Method:** `validate_molecules(molecules: list[dict]) -> list[dict]`
- Check SMILES syntax using RDKit
- Remove duplicates
- Show rejection summary

**Tests:**
- Invalid SMILES rejection
- Duplicate detection
- Empty molecule list handling
- FileNotFoundError on missing dataset

**Quality Gate 1:** Run `pytest tests/cli/ -v --cov=src/grimperium/cli/ --cov-fail-under=85`

---

### Phase 2: High Priority (30 min)

#### Bug #4: Implement Update Button
**File:** `src/grimperium/cli/views/databases_view.py`
**Change:** Make "Refresh Database" button functional
- Verify `data/` directory exists
- Rescan filesystem for available CSVs
- Update molecule counts from actual files
- Show before/after summary

**New Method:** `refresh_databases_from_filesystem() -> int`
- Returns count of databases found
- Displays each database with row count

**Tests:**
- Missing directory handling
- Empty directory handling
- Valid databases found and counted

#### Bug #5: Remove Duplicate Headers
**File:** `src/grimperium/cli/views/settings_view.py`
**Change:** Consolidate redundant headers
- Remove top-level "⚙️ Current Settings" header
- Keep only panel-specific titles
- Ensure single display of each section

**Tests:**
- Count occurrences of each header (should be 1)
- Verify all settings displayed exactly once

#### Bug #6: Add Back Button to Settings
**File:** `src/grimperium/cli/views/settings_view.py`
**Change:** Add "Back to Main Menu" option
- Add to `get_menu_options()` return list
- Handle "back" action in menu loop
- Ensure proper return to main menu

**Tests:**
- Back button visible in menu
- Selecting back exits settings loop
- Loop continues for other actions

**Quality Gate 2:** Run full test suite with coverage check

---

### Phase 3: Medium Priority (30 min)

#### Bug #7: MOPAC Section Alignment
**File:** `src/grimperium/cli/settings_manager.py`
**Change:** Ensure consistent panel structure across all setting sections
- Standardize Panel title/subtitle format
- Align CREST, MOPAC, xTB displays

#### Bug #8: Move xTB Inside CREST
**File:** `src/grimperium/cli/views/settings_view.py`
**Change:** Restructure menu hierarchy
- Remove xTB from top-level menu
- Add xTB toggle inside CREST settings menu
- Update `get_menu_options()` to reflect change

#### Bug #9: Intuitive CREST Toggles
**File:** `src/grimperium/cli/settings_manager.py`
**Change:** Show current state in menu choices
- Format toggles: `"Toggle v3 Algorithm [✓ ON]"` or `"Toggle v3 Algorithm [○ OFF]"`
- Format setters: `"Set Energy Window (current: 6.0)"`

#### Bug #10: Clear MOPAC Settings
**File:** `src/grimperium/cli/settings_manager.py`
**Change:** Same pattern as Bug #9
- Show ON/OFF status for toggles
- Show current values for setters
- Add brief descriptions

#### Bug #11: Universal Visual Feedback
**Files:** Multiple (`settings_manager.py`, `databases_view.py`, `calc_view.py`)
**Change:** Add confirmation messages after actions
- After toggle: `[green]✓ Setting updated: v3 Algorithm ON[/green]`
- After refresh: `[green]✓ Databases refreshed: 2 found[/green]`
- After calculation: `[green]✓ Calculation complete: 3 molecules processed[/green]`

**Pattern:**
- Success: `[green]✓ [action]: [details][/green]`
- Warning: `[yellow]⚠️  [condition]: [details][/yellow]`
- Error: `[red]❌ [error]: [details][/red]`
- Info: `[cyan]ℹ️  [info]: [details][/cyan]`

**Quality Gate 3:** Final validation before commit

---

## ✅ Verification Checklist

### Automated Tests
```bash
# All tests must pass
pytest tests/ -v --cov=src/ --cov-report=term-missing --cov-fail-under=85

# Expected: 248+ tests PASSED, ≥85% coverage
```

### Quality Tools
```bash
# Type hints (must be 0 errors)
mypy --strict src/

# Linting (must be 0 errors)
ruff check src/

# Formatting (must be no changes)
black --check src/
```

### Manual Smoke Test
Start CLI and verify:
- [ ] Database menu loads without "CBS Original"
- [ ] CREST shows correct molecule count (not 0)
- [ ] Calculation flow validates before executing
- [ ] Update button refreshes database list
- [ ] Settings display has no duplicate headers
- [ ] Settings menu has "Back to Main Menu" button
- [ ] MOPAC section aligned with CREST
- [ ] xTB is inside CREST section (not standalone)
- [ ] CREST toggles show current state (✓/○)
- [ ] MOPAC settings show current values
- [ ] All actions show visual feedback (✓/❌/⚠️)

---

## 🎯 Critical Files to Modify

### Primary Changes
1. `src/grimperium/cli/mock_data.py` (Bug #1)
2. `src/grimperium/cli/views/databases_view.py` (Bugs #2, #4)
3. `src/grimperium/cli/views/calc_view.py` (Bug #3)
4. `src/grimperium/cli/views/settings_view.py` (Bugs #5, #6, #8)
5. `src/grimperium/cli/settings_manager.py` (Bugs #7, #9, #10, #11)

### Test Files
1. `tests/cli/test_databases_view.py` (Bugs #1, #2, #4)
2. `tests/cli/test_calc_view.py` (Bug #3)
3. `tests/cli/test_settings_view.py` (Bugs #5, #6, #8)
4. `tests/cli/test_settings_manager.py` (Bugs #7, #9, #10, #11)

---

## 🚀 Implementation Order

1. **Setup:** Create feature branch `fix/batch-12-cli-critical-fixes`
2. **Critical Fixes:** Bugs #1-3 (blockers first)
3. **Quality Gate 1:** Validate critical fixes pass all tests
4. **High Priority:** Bugs #4-6 (UX improvements)
5. **Quality Gate 2:** Validate high priority fixes
6. **Medium Priority:** Bugs #7-11 (polish)
7. **Quality Gate 3:** Full validation (all 11 bugs + 6 aspects)
8. **Manual Test:** Smoke test CLI interactively
9. **Commit:** Atomic commit with detailed message
10. **Push:** Push to remote for review

---

## 📊 Expected Outcomes

### Before BATCH 12
- Tests: 242 passing, 82% coverage
- Bugs: 11 known issues (3 CRITICAL blockers)
- UX: Confusing settings navigation, stale data references

### After BATCH 12
- Tests: 248+ passing, ≥85% coverage
- Bugs: 0 open (all 11 fixed)
- UX: Clean navigation, accurate data, clear visual feedback

### Quality Metrics
- Type hints: 100% (mypy clean)
- Linting: 0 errors (ruff clean)
- Formatting: Consistent (black clean)
- Performance: Optimized (no redundant filesystem/CSV operations)
- Correctness: Validated (all edge cases tested)

---

## 🔗 Reference

**Detailed Specifications:** User provided complete BEFORE/AFTER code for all 11 bugs
**Testing Requirements:** Each bug has 3-6 test cases specified
**Success Criteria:** 6 quality aspects (type hints, tests, linting, formatting, performance, correctness)

**Key Dependencies:**
- `questionary` (interactive prompts)
- `rich` (console styling)
- `pandas` (CSV operations)
- `rdkit` (SMILES validation for Bug #3)
- `pytest` + `pytest-cov` (testing)
- `mypy`, `ruff`, `black` (quality tools)

---

**Estimated Total Time:** 90-120 minutes
**Branch:** `fix/batch-12-cli-critical-fixes`
**Commit Message Template:**
```
fix(batch-12): resolve 11 critical CLI bugs

- Fix #1: Remove CBS Original database reference (deleted file)
- Fix #2: Make CREST molecule count dynamic (not hardcoded)
- Fix #3: Add multi-stage calculation flow validation
- Fix #4: Implement Update button (scan filesystem)
- Fix #5: Remove duplicate section headers
- Fix #6: Add Back button to settings menu
- Fix #7: Align MOPAC section with CREST format
- Fix #8: Move xTB inside CREST hierarchy
- Fix #9: Show current state in CREST toggles (✓/○)
- Fix #10: Show current values in MOPAC settings
- Fix #11: Add universal visual feedback pattern

Quality: 248+ tests, ≥85% coverage, mypy/ruff/black clean
Related: Phase C BATCH 12 complete
```
