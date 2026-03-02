# Plan: Migrate MOPAC Descriptor Parsing from .aux to .out

## Context

HOMO/LUMO/GAP extraction via `.aux` files is fragile (Fortran D notation, EIGENVALUES parsing). The `.out` file already provides these values in human-readable format. The project already parses HOF from `.out` via `energy_extractor.py`, so this migration aligns with existing patterns.

**Goal:** Remove all `.aux` file dependencies and extract all 11 MOPAC descriptors from `.out` only.

## Files to Modify

| File | Change |
|------|--------|
| `src/grimperium/crest_pm7/mopac_optimizer.py` | Remove `"AUX"` from keywords (line 83) |
| `src/grimperium/crest_pm7/mopac_descriptors.py` | Rewrite parser: .aux → .out regex |
| `src/grimperium/crest_pm7/paths.py` | Remove `"auxiliary"` key (line 208) |
| `src/grimperium/crest_pm7/csv_enhancements.py` | Update comments only (code already passes .out path) |
| `tests/unit/test_pm7_pipeline_refactor.py` | Replace .aux fixtures/tests with .out equivalents |

**NOT modified:** `artifact_manager.py` — its "AUX" reference (line 164) is a Windows reserved device name, unrelated to MOPAC.

## Implementation Steps

### Step 1: Remove AUX keyword from MOPAC input
**File:** `mopac_optimizer.py` line 83
- Change `["PM7", "EF", "AUX"]` → `["PM7", "EF"]`

### Step 2: Rewrite mopac_descriptors.py parser
**File:** `mopac_descriptors.py` (full rewrite of parser logic)

Keep unchanged:
- `_DESCRIPTOR_KEYS` list (11 keys)
- `_empty_descriptors()` helper
- `extract_mopac_descriptors()` public signature (already receives .out path via `mopac_base_path`)

Remove:
- `parse_fortran_float()` function
- `_parse_aux_file()` function
- All Fortran D notation regex

Add:
- `_parse_out_file(out_content: str) -> dict[str, Any]` with regex for .out format

Update `extract_mopac_descriptors()`:
- Read `.out` file directly (not `.aux` derived from it)
- Call `_parse_out_file()` instead of `_parse_aux_file()`

**11 descriptors to extract from .out:**
1. HOF → `FINAL HEAT OF FORMATION = X KCAL`
2. Dipole → `DIPOLE = X DEBYE` (from summary line, not component table)
3. IP → `IONIZATION POTENTIAL = X EV`
4. HOMO/LUMO → `HOMO LUMO ENERGIES (EV) = X Y` (3 format variants)
5. GAP → computed from HOMO/LUMO
6. COSMO area → `COSMO AREA = X SQUARE ANGSTROMS`
7. COSMO volume → `COSMO VOLUME = X CUBIC ANGSTROMS`
8. Gradient → `GRADIENT NORM = X`
9. SCF cycles → `SCF CALCULATIONS = X`
10. Point group → `POINT GROUP: X`
11. Time → `WALL-CLOCK TIME = X SECONDS`

### Step 3: Remove "auxiliary" key from paths.py
**File:** `paths.py` line 208
- Remove `"auxiliary": mol_dir / f"mopac_conf_{conformer_idx}.aux"` from dict
- Update docstring (line 194): remove "auxiliary" from key list

### Step 4: Update comments in csv_enhancements.py
**File:** `csv_enhancements.py`
- Line 220: change ".aux file" → ".out file" in comment
- Line 242: same
- Code logic unchanged (already passes `.out` path)

### Step 5: Update tests
**File:** `test_pm7_pipeline_refactor.py`

Remove:
- `H2O_AUX_CONTENT` fixture
- `TestAuxFileParsing` class (all 7 tests)
- Import of `_parse_aux_file`, `parse_fortran_float`
- `test_mopac_input_includes_ef_aux_keywords` (line 68)

Add:
- `H2O_OUT_CONTENT` fixture (realistic .out snippet)
- `TestOutFileParsing` class with tests:
  - `test_out_parsing_extracts_all_fields` — validates all 11 descriptors
  - `test_out_parsing_homo_lumo_variants` — parametrized for 3 HOMO/LUMO formats
  - `test_out_parsing_missing_fields_returns_nan` — incomplete .out → NaN
  - `test_extract_mopac_descriptors_missing_file` — missing .out → empty dict
- Update `test_mopac_input_includes_ef_aux_keywords` → `test_mopac_input_keywords_no_aux`

### Step 6: CHANGELOG entry

## Verification

```bash
# Lint
ruff check src/grimperium/crest_pm7/mopac_descriptors.py src/grimperium/crest_pm7/mopac_optimizer.py src/grimperium/crest_pm7/paths.py

# Type check
mypy --strict src/grimperium/crest_pm7/mopac_descriptors.py

# Tests
pytest tests/unit/test_pm7_pipeline_refactor.py -v

# Full suite
pytest tests/ -v

# Grep for remaining .aux references (should be zero in src/, except artifact_manager Windows reserved)
grep -rn "\.aux" src/grimperium/crest_pm7/ --include="*.py"
```
