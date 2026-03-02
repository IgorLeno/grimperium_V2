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

## Scope and Compatibility

This is a **pure internal implementation change** — no public API is altered.

| Concern | Impact |
|---------|--------|
| Public API | None. `extract_mopac_descriptors()` signature and return keys unchanged. |
| `.aux` file generation | MOPAC will no longer produce `.aux` files. Existing `.aux` files in batch directories are ignored and can be archived or deleted. |
| `.mop` input files | The `AUX` keyword is removed from generated inputs; recalculating a molecule will not regenerate `.aux` files. |
| External tools/scripts | Run `grep -r "\.aux" . --include="*.py" --include="*.sh"` to check for dependencies. The only surviving `.aux` reference in `src/` is the Windows device-name note in `artifact_manager.py`. |
| Dataset reprocessing | Not required. The same 11 descriptor columns are populated; values may differ slightly due to improved HOMO/LUMO extraction precision. |
| `artifact_manager.py` "AUX" note | This is a Windows reserved device-name caveat (`AUX`, `CON`, `PRN`, etc.), not related to MOPAC auxiliary output. It is intentionally left unchanged. |

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
4. HOMO/LUMO → three format variants (see below); GAP computed from HOMO/LUMO
5. GAP → `LUMO - HOMO`; NaN if either is missing or if `LUMO <= HOMO`

### HOMO/LUMO Format Variants

| # | Condition | Example line |
|---|-----------|-------------|
| 1 | Standard closed-shell | `HOMO LUMO ENERGIES (EV) =  -10.096    -0.351` |
| 2 | Open-shell / unrestricted (SOMO) | `HOMO (SOMO) / LUMO ENERGIES (EV) = -10.096 / -0.351` |
| 3 | Multiline (edge case, unusual point groups) | `HOMO ENERGY (EV) = -10.096` on one line, `LUMO ENERGY (EV) = -0.351` on another |

Parsing must try variants in order (1 → 2 → 3) and stop at the first match.

**Error-handling rules for the HOMO/LUMO extraction routine:**
- If the regex fails (no match), do NOT raise an exception; leave `mopac_homo_ev`
  and `mopac_lumo_ev` as `np.nan`.
- If the regex matches but `float()` conversion fails, log a WARNING with the
  matched string and leave the field as `np.nan`.
- GAP computation: `if LUMO is NaN or HOMO is NaN → GAP = NaN (no log needed)`.
  `if LUMO <= HOMO → GAP = NaN and log WARNING("LUMO (%.3f) <= HOMO (%.3f)")`.
- Downstream code must never assume GAP is positive; always check for NaN first.
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

**Parity test (must pass before removing `.aux` code):**

Add `test_aux_out_parsing_parity` that loads both `H2O_AUX_CONTENT` and
`H2O_OUT_CONTENT` fixtures, parses each with their respective parsers
(`_parse_aux_file` for `.aux`, `_parse_out_file` for `.out`), and asserts all
11 descriptor fields match within a tolerance (e.g., `abs=1e-4` for floats;
exact equality for strings and ints). This test must **pass before** the
`_parse_aux_file` function and `.aux` fixtures are deleted from the suite.
Once it passes, the test itself and all `.aux` artifacts can be removed in a
single follow-up commit.

### Step 6: CHANGELOG entry

Add the following under the `### Changed` section in `CHANGELOG.md`:

```markdown
- **MOPAC descriptor extraction migrated from `.aux` to `.out`** — the internal
  parser now reads all 11 descriptors directly from the human-readable `.out`
  summary block instead of the Fortran binary-style `.aux` file.
  - Removed internal helpers: `parse_fortran_float()`, `_parse_aux_file()`, and
    all Fortran D-notation (`1.234D+02`) conversion logic.
  - HOMO/LUMO energies are now read explicitly from the `.out` summary line,
    eliminating the fragile eigenvalue-index guessing required by `.aux`.
  - HOF extraction is unchanged; it was already sourced from `.out` via
    `energy_extractor.py` — this migration aligns the remaining 10 descriptors.
  - **No user-facing API changes**: `extract_mopac_descriptors()` signature,
    return keys, and column names in the output CSV are all unchanged.
  - Improved robustness: no more truncated `EIGENVALUES` array failures for
    molecules with 50+ electrons.
```

## Validation Strategy

Before promoting to production, validate the new `.out` parser against a
representative dataset:

1. **Parity check (≥50 molecules):** Run both the old `.aux` parser and the new
   `.out` parser over the same molecule set. Compare all 11 descriptor fields
   field-by-field: numeric tolerances `<0.01%` relative error for energies (HOF,
   HOMO, LUMO, IP), `<0.1%` for areas/volumes; exact match for point group and
   SCF count. Document any systematic offsets.

2. **Reference values:** Where published MOPAC reference values exist (e.g., from
   the MOPAC documentation or the Chemperium benchmark set), cross-check at
   least 5 molecules.

3. **Edge-case tests:**
   - Failed geometry optimization (non-convergence): expect all numeric fields
     to be `NaN`, point group to be `None`.
   - Unconverged SCF: `.out` may lack the summary block; parser must return
     `_empty_descriptors()` without raising.
   - Unusual point groups (e.g., `C∞v`, `D∞h`, `T`): verify the point-group
     regex captures the full symbol.
   - Open-shell molecules (doublet/triplet): use SOMO format variant.

4. **Pass/fail criteria for promotion:**
   - ≥95% of molecules have all 11 descriptors non-NaN (same threshold as old
     `.aux` parser).
   - No systematic offset > 0.01 eV for HOMO/LUMO energies.
   - All edge-case tests return graceful `NaN`/`None` without exceptions.

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
