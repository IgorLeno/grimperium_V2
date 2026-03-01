# Plan: Refatorar Pipeline de Seleção PM7-only + Parsing MOPAC

## Context

The current pipeline selects conformers by minimizing `|H298_cbs - HOF|` (CBS-based selection via `DeltaCalculations.calculate_deltas_and_select()`). This must change to **PM7-only selection**: pick the conformer with the lowest `FINAL HEAT OF FORMATION` from MOPAC, regardless of CBS. Additionally, MOPAC inputs need `EF` + `AUX` keywords, and electronic descriptors must be extracted from the selected conformer's `.aux` file.

**Key invariants:**
- `H298_cbs` is already in kcal/mol (confirmed from `data/thermo_cbs_chon.csv`). No conversion needed.
- `.aux` format is MOPAC2016.22.234L with Fortran `D` notation (e.g. `+0.577997D+02`).

---

## Step 1: Add `crest_rank` to `ConformerData`

**File:** `src/grimperium/crest_pm7/molecule_processor.py`

- Add field `crest_rank: int = 0` after `index` (line ~54)
- Update docstring to document it
- Update `to_dict()` to include `"crest_rank": self.crest_rank`
- In `process()` (~line 545), set `crest_rank=idx + 1` when creating `ConformerData`

CREST output is already sorted by energy ascending, so `idx + 1` gives the correct 1-based rank.

---

## Step 2: Add `EF` + `AUX` to MOPAC keywords, remove `1SCF`

**File:** `src/grimperium/crest_pm7/mopac_optimizer.py`

In `_create_mopac_input()` (lines 83-90):

```python
# BEFORE:
keywords = ["PM7"]
...
else:
    keywords.append("1SCF")

# AFTER:
keywords = ["PM7", "EF", "AUX"]
precise_scf = config.mopac_precise_scf if config is not None else True
if precise_scf:
    keywords.append("PRECISE")
    if config is not None:
        keywords.append(f"SCFCRT={format(config.mopac_scf_threshold, '.1e')}")
# No else — 1SCF removed entirely
```

---

## Step 3: Add `k_selected_pm7` + `get_selected_conformer()` to `PM7Result`

**File:** `src/grimperium/crest_pm7/molecule_processor.py`

1. Add field `k_selected_pm7: int | None = None` to `PM7Result` (after `num_conformers_selected`)
2. Add method:
   ```python
   def get_selected_conformer(self) -> ConformerData | None:
       successful = self.successful_conformers
       if not successful:
           return None
       return min(successful, key=lambda c: c.energy_hof)  # type: ignore[arg-type]
   ```
3. Update `to_dict()` to include `"k_selected_pm7": self.k_selected_pm7`
4. In `process()`, after energy difference calc (~line 593):
   ```python
   selected = result.get_selected_conformer()
   if selected is not None:
       result.k_selected_pm7 = selected.crest_rank
   ```

---

## Step 4: Create MOPAC `.aux` descriptor parser

**New file:** `src/grimperium/crest_pm7/mopac_descriptors.py`

Based on **real MOPAC2016.22.234L output** (H2O with `PM7 PRECISE EF AUX`).

### Key format details (from real `.aux`):
- Values use Fortran `D` notation: `+0.577997D+02` → need `parse_fortran_float()` helper
- `DIPOLE:DEBYE=` is scalar magnitude (not vector)
- HOMO/LUMO from `NUM_ELECTRONS` + `EIGENVALUES[n]=` array (no dedicated key)
- COSMO area/volume: `AREA:SQUARE ANGSTROMS=`, `VOLUME:CUBIC ANGSTROMS=`
- `EIGENVALUES` are plain floats (not Fortran D notation)

### Functions:

1. **`parse_fortran_float(s: str) -> float`** — Convert `D` notation to Python float

2. **`_parse_aux_file(aux_file: Path) -> dict[str, Any]`** — Parse `.aux` with these regex patterns:
   | Key in `.aux` | Regex | Descriptor |
   |---|---|---|
   | `DIPOLE:DEBYE=+0.212919D+01` | `r"DIPOLE:DEBYE=([+-]?\d+\.\d+D[+-]\d+)"` | `mopac_dipole_debye` |
   | `IONIZATION_POTENTIAL:EV=+0.121155D+02` | `r"IONIZATION_POTENTIAL:EV=([+-]?\d+\.\d+D[+-]\d+)"` | `mopac_ionization_potential_ev` |
   | `AREA:SQUARE ANGSTROMS=+0.424300D+02` | `r"AREA:SQUARE ANGSTROMS=([+-]?\d+\.\d+D[+-]\d+)"` | `mopac_cosmo_area_a2` |
   | `VOLUME:CUBIC ANGSTROMS=+0.251718D+02` | `r"VOLUME:CUBIC ANGSTROMS=([+-]?\d+\.\d+D[+-]\d+)"` | `mopac_cosmo_volume_a3` |
   | `GRADIENT_NORM:KCAL/MOL/ANGSTROM=+0.883541D-02` | `r"GRADIENT_NORM:KCAL/MOL/ANGSTROM=([+-]?\d+\.\d+D[+-]\d+)"` | `mopac_gradient_norm` |
   | `NUMBER_SCF_CYCLES=7` | `r"NUMBER_SCF_CYCLES=(\d+)"` | `mopac_num_scf_cycles` |
   | `POINT_GROUP=C2v` | `r"POINT_GROUP=(\S+)"` | `mopac_point_group` |
   | `CPU_TIME:SECONDS=0.09` | `r"CPU_TIME:SECONDS=\s*([\d.]+)"` | `mopac_time_s` |
   | `NUM_ELECTRONS=08` + `EIGENVALUES[6]=...` | Combined | `mopac_homo_ev`, `mopac_lumo_ev`, `mopac_gap_ev` |

   **HOMO/LUMO algorithm:**
   ```python
   homo_idx = (num_electrons // 2) - 1  # 0-indexed
   lumo_idx = num_electrons // 2
   gap = lumo - homo
   ```

3. **`extract_mopac_descriptors(mopac_base_path: Path) -> dict[str, Any]`** — Main entry:
   - Build `.aux` path from base path
   - If `.aux` exists: parse it, log INFO
   - If not: return dict of `np.nan` values, log WARNING
   - **No `.out` fallback** (`.aux` is always generated with `AUX` keyword)

### 11 descriptors extracted:
| Column | Type | Source |
|--------|------|--------|
| `mopac_dipole_debye` | float | Fortran D |
| `mopac_ionization_potential_ev` | float | Fortran D |
| `mopac_homo_ev` | float | EIGENVALUES array |
| `mopac_lumo_ev` | float | EIGENVALUES array |
| `mopac_gap_ev` | float | Computed: lumo - homo |
| `mopac_cosmo_area_a2` | float | Fortran D |
| `mopac_cosmo_volume_a3` | float | Fortran D |
| `mopac_gradient_norm` | float | Fortran D |
| `mopac_num_scf_cycles` | int | Plain int |
| `mopac_point_group` | str | Plain string |
| `mopac_time_s` | float | Plain float |

---

## Step 5: Refactor `csv_enhancements.py`

**File:** `src/grimperium/crest_pm7/csv_enhancements.py`

### 5a. Remove `DeltaCalculations.calculate_deltas_and_select()` (lines 106-164)
- Keep `calculate_abs_diff()` and `calculate_abs_diff_pct()` unchanged

### 5b. Change `update_molecule_with_mopac_results()` signature (lines 286-293)
- **Remove:** `mopac_hof_values: list[float]`
- **Add:** `selected_conformer: ConformerData | None`, `k_selected_pm7: int | None`

### 5c. Replace body logic (lines 311-345)
- Remove delta_1/2/3 calculation + conformer_selected
- Add `target_delta_kcalmol = h298_cbs - h298_pm7` (signed, or None)
- Call `extract_mopac_descriptors()` for selected conformer's output path
- Update `updates` dict:

**Remove:** `delta_1`, `delta_2`, `delta_3`, `conformer_selected`

**Add:**
```python
"target_delta_kcalmol": target_delta_kcalmol,
"k_selected_pm7": k_selected_pm7,
"mopac_dipole_debye": descriptors.get("mopac_dipole_debye"),
"mopac_ionization_potential_ev": descriptors.get("mopac_ionization_potential_ev"),
"mopac_homo_ev": descriptors.get("mopac_homo_ev"),
"mopac_lumo_ev": descriptors.get("mopac_lumo_ev"),
"mopac_gap_ev": descriptors.get("mopac_gap_ev"),
"mopac_cosmo_area_a2": descriptors.get("mopac_cosmo_area_a2"),
"mopac_cosmo_volume_a3": descriptors.get("mopac_cosmo_volume_a3"),
"mopac_gradient_norm": descriptors.get("mopac_gradient_norm"),
"mopac_num_scf_cycles": descriptors.get("mopac_num_scf_cycles"),
"mopac_point_group": descriptors.get("mopac_point_group"),
"mopac_time_s": descriptors.get("mopac_time_s"),
```

---

## Step 6: Update `execution_manager.py`

**File:** `src/grimperium/crest_pm7/batch/execution_manager.py`

Replace lines 341-371 in `_process_molecule()`:

```python
# BEFORE: build mopac_hof_values list (lines 342-358)
# AFTER:
selected_conformer = pm7_result.get_selected_conformer()
k_selected_pm7 = pm7_result.k_selected_pm7

success = CSVManagerExtensions.update_molecule_with_mopac_results(
    csv_manager=self.csv_manager,
    mol_id=mol_id,
    h298_cbs=h298_cbs,
    h298_pm7=h298_pm7,
    selected_conformer=selected_conformer,
    k_selected_pm7=k_selected_pm7,
    batch_settings=batch_settings,
)
```

Update log message (line 374): "deltas" → "descriptors".

---

## Step 7: Update CSV schema in `csv_manager.py`

**File:** `src/grimperium/crest_pm7/batch/csv_manager.py`

In `RESULT_COLUMNS`:
- **Remove:** `delta_1`, `delta_2`, `delta_3`, `conformer_selected`
- **Add:** `target_delta_kcalmol`, `k_selected_pm7`, `mopac_dipole_debye`, `mopac_ionization_potential_ev`, `mopac_homo_ev`, `mopac_lumo_ev`, `mopac_gap_ev`, `mopac_cosmo_area_a2`, `mopac_cosmo_volume_a3`, `mopac_gradient_norm`, `mopac_num_scf_cycles`, `mopac_point_group`, `mopac_time_s`

In `pm7result_to_csv_update()`: replace removed keys with new keys initialized to `None`.

---

## Step 8: Tests (11 new + update existing)

**New file:** `tests/unit/test_pm7_pipeline_refactor.py`

| # | Test | Validates |
|---|------|-----------|
| 1 | `test_mopac_input_includes_ef_aux_keywords` | Step 2: EF+AUX in keywords |
| 2 | `test_mopac_input_no_1scf` | Step 2: 1SCF removed |
| 3 | `test_pm7result_k_selected_pm7_uses_crest_rank` | Step 3: k = crest_rank of min HOF |
| 4 | `test_pm7result_get_selected_conformer_returns_lowest_hof` | Step 3: selection logic |
| 5 | `test_aux_file_parsing_water_molecule` | Step 4: real H2O .aux fixture |
| 6 | `test_aux_parsing_fortran_d_notation` | Step 4: Fortran D→float |
| 7 | `test_aux_parsing_homo_lumo_from_eigenvalues` | Step 4: NUM_ELECTRONS + EIGENVALUES |
| 8 | `test_parsing_handles_missing_fields_returns_nan` | Step 4: robustness |
| 9 | `test_csv_update_calculates_target_delta_signed` | Step 5: signed delta |
| 10 | `test_csv_update_includes_k_selected_pm7` | Step 5: k in updates |
| 11 | `test_execution_manager_passes_selected_conformer_and_k` | Step 6: integration |

**Test fixture:** Real H2O `.aux` content from MOPAC2016.22.234L

**Update existing tests:**
- `tests/test_csv_enhancements_batch_settings.py` — update calls to new signature
- `tests/unit/test_crest_pm7_models.py` — add `crest_rank` to expected `to_dict()` output

---

## Implementation Order

```
Step 1 (crest_rank)       ─┐
Step 2 (EF+AUX keywords)   │── independent, can parallelize
Step 4 (descriptor parser)  │
Step 7 (CSV schema)        ─┘
         │
         ▼
Step 3 (k_selected_pm7)   ── depends on Step 1
         │
         ▼
Step 5 (csv_enhancements)  ── depends on Steps 3, 4
         │
         ▼
Step 6 (execution_manager) ── depends on Steps 3, 5, 7
         │
         ▼
Step 8 (tests)             ── run after all code changes
```

---

## Files Modified

| File | Action |
|------|--------|
| `src/grimperium/crest_pm7/molecule_processor.py` | Edit (Steps 1, 3) |
| `src/grimperium/crest_pm7/mopac_optimizer.py` | Edit (Step 2) |
| `src/grimperium/crest_pm7/mopac_descriptors.py` | **Create** (Step 4) |
| `src/grimperium/crest_pm7/csv_enhancements.py` | Edit (Step 5) |
| `src/grimperium/crest_pm7/batch/execution_manager.py` | Edit (Step 6) |
| `src/grimperium/crest_pm7/batch/csv_manager.py` | Edit (Step 7) |
| `tests/unit/test_pm7_pipeline_refactor.py` | **Create** (Step 8) |
| `tests/test_csv_enhancements_batch_settings.py` | Edit (Step 8) |
| `tests/unit/test_crest_pm7_models.py` | Edit (Step 8) |

**Not modified:** `src/grimperium/core/models.py`, `src/grimperium/crest_pm7/config.py`, `energy_extractor.py`

---

## Verification

1. `pytest tests/ -v` — all tests pass
2. `mypy src/grimperium/crest_pm7/ --strict` — zero new errors
3. `ruff check src/` — clean
4. Manual check: `k_selected_pm7` always 1..5 in output CSV
5. Manual check: `target_delta_kcalmol` is signed (can be negative)
6. Manual check: descriptors extracted from selected conformer only
7. Verify `parse_fortran_float("+0.577997D+02")` returns `57.7997`
