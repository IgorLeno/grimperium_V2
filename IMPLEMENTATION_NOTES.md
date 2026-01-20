# Implementation Notes

## CSV Writer Location

- **File:** `src/grimperium/crest_pm7/batch/csv_manager.py`
- **Function:** `CSVManager._build_row()` (approx line 700-800)
- **CSV write call:** Line 226 - `self.df.to_csv(self.csv_path, index=False)`
- **Column definition:** Line 50 - includes "most_stable_hof" in EXPECTED_COLUMNS
- **Current column count:** 43 (needs expansion to 49)

### Key Findings

1. **Primary field to rename:** "most_stable_hof" → "H298_pm7"
   - Defined in: `batch/models.py:226, 335`
   - Written at: `csv_manager.py:746-747`
   - Used in: `molecule_processor.py:187` (property)

2. **Related data model:**
   - `MoleculeResult` (batch/models.py) contains most_stable_hof field
   - `MoleculeProcessor` (molecule_processor.py) calculates this value
   - `DetailManager` (batch/detail_manager.py) aggregates HOF values

3. **Current CSV structure:**
   - Written by: `CSVManager._build_row()` method
   - Columns managed by: `EXPECTED_COLUMNS` constant (line 50)
   - DataFrame operations: Lines 121, 148, 229, 454-455, 786

### Files to Modify (Phase A)

1. `csv_manager.py` - Update column mapping (Task 6-7)
2. `settings_manager.py` - Add dropdown fields (Task 2-5)
3. `molecule_processor.py` - Fix delta calculations (Task 8-9)
4. `config.py` - Verify status enum (Task 10)

### CSV Data vs Code Mismatch (CRITICAL)

**CSV File:** `data/thermo_pm7.csv` has **41 columns** with Phase A schema:
- ✅ Already uses "H298_pm7" (not most_stable_hof)
- ✅ Already has abs_diff, abs_diff_%
- ✅ Already has delta_1, delta_2, delta_3, conformer_selected
- ✅ Already has settings columns (xtb, v3, qm, nci, c_method, etc.)
- ✅ Already has reruns field
- ✅ Has multiplicity, charge (not in code schema)
- ✅ Has H298_cbs (not in code schema)
- ✅ Has mopac_status (not in code schema)

**Code:** `csv_manager.py` has **45 columns** with old schema:
- ❌ Still uses "most_stable_hof"
- ❌ Missing: abs_diff, abs_diff_%
- ❌ Missing: delta_1, delta_2, delta_3, conformer_selected
- ❌ Missing: reruns
- ❌ Uses old config column names (crest_v3, crest_quick, crest_gfnff vs v3, qm, c_method)

**Action Required:**
The code needs to be updated to match the CSV schema, not the other way around.
Phase A is about **synchronizing code with existing CSV structure**.

### Next Steps

- Update csv_manager.py to match CSV schema (41→49 columns)
- Begin Task 2: Add new CRESTSettings fields to match CSV
- Ensure code writes to correct column names
