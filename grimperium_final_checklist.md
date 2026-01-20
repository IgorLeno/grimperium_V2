# 📋 DOCUMENTO 4: GRIMPERIUM FINAL CHECKLIST

```markdown
# GRIMPERIUM V2 - FINAL IMPLEMENTATION CHECKLIST

**Prepared:** 2026-01-20  
**Target:** Complete Phase A validation fixes  
**Status:** Ready for execution  

---

## PRE-IMPLEMENTATION VERIFICATION

Before starting, verify:

### Repository State
- [ ] Repository synchronized: `git status` shows no uncommitted changes
- [ ] Latest commit visible: `git log --oneline -1`
- [ ] Branch is correct (main or feature branch)
- [ ] Can list files: `git ls-tree HEAD src/grimperium/cli/settings_manager.py`

### Key Files Located
- [ ] `src/grimperium/cli/settings_manager.py` exists
- [ ] `src/grimperium/crest_pm7/molecule_processor.py` exists
- [ ] `src/grimperium/crest_pm7/config.py` exists
- [ ] `docs/CREST_INTEGRATION.md` exists (for template)
- [ ] `docs/MOPAC_INTEGRATION.md` exists (for template)
- [ ] CSV writer file: _________________________ (to be located)

### Test Data Available
- [ ] 3-molecule test batch results exist
- [ ] `thermo_pm7.csv` accessible
- [ ] mol_00001, mol_00002, mol_00003 data visible

---

## TASK 1: SETTINGS SYSTEM OVERHAUL (30 MIN)

### Subtask 1.1: Add New Fields to CRESTSettings

**File:** `src/grimperium/cli/settings_manager.py`  
**Location:** `@dataclass class CRESTSettings`

- [ ] Locate `class CRESTSettings` block
- [ ] Add `CREST_METHOD_OPTIONS = ["gfn2", "gfnff", "gfn2//gfnff"]`
- [ ] Add `QUICK_MODE_OPTIONS = ["off", "quick", "squick", "mquick"]`
- [ ] Add `crest_method: str = "gfn2"` field
- [ ] Add `quick_mode: str = "off"` field
- [ ] REMOVE `gfnff: bool` field
- [ ] REMOVE `quick: bool` field

**Test:**
```python
settings = CRESTSettings()
assert settings.crest_method == "gfn2"
assert settings.quick_mode == "off"
```

### Subtask 1.2: Update `to_dict()` Method

**File:** `src/grimperium/cli/settings_manager.py`  
**Location:** `def to_dict(self) -> dict:`

- [ ] Find line with `"crest_gfnff": str(self.crest.gfnff)` → DELETE
- [ ] Find line with `"crest_quick": str(self.crest.quick)` → DELETE
- [ ] Add `"crest_method": self.crest.crest_method,`
- [ ] Add `"crest_quick_mode": self.crest.quick_mode,`

**Test:**
```python
sm = SettingsManager()
d = sm.to_dict()
assert "crest_method" in d
assert "crest_quick_mode" in d
assert "crest_gfnff" not in d
assert "crest_quick" not in d
```

### Subtask 1.3: Update `from_dict()` Method

**File:** `src/grimperium/cli/settings_manager.py`  
**Location:** `def from_dict(self, data: dict):`

- [ ] Add backward compatibility block BEFORE conversion loop:
```python
# Backward compatibility for old toggles
if "crest_gfnff" in data and "crest_method" not in data:
    data["crest_method"] = "gfnff" if self._parse_bool(data["crest_gfnff"]) else "gfn2"

if "crest_quick" in data and "crest_quick_mode" not in data:
    data["crest_quick_mode"] = "quick" if self._parse_bool(data["crest_quick"]) else "off"
```

- [ ] In conversion map, REMOVE entries for `"crest_gfnff"` and `"crest_quick"`
- [ ] Add entries for `"crest_method"` and `"crest_quick_mode"` (string values, not boolean)

**Test:**
```python
# Old format
sm = SettingsManager()
sm.from_dict({"crest_gfnff": "true", "crest_quick": "true"})
assert sm.crest.crest_method == "gfnff"
assert sm.crest.quick_mode == "quick"

# New format
sm2 = SettingsManager()
sm2.from_dict({"crest_method": "gfn2//gfnff", "crest_quick_mode": "squick"})
assert sm2.crest.crest_method == "gfn2//gfnff"
assert sm2.crest.quick_mode == "squick"
```

### Subtask 1.4: Update HELP_TEXT Dictionary

**File:** `src/grimperium/cli/settings_manager.py`  
**Location:** `HELP_TEXT = {`

- [ ] Find and REMOVE entries:
  - `"crest_gfnff": ...`
  - `"crest_quick": ...`

- [ ] Add new entries:
  - `"crest_method": "Choose CREST quantum method: gfn2 (balanced), gfnff (faster), gfn2//gfnff (two-step)"`
  - `"crest_quick_mode": "Choose speed/accuracy tradeoff: off (full), quick (fast), squick (super-fast), mquick (fastest)"`

### Subtask 1.5: Update show_crest_summary()

**File:** `src/grimperium/cli/settings_manager.py`  
**Location:** `def show_crest_summary(self):`

- [ ] Find line: `table.add_row("GFN-FF Force Field", ...)`
- [ ] REPLACE WITH: `table.add_row("CREST Method", self.crest.crest_method.upper())`
- [ ] Find line: `table.add_row("Quick Mode", ...)`
- [ ] REPLACE WITH: `table.add_row("Quick Mode", self.crest.quick_mode)`

### Subtask 1.6: Replace display_crest_menu()

**File:** `src/grimperium/cli/settings_manager.py`  
**Location:** `def display_crest_menu(self):`

- [ ] Keep all v3, xtb, nci toggles AS-IS
- [ ] REPLACE gfnff toggle with dropdown:
```python
elif choice == "Set CREST Method":
    self.crest.crest_method = questionary.select(
        "Choose CREST method:",
        choices=["gfn2", "gfnff", "gfn2//gfnff"]
    ).ask()
```

- [ ] REPLACE quick toggle with dropdown:
```python
elif choice == "Set Quick Mode":
    self.crest.quick_mode = questionary.select(
        "Choose quick mode:",
        choices=["off", "quick", "squick", "mquick"]
    ).ask()
```

**Test (Manual):**
```bash
# Run CLI, navigate to CREST settings
# Verify:
# - "Set CREST Method" shows with 3 options
# - "Set Quick Mode" shows with 4 options
# - Selections persist after saving
```

---

## TASK 2: LOCATE CSV WRITER (10 MIN)

### Subtask 2.1: Search for CSV Output

**Run:**
```bash
grep -r "\.to_csv\|thermo_pm7" src/grimperium/ --include="*.py" -n
grep -r "DataFrame\|pd\.DataFrame" src/grimperium/crest_pm7/ --include="*.py" -n
grep -r "most_stable_hof\|H298_pm7" src/ --include="*.py" -n
```

### Subtask 2.2: Document Location

- [ ] CSV writer file: `_________________________________`
- [ ] Function name: `_________________________________`
- [ ] Approximate line: `_________________________________`
- [ ] Current columns count: `_________________________________`

**Expected:** Function that converts molecule results to pandas DataFrame, then saves to CSV

---

## TASK 3: CSV SCHEMA UPDATE (45 MIN)

**File:** `[CSV Writer Location - Found in Task 2]`

### Subtask 3.1: Remove Old Columns

- [ ] Remove all code related to `crest_timeout` column
- [ ] Remove all code related to `mopac_timeout` column

**Verify:**
```python
assert "crest_timeout" not in df.columns
assert "mopac_timeout" not in df.columns
```

### Subtask 3.2: Rename Column

- [ ] Find: `"most_stable_hof": ...`
- [ ] Replace with: `"H298_pm7": ...`

**Verify:**
```python
assert "most_stable_hof" not in df.columns
assert "H298_pm7" in df.columns
```

### Subtask 3.3: Add Metrics Columns

Add after H298_pm7:

- [ ] `"abs_diff": abs(molecule.H298_cbs - molecule.H298_pm7),`
- [ ] `"abs_diff_%": (abs(molecule.H298_cbs - molecule.H298_pm7) / abs(molecule.H298_cbs)) * 100,`

**Verification Formula:**
```
H298_cbs = 0.15234
H298_pm7 = -3.68
abs_diff = abs(0.15234 - (-3.68)) = 3.83234
abs_diff_% = (3.83234 / 0.15234) * 100 = 2516.84%
```

### Subtask 3.4: Add Descriptor Columns

Add to row dict:

- [ ] `"nrotbonds": molecule.nrotbonds,`
- [ ] `"tpsa": molecule.tpsa,`
- [ ] `"aromatic_rings": molecule.aromatic_rings,`

### Subtask 3.5: Add CREST Settings Columns

- [ ] `"xtb": molecule.settings.crest.xtb,`
- [ ] `"v3": molecule.settings.crest.v3,`
- [ ] `"qm": molecule.settings.crest.quick_mode,`
- [ ] `"nci": molecule.settings.crest.nci,`
- [ ] `"c_method": molecule.settings.crest.crest_method,`
- [ ] `"energy_window": molecule.settings.crest.energy_window,`
- [ ] `"rmsd_threshold": molecule.settings.crest.rmsd_threshold,`
- [ ] `"threads": molecule.settings.crest.threads,`

### Subtask 3.6: Add MOPAC Settings Columns

- [ ] `"precise_scf": molecule.settings.mopac.precise_scf,`
- [ ] `"scf_threshold": molecule.settings.mopac.scf_threshold,`

### Subtask 3.7: Add Delta Columns

- [ ] `"delta_1": molecule.deltas.get("delta_1") if hasattr(molecule, 'deltas') else None,`
- [ ] `"delta_2": molecule.deltas.get("delta_2") if hasattr(molecule, 'deltas') else None,`
- [ ] `"delta_3": molecule.deltas.get("delta_3") if hasattr(molecule, 'deltas') else None,`
- [ ] `"conformer_selected": molecule.conformer_selected if hasattr(molecule, 'conformer_selected') else None,`

### Subtask 3.8: Add Rerun Counter

- [ ] `"reruns": molecule.reruns if hasattr(molecule, 'reruns') else 0,`

### Subtask 3.9: Fix mopac_time

- [ ] `"mopac_time": molecule.total_mopac_time if hasattr(molecule, 'total_mopac_time') else 0,`

### Subtask 3.10: Set Column Order

After creating DataFrame, reorder columns to exact sequence:

- [ ] Create list of 49 column names (exact order from spec)
- [ ] Reorder: `df = df[column_order]`

**Expected column order:**
```python
column_order = [
    "mol_id", "status", "smiles", "multiplicity", "charge", "nheavy",
    "H298_cbs", "H298_pm7", "abs_diff", "abs_diff_%",
    "batch_id", "timestamp", "reruns",
    "nrotbonds", "tpsa", "aromatic_rings",
    "crest_status", "xtb", "v3", "qm", "nci", "c_method",
    "energy_window", "rmsd_threshold", "threads",
    "crest_conformers_generated", "crest_time", "num_conformers_selected",
    "mopac_status", "precise_scf", "scf_threshold", "mopac_time",
    "delta_1", "delta_2", "delta_3", "conformer_selected",
    "error_message", "batch_order", "batch_failure_policy",
    "assigned_crest_timeout", "assigned_mopac_timeout"
]
df = df[column_order]
```

**Verification:**
```python
df = pd.read_csv("thermo_pm7.csv")
assert len(df.columns) == 49, f"Got {len(df.columns)}, expected 49"
assert list(df.columns) == column_order
assert "crest_timeout" not in df.columns
assert "H298_pm7" in df.columns
```

---

## TASK 4: MOLECULE PROCESSOR FIXES (20 MIN)

**File:** `src/grimperium/crest_pm7/molecule_processor.py`

### Subtask 4.1: Aggregate MOPAC Time

- [ ] Find where `mopac_time` is set (or initialize)
- [ ] Replace with aggregation:

```python
total_mopac_time = 0
for conformer in self.conformers:
    if hasattr(conformer, 'mopac_execution_time') and conformer.mopac_execution_time > 0:
        total_mopac_time += conformer.mopac_execution_time

result["mopac_time"] = total_mopac_time
```

**Verification:**
```python
result = processor.process(molecule)
result_dict = result.to_dict()
assert result_dict["mopac_time"] > 0
```

### Subtask 4.2: Fix Delta Calculations

- [ ] Find function that calculates deltas
- [ ] Replace with correct logic:

```python
def calculate_deltas(self, hofs: List[float], h298_cbs: float):
    # Sort HOFs by energy (lowest first) and take top 3
    sorted_hofs = sorted(hofs)[:3]
    
    deltas = {}
    for i, hof in enumerate(sorted_hofs):
        deltas[f"delta_{i+1}"] = abs(h298_cbs - hof)
    
    # Find conformer with minimum delta
    min_delta = min(deltas.values())
    best_idx = list(deltas.values()).index(min_delta)
    
    return {
        "delta_1": deltas.get("delta_1"),
        "delta_2": deltas.get("delta_2"),
        "delta_3": deltas.get("delta_3"),
        "conformer_selected": best_idx + 1,  # 1-indexed
        "best_hof": sorted_hofs[best_idx]
    }
```

**Verification:**
```python
hofs = [-3.68, -3.52, -3.45, -2.90, -2.80]
h298_cbs = 0.15234

result = processor.calculate_deltas(hofs, h298_cbs)

assert result["delta_1"] == pytest.approx(3.83234)
assert result["delta_2"] == pytest.approx(3.67234)
assert result["delta_3"] == pytest.approx(3.60234)
assert result["conformer_selected"] == 3
assert result["best_hof"] == -3.45
```

---

## TASK 5: CONFIG & STATUS (15 MIN)

**File:** `src/grimperium/crest_pm7/config.py`

### Subtask 5.1: Verify MoleculeStatus Enum

- [ ] Search for `class MoleculeStatus` or `MoleculeStatus = Enum`
- [ ] Verify it has values: PENDING, CURRENT_BATCH, RUNNING, OK, RERUN, SKIP
- [ ] If missing, add:

```python
from enum import Enum

class MoleculeStatus(Enum):
    PENDING = "Pending"
    CURRENT_BATCH = "Current Batch"
    RUNNING = "Running"
    OK = "OK"
    RERUN = "Rerun"
    SKIP = "Skip"
```

### Subtask 5.2: Add Rerun Counter Field

- [ ] Add to molecule dataclass:
```python
reruns: int = 0  # 0-3, max 3 retries
```

### Subtask 5.3: Implement Rerun Logic

In batch processor, when molecule fails:

- [ ] Find failure handling code
- [ ] Add:

```python
if not success:
    molecule.reruns += 1
    
    if molecule.reruns >= 3:
        molecule.status = MoleculeStatus.SKIP
        molecule.error_message = f"Failed after 3 retries"
    else:
        molecule.status = MoleculeStatus.RERUN
        molecule.error_message = f"Retry attempt {molecule.reruns}"
```

---

## TASK 6: CREATE RDKIT DOCUMENTATION (25 MIN)

**File:** `docs/RDKIT_INTEGRATION.md` (NEW)

- [ ] Create new file in `docs/` directory
- [ ] Copy structure from `docs/CREST_INTEGRATION.md`
- [ ] Add these sections:

1. **Overview**
   - [ ] What is RDKit
   - [ ] Why used in Phase A

2. **Key Capabilities**
   - [ ] Descriptor computation
   - [ ] SMILES parsing
   - [ ] Aromaticity analysis

3. **Phase A Integration**
   - [ ] Entry point in pipeline
   - [ ] 3 descriptors computed

4. **Descriptors Detailed**
   - [ ] nrotbonds (definition, range, example)
   - [ ] tpsa (definition, range, example)
   - [ ] aromatic_rings (definition, range, example)

5. **Workflow Diagram**
   - [ ] Where RDKit fits in pipeline

6. **Implementation**
   - [ ] File location: `src/grimperium/rdkit_handler.py`
   - [ ] Key function signature

7. **Performance**
   - [ ] Execution time (<1ms)
   - [ ] Scalability

8. **Limitations & Future**
   - [ ] Current limitations
   - [ ] Phase B plans

9. **References**
   - [ ] RDKit docs link
   - [ ] Scientific papers

10. **Status**
    - [ ] Phase A: Integrated
    - [ ] Phase B: Expansion planned

**Verification:**
```bash
test -f docs/RDKIT_INTEGRATION.md
wc -l docs/RDKIT_INTEGRATION.md  # Should be ~300-400 lines
```

---

## TESTING PHASE: VALIDATION

### TEST 1: Settings Roundtrip

```python
from src.grimperium.cli.settings_manager import SettingsManager

# Create and modify
sm = SettingsManager()
sm.crest.crest_method = "gfnff"
sm.crest.quick_mode = "quick"

# Convert to dict
d = sm.to_dict()
assert d["crest_method"] == "gfnff"
assert d["crest_quick_mode"] == "quick"

# Create new manager and load
sm2 = SettingsManager()
sm2.from_dict(d)
assert sm2.crest.crest_method == "gfnff"
assert sm2.crest.quick_mode == "quick"

print("✓ TEST 1 PASSED: Settings Roundtrip")
```

- [ ] Run test
- [ ] Verify passes

### TEST 2: Backward Compatibility

```python
sm = SettingsManager()

# Old format
old_dict = {"crest_gfnff": "true", "crest_quick": "false"}
sm.from_dict(old_dict)

assert sm.crest.crest_method == "gfnff"
assert sm.crest.quick_mode == "off"

print("✓ TEST 2 PASSED: Backward Compatibility")
```

- [ ] Run test
- [ ] Verify passes

### TEST 3: CSV Schema

```python
import pandas as pd

df = pd.read_csv("data/molecules_pm7/computed/thermo_pm7.csv")

# Check column count
assert len(df.columns) == 49, f"Expected 49, got {len(df.columns)}"

# Check column order
expected = ["mol_id", "status", "smiles", ..., "assigned_mopac_timeout"]
assert list(df.columns) == expected

# Check no old columns
assert "crest_timeout" not in df.columns
assert "mopac_timeout" not in df.columns
assert "most_stable_hof" not in df.columns

# Check new columns exist
assert "abs_diff" in df.columns
assert "delta_1" in df.columns
assert "conformer_selected" in df.columns

print("✓ TEST 3 PASSED: CSV Schema (49 columns)")
```

- [ ] Run test
- [ ] Verify passes

### TEST 4: Molecule Values

```python
import pandas as pd
import pytest

df = pd.read_csv("data/molecules_pm7/computed/thermo_pm7.csv")

# mol_00001
row1 = df[df["mol_id"] == "mol_00001"].iloc
assert row1["nrotbonds"] == 1
assert 10 < row1["tpsa"] < 15
assert row1["aromatic_rings"] == 1
assert row1["H298_cbs"] == pytest.approx(0.15234)
assert row1["abs_diff"] == pytest.approx(3.60234, rel=0.01)
assert row1["conformer_selected"] == 3
assert row1["mopac_time"] > 0

# mol_00002
row2 = df[df["mol_id"] == "mol_00002"].iloc
assert row2["nrotbonds"] == 0
assert row2["aromatic_rings"] == 1
assert row2["H298_cbs"] == pytest.approx(-12.64087)

# mol_00003
row3 = df[df["mol_id"] == "mol_00003"].iloc
assert row3["nrotbonds"] == 0
assert row3["aromatic_rings"] == 2
assert row3["H298_cbs"] == pytest.approx(-0.57999)

print("✓ TEST 4 PASSED: Molecule Values")
```

- [ ] Run test
- [ ] Verify passes

### TEST 5: Settings Persistence

```python
import pandas as pd

df = pd.read_csv("data/molecules_pm7/computed/thermo_pm7.csv")

for idx, row in df.iterrows():
    # All molecules should have settings
    assert pd.notna(row["xtb"])
    assert pd.notna(row["v3"])
    assert pd.notna(row["qm"])
    assert row["qm"] in ["off", "quick", "squick", "mquick"]
    assert pd.notna(row["c_method"])
    assert row["c_method"] in ["gfn2", "gfnff", "gfn2//gfnff"]
    assert pd.notna(row["precise_scf"])

print("✓ TEST 5 PASSED: Settings Persistence")
```

- [ ] Run test
- [ ] Verify passes

### TEST 6: Delta Calculations

```python
import pandas as pd
import pytest

df = pd.read_csv("data/molecules_pm7/computed/thermo_pm7.csv")

for idx, row in df.iterrows():
    if row["num_conformers_selected"] >= 1:
        assert pd.notna(row["delta_1"])
        assert row["conformer_selected"] in [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/923ec5e4-a41c-424c-8ee7-10a942385a05/WORKFLOW.md)
    
    if row["num_conformers_selected"] >= 2:
        assert pd.notna(row["delta_2"])
    
    if row["num_conformers_selected"] >= 3:
        assert pd.notna(row["delta_3"])

print("✓ TEST 6 PASSED: Delta Calculations")
```

- [ ] Run test
- [ ] Verify passes

### TEST 7: CLI Menu (Manual)

- [ ] Run Grimperium CLI
- [ ] Navigate to CREST settings
- [ ] Verify "Set CREST Method" option appears with 3 choices: gfn2, gfnff, gfn2//gfnff
- [ ] Verify "Set Quick Mode" option appears with 4 choices: off, quick, squick, mquick
- [ ] Select gfnff and quick
- [ ] Return to main menu
- [ ] Re-enter CREST settings
- [ ] Verify selections persisted

- [ ] Manual test passed

### TEST 8: Documentation

- [ ] File exists: `docs/RDKIT_INTEGRATION.md`
- [ ] File is readable and has ~300-400 lines
- [ ] Contains all 10 required sections
- [ ] Formatting consistent with other docs
- [ ] No broken links
- [ ] Valid markdown syntax

- [ ] Documentation check passed

---

## FINAL VALIDATION

### Pre-Integration Checklist

**Code Quality:**
- [ ] No syntax errors
- [ ] No undefined variables
- [ ] Proper imports

**Functionality:**
- [ ] Settings roundtrip works
- [ ] Backward compatibility verified
- [ ] CSV has 49 columns in correct order
- [ ] All values populated for test molecules
- [ ] No old column names remain

**Performance:**
- [ ] Settings changes don't impact speed
- [ ] CSV export still completes quickly
- [ ] No memory leaks

**Documentation:**
- [ ] RDKit doc created
- [ ] Help text updated
- [ ] No dead links

### Run Full Test Suite

```bash
pytest tests/test_phase_a_validation.py -v
```

Expected output:
```
test_settings_roundtrip PASSED
test_settings_backward_compat PASSED
test_csv_schema_49_columns PASSED
test_csv_column_order PASSED
test_molecule_00001_values PASSED
test_molecule_00002_values PASSED
test_molecule_00003_values PASSED
test_delta_calculations PASSED
test_mopac_time_aggregation PASSED
test_settings_persistence PASSED
test_rdkit_doc_exists PASSED

======================== 11 passed in 3.42s ========================
```

- [ ] All tests pass

### Run 3-Molecule Batch

```bash
python -m grimperium.cli.main
# Select: Calculate PM7 Values
# Set: 3 molecules, 20 min timeout
# Run batch
```

Expected:
- [ ] Batch completes successfully
- [ ] 3 molecules show status "OK"
- [ ] Output CSV matches schema
- [ ] No errors or warnings

### Compare Output CSV

After batch:
- [ ] Read `thermo_pm7.csv`
- [ ] Verify 49 columns
- [ ] Verify all values present
- [ ] Verify metrics calculated
- [ ] Verify settings tracked
- [ ] No NaN values (except where expected)

---

## COMPLETION CHECKLIST

### Settings System
- [ ] Dropdowns implemented for method and quick-mode
- [ ] Old toggles removed
- [ ] Backward compatibility works
- [ ] to_dict/from_dict symmetric
- [ ] Help text updated

### CSV Schema
- [ ] 49 columns total
- [ ] Exact order verified
- [ ] Old columns removed
- [ ] New columns added
- [ ] Metrics calculated
- [ ] Settings persisted

### Data Validation
- [ ] mopac_time aggregated
- [ ] Delta calculations correct
- [ ] Conformer selection logic works
- [ ] Best conformer identified by minimum delta
- [ ] H298_pm7 uses best conformer HOF

### Documentation
- [ ] RDKIT_INTEGRATION.md created
- [ ] 10 sections completed
- [ ] Formatting consistent
- [ ] No broken references

### Testing
- [ ] All 11 unit tests pass
- [ ] 3-molecule batch successful
- [ ] Output CSV valid
- [ ] Manual CLI test passed

### No Regressions
- [ ] Existing functionality unchanged
- [ ] No new bugs introduced
- [ ] Performance maintained
- [ ] Error handling intact

---

## SIGN-OFF

**Implementation Complete:** [ ]  
**All Tests Pass:** [ ]  
**Ready for Phase B:** [ ]  

**Completed By:** _______________________  
**Date:** _______________________  
**Commit Hash:** _______________________  

---

## IF SOMETHING BREAKS

### "Settings won't save"
→ Check from_dict() is handling new keys  
→ Verify settings.py imports are correct  
→ Check JSON serialization

### "CSV has wrong column order"
→ Verify column_order list matches spec exactly  
→ Check df reorder happens before to_csv()  
→ Count columns: should be 49

### "Deltas are wrong"
→ Verify formula: delta = |H298_cbs - hof|  
→ Check conformer_selected picks minimum delta  
→ Verify 1-indexed (not 0-indexed)

### "mopac_time is zero"
→ Check aggregation loop runs  
→ Verify conformer.mopac_execution_time exists  
→ Check values aren't being reset

### "Old CSV won't load"
→ Verify backward compat code in from_dict()  
→ Check boolean parsing (_parse_bool)  
→ Test with actual old CSV

### "Tests fail"
→ Run individually: `pytest tests/test_X.py -v`  
→ Check test data exists  
→ Verify CSV path correct  
→ Look at exact error message

---

**Status:** ✅ Complete Checklist  
**Ready for:** Implementation  
**Expected Time:** 3-4 hours  
```
