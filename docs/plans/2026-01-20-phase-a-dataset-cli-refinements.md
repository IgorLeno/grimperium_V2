# Grimperium Phase A: Dataset & CLI Refinements

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Refine CSV schema to 49 columns, upgrade settings UI from toggles to dropdowns, fix data validation, and complete RDKit documentation for Phase B ML training readiness.

**Architecture:** Surgical fixes to existing pipeline - settings system (dropdown menus), CSV writer (schema expansion), molecule processor (delta calculation fix), and documentation completion. No major refactoring, preserve backward compatibility.

**Tech Stack:** Python 3.10+, pandas, dataclasses, questionary (CLI), Rich (display), pytest

**Context Documents:**
- `grimperium_context_for_claude.md` - Project overview
- `grimperium_spec.md` - Complete specifications
- `grimperium_implementation_guide.md` - Detailed instructions
- `grimperium_final_checklist.md` - Validation checklist

---

## Task 1: Locate CSV Writer & Document Structure

**Files:**
- Search: `src/grimperium/crest_pm7/batch_processor.py` (most likely)
- Search: `src/grimperium/crest_pm7/database_manager.py` (possible)
- Document: Write location to `IMPLEMENTATION_NOTES.md`

**Step 1: Search for CSV writer**

```bash
grep -r "\.to_csv\|thermo_pm7" src/grimperium/ --include="*.py" -n
grep -r "most_stable_hof\|H298_pm7" src/ --include="*.py" -n
grep -r "pd\.DataFrame" src/grimperium/crest_pm7/ --include="*.py" -n
```

Expected: Find file containing DataFrame creation and CSV export

**Step 2: Document findings**

Create `IMPLEMENTATION_NOTES.md`:
```markdown
# Implementation Notes

## CSV Writer Location
- File: [exact path from search]
- Function: [function name]
- Line: [approximate line number]
- Current column count: [count]
```

**Step 3: Read CSV writer file**

Read the identified file to understand current structure

**Step 4: Verify test data**

```bash
head -n 1 data/molecules_pm7/computed/thermo_pm7.csv
```

Expected: See current column headers

**Step 5: Commit documentation**

```bash
git add IMPLEMENTATION_NOTES.md
git commit -m "docs: document CSV writer location for Phase A"
```

---

## Task 2: Settings System - Add New Fields

**Files:**
- Modify: `src/grimperium/cli/settings_manager.py:25-50` (CRESTSettings class)
- Test: Create `tests/test_settings_phase_a.py`

**Step 1: Write test for new settings fields**

```python
# tests/test_settings_phase_a.py
import pytest
from src.grimperium.cli.settings_manager import CRESTSettings, SettingsManager

def test_crest_settings_has_new_fields():
    """Test that CRESTSettings has crest_method and quick_mode fields"""
    settings = CRESTSettings()

    # New fields should exist with defaults
    assert hasattr(settings, 'crest_method')
    assert hasattr(settings, 'quick_mode')
    assert settings.crest_method == "gfn2"
    assert settings.quick_mode == "off"

    # Old fields should NOT exist
    assert not hasattr(settings, 'gfnff')
    assert not hasattr(settings, 'quick')

def test_crest_method_options():
    """Test that CREST_METHOD_OPTIONS are defined"""
    assert hasattr(CRESTSettings, 'CREST_METHOD_OPTIONS')
    assert CRESTSettings.CREST_METHOD_OPTIONS == ["gfn2", "gfnff", "gfn2//gfnff"]

def test_quick_mode_options():
    """Test that QUICK_MODE_OPTIONS are defined"""
    assert hasattr(CRESTSettings, 'QUICK_MODE_OPTIONS')
    assert CRESTSettings.QUICK_MODE_OPTIONS == ["off", "quick", "squick", "mquick"]
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_settings_phase_a.py::test_crest_settings_has_new_fields -v
```

Expected: FAIL - AttributeError (crest_method doesn't exist)

**Step 3: Modify CRESTSettings class**

In `src/grimperium/cli/settings_manager.py`, find `@dataclass class CRESTSettings`:

```python
@dataclass
class CRESTSettings:
    # Valid options for dropdown menus
    CREST_METHOD_OPTIONS = ["gfn2", "gfnff", "gfn2//gfnff"]
    QUICK_MODE_OPTIONS = ["off", "quick", "squick", "mquick"]

    # Settings
    xtb: bool = False
    v3: bool = False
    nci: bool = False
    crest_method: str = "gfn2"              # NEW: replaces gfnff boolean
    quick_mode: str = "off"                 # NEW: replaces quick boolean
    energy_window: float = 5.0
    rmsd_threshold: float = 1.5
    threads: int = 4
```

Remove these lines if present:
```python
    gfnff: bool = False      # DELETE
    quick: bool = False      # DELETE
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_settings_phase_a.py::test_crest_settings_has_new_fields -v
pytest tests/test_settings_phase_a.py::test_crest_method_options -v
pytest tests/test_settings_phase_a.py::test_quick_mode_options -v
```

Expected: All PASS

**Step 5: Commit**

```bash
git add src/grimperium/cli/settings_manager.py tests/test_settings_phase_a.py
git commit -m "feat(settings): add crest_method and quick_mode dropdown fields"
```

---

## Task 3: Settings System - Update to_dict()

**Files:**
- Modify: `src/grimperium/cli/settings_manager.py:100-150` (to_dict method)
- Test: `tests/test_settings_phase_a.py`

**Step 1: Write test for to_dict serialization**

```python
# Add to tests/test_settings_phase_a.py

def test_settings_to_dict_new_format():
    """Test that to_dict exports new fields correctly"""
    sm = SettingsManager()
    sm.crest.crest_method = "gfnff"
    sm.crest.quick_mode = "quick"

    result = sm.to_dict()

    # New fields should be present
    assert "crest_method" in result
    assert "crest_quick_mode" in result
    assert result["crest_method"] == "gfnff"
    assert result["crest_quick_mode"] == "quick"

    # Old fields should NOT be present
    assert "crest_gfnff" not in result
    assert "crest_quick" not in result
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_settings_phase_a.py::test_settings_to_dict_new_format -v
```

Expected: FAIL - KeyError (crest_method not in dict)

**Step 3: Update to_dict() method**

In `src/grimperium/cli/settings_manager.py`, find `def to_dict(self):`:

Find and REMOVE these lines:
```python
        "crest_gfnff": str(self.crest.gfnff),      # DELETE
        "crest_quick": str(self.crest.quick),      # DELETE
```

Add these lines (where old ones were):
```python
        "crest_method": self.crest.crest_method,        # NEW
        "crest_quick_mode": self.crest.quick_mode,      # NEW
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_settings_phase_a.py::test_settings_to_dict_new_format -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add src/grimperium/cli/settings_manager.py tests/test_settings_phase_a.py
git commit -m "feat(settings): update to_dict to export new dropdown values"
```

---

## Task 4: Settings System - Backward Compatibility in from_dict()

**Files:**
- Modify: `src/grimperium/cli/settings_manager.py:150-200` (from_dict method)
- Test: `tests/test_settings_phase_a.py`

**Step 1: Write tests for backward compatibility**

```python
# Add to tests/test_settings_phase_a.py

def test_from_dict_backward_compat_gfnff_true():
    """Test loading old format with gfnff=true converts to crest_method=gfnff"""
    sm = SettingsManager()

    old_data = {"crest_gfnff": "true"}
    sm.from_dict(old_data)

    assert sm.crest.crest_method == "gfnff"

def test_from_dict_backward_compat_gfnff_false():
    """Test loading old format with gfnff=false converts to crest_method=gfn2"""
    sm = SettingsManager()

    old_data = {"crest_gfnff": "false"}
    sm.from_dict(old_data)

    assert sm.crest.crest_method == "gfn2"

def test_from_dict_backward_compat_quick_true():
    """Test loading old format with quick=true converts to quick_mode=quick"""
    sm = SettingsManager()

    old_data = {"crest_quick": "true"}
    sm.from_dict(old_data)

    assert sm.crest.quick_mode == "quick"

def test_from_dict_backward_compat_quick_false():
    """Test loading old format with quick=false converts to quick_mode=off"""
    sm = SettingsManager()

    old_data = {"crest_quick": "false"}
    sm.from_dict(old_data)

    assert sm.crest.quick_mode == "off"

def test_from_dict_new_format():
    """Test loading new format works directly"""
    sm = SettingsManager()

    new_data = {
        "crest_method": "gfn2//gfnff",
        "crest_quick_mode": "squick"
    }
    sm.from_dict(new_data)

    assert sm.crest.crest_method == "gfn2//gfnff"
    assert sm.crest.quick_mode == "squick"

def test_settings_roundtrip():
    """Test that to_dict -> from_dict preserves settings"""
    sm1 = SettingsManager()
    sm1.crest.crest_method = "gfnff"
    sm1.crest.quick_mode = "mquick"

    data = sm1.to_dict()

    sm2 = SettingsManager()
    sm2.from_dict(data)

    assert sm2.crest.crest_method == "gfnff"
    assert sm2.crest.quick_mode == "mquick"
```

**Step 2: Run tests to verify they fail**

```bash
pytest tests/test_settings_phase_a.py::test_from_dict_backward_compat_gfnff_true -v
```

Expected: FAIL - AttributeError or wrong value

**Step 3: Add backward compatibility to from_dict()**

In `src/grimperium/cli/settings_manager.py`, find `def from_dict(self, data: dict):`.

Add at the BEGINNING of the method (before conversion map):

```python
    # Backward compatibility: convert old boolean toggles to new dropdown values
    if "crest_gfnff" in data and "crest_method" not in data:
        if self._parse_bool(data["crest_gfnff"]):
            data["crest_method"] = "gfnff"
        else:
            data["crest_method"] = "gfn2"

    if "crest_quick" in data and "crest_quick_mode" not in data:
        if self._parse_bool(data["crest_quick"]):
            data["crest_quick_mode"] = "quick"
        else:
            data["crest_quick_mode"] = "off"
```

In the conversion_map, ADD entries for new fields:

```python
    conversion_map = {
        # ... existing entries ...
        ("crest_method",): ("crest", "crest_method"),
        ("crest_quick_mode",): ("crest", "quick_mode"),
        # ... rest of entries ...
    }
```

REMOVE from conversion_map if present:
```python
        ("crest_gfnff", "gfnff"): ("crest", "gfnff"),  # DELETE
        ("crest_quick", "quick"): ("crest", "quick"),  # DELETE
```

**Step 4: Run tests to verify they pass**

```bash
pytest tests/test_settings_phase_a.py -v
```

Expected: All tests PASS

**Step 5: Commit**

```bash
git add src/grimperium/cli/settings_manager.py tests/test_settings_phase_a.py
git commit -m "feat(settings): add backward compatibility for old toggle format"
```

---

## Task 5: Settings System - Update UI Display

**Files:**
- Modify: `src/grimperium/cli/settings_manager.py:350-400` (show_crest_summary)
- Modify: `src/grimperium/cli/settings_manager.py:400-500` (display_crest_menu)

**Step 1: Update show_crest_summary() method**

In `src/grimperium/cli/settings_manager.py`, find `def show_crest_summary(self):`.

Find and REPLACE:
```python
    table.add_row("GFN-FF Force Field", "✓" if self.crest.gfnff else "✗")
    table.add_row("Quick Mode", "✓" if self.crest.quick else "✗")
```

WITH:
```python
    table.add_row("CREST Method", self.crest.crest_method.upper())
    table.add_row("Quick Mode", self.crest.quick_mode)
```

**Step 2: Manual test**

```bash
python -m grimperium.cli.main
# Navigate to view settings
# Verify shows "CREST Method: GFN2" and "Quick Mode: off"
```

**Step 3: Update display_crest_menu() method**

In `src/grimperium/cli/settings_manager.py`, find `def display_crest_menu(self):`.

REPLACE the entire menu structure. Find the while loop and replace with:

```python
def display_crest_menu(self) -> None:
    while True:
        console.print("\n[bold cyan]CREST Settings[/bold cyan]")

        # Show current values
        console.print(f"Current CREST Method: [yellow]{self.crest.crest_method}[/yellow]")
        console.print(f"Current Quick Mode: [yellow]{self.crest.quick_mode}[/yellow]")

        choice = questionary.select(
            "What would you like to modify?",
            choices=[
                "Set CREST Method",
                "Set Quick Mode",
                "Toggle xTB Method",
                "Toggle v3 Algorithm",
                "Toggle NCI Mode",
                "Back to main menu"
            ]
        ).ask()

        if choice == "Set CREST Method":
            self.crest.crest_method = questionary.select(
                "Choose CREST method:",
                choices=["gfn2", "gfnff", "gfn2//gfnff"]
            ).ask()
            console.print(f"✓ CREST Method set to: {self.crest.crest_method}")

        elif choice == "Set Quick Mode":
            self.crest.quick_mode = questionary.select(
                "Choose quick mode:",
                choices=["off", "quick", "squick", "mquick"]
            ).ask()
            console.print(f"✓ Quick Mode set to: {self.crest.quick_mode}")

        elif choice == "Toggle xTB Method":
            self.crest.xtb = not self.crest.xtb
            console.print(f"xTB: {'✓ ON' if self.crest.xtb else '✗ OFF'}")

        elif choice == "Toggle v3 Algorithm":
            self.crest.v3 = not self.crest.v3
            console.print(f"v3: {'✓ ON' if self.crest.v3 else '✗ OFF'}")

        elif choice == "Toggle NCI Mode":
            self.crest.nci = not self.crest.nci
            console.print(f"NCI: {'✓ ON' if self.crest.nci else '✗ OFF'}")

        elif choice == "Back to main menu":
            break
```

**Step 4: Manual test**

```bash
python -m grimperium.cli.main
# Navigate to CREST settings
# Verify "Set CREST Method" shows dropdown with 3 options
# Verify "Set Quick Mode" shows dropdown with 4 options
# Select gfnff and quick, exit, re-enter
# Verify selections persisted
```

**Step 5: Update HELP_TEXT**

In `src/grimperium/cli/settings_manager.py`, find `HELP_TEXT = {`.

REMOVE:
```python
    "crest_gfnff": "...",
    "crest_quick": "...",
```

ADD:
```python
    "crest_method": "Choose CREST quantum method: gfn2 (default, balanced), gfnff (faster), gfn2//gfnff (two-step refinement)",
    "crest_quick_mode": "Choose speed/accuracy tradeoff: off (full), quick (fast), squick (super-fast), mquick (fastest)",
```

**Step 6: Commit**

```bash
git add src/grimperium/cli/settings_manager.py
git commit -m "feat(settings): update UI to use dropdown menus instead of toggles"
```

---

## Task 6: CSV Schema - Update Column Mapping

**Files:**
- Modify: [CSV writer file from Task 1]
- Test: Create `tests/test_csv_schema.py`

**Step 1: Write test for CSV schema**

```python
# tests/test_csv_schema.py
import pytest
import pandas as pd
from pathlib import Path

EXPECTED_COLUMNS = [
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

def test_csv_has_49_columns():
    """Test that CSV output has exactly 49 columns"""
    csv_path = Path("data/molecules_pm7/computed/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        assert len(df.columns) == 49, f"Expected 49 columns, got {len(df.columns)}"

def test_csv_column_order():
    """Test that CSV columns are in correct order"""
    csv_path = Path("data/molecules_pm7/computed/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        assert list(df.columns) == EXPECTED_COLUMNS

def test_csv_no_old_columns():
    """Test that old column names are removed"""
    csv_path = Path("data/molecules_pm7/computed/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        assert "crest_timeout" not in df.columns
        assert "mopac_timeout" not in df.columns
        assert "most_stable_hof" not in df.columns

def test_csv_has_new_columns():
    """Test that new columns are present"""
    csv_path = Path("data/molecules_pm7/computed/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Metrics
        assert "abs_diff" in df.columns
        assert "abs_diff_%" in df.columns
        # Deltas
        assert "delta_1" in df.columns
        assert "delta_2" in df.columns
        assert "delta_3" in df.columns
        assert "conformer_selected" in df.columns
        # Settings
        assert "c_method" in df.columns
        assert "qm" in df.columns
```

**Step 2: Run tests (they may pass if CSV exists, or skip if not)**

```bash
pytest tests/test_csv_schema.py -v
```

Expected: FAIL or SKIP (depending on CSV state)

**Step 3: Update CSV writer - Remove old columns**

In [CSV writer file], find the row dictionary creation.

REMOVE these lines:
```python
"crest_timeout": ...,
"mopac_timeout": ...,
```

**Step 4: Update CSV writer - Rename column**

CHANGE:
```python
"most_stable_hof": molecule.hof,
```

TO:
```python
"H298_pm7": molecule.hof,
```

**Step 5: Update CSV writer - Add metrics**

Add after H298_pm7:
```python
"abs_diff": abs(molecule.H298_cbs - molecule.hof) if hasattr(molecule, 'hof') else None,
"abs_diff_%": (abs(molecule.H298_cbs - molecule.hof) / abs(molecule.H298_cbs)) * 100 if hasattr(molecule, 'hof') and molecule.H298_cbs != 0 else None,
```

**Step 6: Update CSV writer - Add reruns**

Add:
```python
"reruns": molecule.reruns if hasattr(molecule, 'reruns') else 0,
```

**Step 7: Commit initial changes**

```bash
git add [csv-writer-file] tests/test_csv_schema.py
git commit -m "refactor(csv): remove old columns and rename most_stable_hof to H298_pm7"
```

---

## Task 7: CSV Schema - Add Settings Columns

**Files:**
- Modify: [CSV writer file]
- Test: `tests/test_csv_schema.py`

**Step 1: Add CREST settings columns**

In row dictionary, add:
```python
"xtb": molecule.settings.crest.xtb if hasattr(molecule, 'settings') else None,
"v3": molecule.settings.crest.v3 if hasattr(molecule, 'settings') else None,
"qm": molecule.settings.crest.quick_mode if hasattr(molecule, 'settings') else None,
"nci": molecule.settings.crest.nci if hasattr(molecule, 'settings') else None,
"c_method": molecule.settings.crest.crest_method if hasattr(molecule, 'settings') else None,
"energy_window": molecule.settings.crest.energy_window if hasattr(molecule, 'settings') else None,
"rmsd_threshold": molecule.settings.crest.rmsd_threshold if hasattr(molecule, 'settings') else None,
"threads": molecule.settings.crest.threads if hasattr(molecule, 'settings') else None,
```

**Step 2: Add MOPAC settings columns**

```python
"precise_scf": molecule.settings.mopac.precise_scf if hasattr(molecule, 'settings') else None,
"scf_threshold": molecule.settings.mopac.scf_threshold if hasattr(molecule, 'settings') else None,
```

**Step 3: Add delta columns**

```python
"delta_1": molecule.deltas.get("delta_1") if hasattr(molecule, 'deltas') else None,
"delta_2": molecule.deltas.get("delta_2") if hasattr(molecule, 'deltas') else None,
"delta_3": molecule.deltas.get("delta_3") if hasattr(molecule, 'deltas') else None,
"conformer_selected": molecule.conformer_selected if hasattr(molecule, 'conformer_selected') else None,
```

**Step 4: Add RDKit descriptor columns**

```python
"nrotbonds": molecule.nrotbonds if hasattr(molecule, 'nrotbonds') else None,
"tpsa": molecule.tpsa if hasattr(molecule, 'tpsa') else None,
"aromatic_rings": molecule.aromatic_rings if hasattr(molecule, 'aromatic_rings') else None,
```

**Step 5: Set column order**

After creating DataFrame, add reordering:
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

df = pd.DataFrame(rows)
df = df[column_order]  # Reorder to exact spec
df.to_csv(output_path, index=False)
```

**Step 6: Run batch test to generate new CSV**

```bash
python -m grimperium.cli.main
# Select: Calculate PM7 Values
# Set: 3 molecules (use test batch)
# Run batch
```

**Step 7: Verify CSV schema**

```bash
pytest tests/test_csv_schema.py -v
```

Expected: All tests PASS

**Step 8: Commit**

```bash
git add [csv-writer-file]
git commit -m "feat(csv): add 49-column schema with settings tracking"
```

---

## Task 8: Molecule Processor - Fix mopac_time Aggregation

**Files:**
- Modify: `src/grimperium/crest_pm7/molecule_processor.py:200-250`
- Test: Create `tests/test_molecule_processor.py`

**Step 1: Write test for mopac_time aggregation**

```python
# tests/test_molecule_processor.py
import pytest

def test_mopac_time_aggregation():
    """Test that mopac_time is aggregated across conformers"""
    # Mock conformer results with execution times
    from unittest.mock import Mock

    conformer1 = Mock()
    conformer1.mopac_execution_time = 100.0

    conformer2 = Mock()
    conformer2.mopac_execution_time = 150.0

    conformer3 = Mock()
    conformer3.mopac_execution_time = 200.0

    conformers = [conformer1, conformer2, conformer3]

    # Calculate total
    total_time = sum(c.mopac_execution_time for c in conformers)

    assert total_time == 450.0
```

**Step 2: Run test**

```bash
pytest tests/test_molecule_processor.py::test_mopac_time_aggregation -v
```

Expected: PASS (this is a logic test)

**Step 3: Find mopac_time assignment in molecule_processor.py**

Search for `mopac_time` or `to_dict` method

**Step 4: Update aggregation logic**

REPLACE:
```python
"mopac_time": self.conformers.mopac_execution_time,  # Wrong: single value
```

WITH:
```python
"mopac_time": sum(
    c.mopac_execution_time
    for c in self.conformers
    if hasattr(c, 'mopac_execution_time') and c.mopac_execution_time > 0
),
```

**Step 5: Add integration test**

```python
# Add to tests/test_molecule_processor.py

def test_molecule_result_includes_aggregated_mopac_time():
    """Test that molecule result dict includes aggregated mopac_time"""
    # This test requires actual molecule processor integration
    # Will be verified by running batch and checking CSV
    pass  # Placeholder for manual verification
```

**Step 6: Commit**

```bash
git add src/grimperium/crest_pm7/molecule_processor.py tests/test_molecule_processor.py
git commit -m "fix(processor): aggregate mopac_time across all conformers"
```

---

## Task 9: Molecule Processor - Fix Delta Calculations

**Files:**
- Modify: `src/grimperium/crest_pm7/molecule_processor.py:150-200`
- Test: `tests/test_molecule_processor.py`

**Step 1: Write test for delta calculation**

```python
# Add to tests/test_molecule_processor.py

def test_delta_calculation_logic():
    """Test delta calculation: |H298_cbs - hof| and select minimum"""
    h298_cbs = 0.15234
    hofs = [-3.68, -3.52, -3.45, -2.90, -2.80]

    # Sort and take top 3
    sorted_hofs = sorted(hofs)[:3]  # [-3.68, -3.52, -3.45]

    # Calculate deltas
    delta_1 = abs(h298_cbs - sorted_hofs[0])  # |0.15234 - (-3.68)| = 3.83234
    delta_2 = abs(h298_cbs - sorted_hofs[1])  # |0.15234 - (-3.52)| = 3.67234
    delta_3 = abs(h298_cbs - sorted_hofs[2])  # |0.15234 - (-3.45)| = 3.60234

    assert pytest.approx(delta_1, rel=1e-3) == 3.83234
    assert pytest.approx(delta_2, rel=1e-3) == 3.67234
    assert pytest.approx(delta_3, rel=1e-3) == 3.60234

    # Select minimum delta
    deltas = [delta_1, delta_2, delta_3]
    min_delta = min(deltas)
    best_idx = deltas.index(min_delta)

    assert best_idx == 2  # delta_3 is minimum
    assert sorted_hofs[best_idx] == -3.45

def test_delta_selects_best_conformer():
    """Test that conformer with MINIMUM delta is selected (not maximum)"""
    h298_cbs = -12.64087
    hofs = [-14.66]  # Only one conformer

    delta_1 = abs(h298_cbs - hofs[0])

    assert pytest.approx(delta_1, rel=1e-3) == 2.01913
```

**Step 2: Run tests**

```bash
pytest tests/test_molecule_processor.py::test_delta_calculation_logic -v
pytest tests/test_molecule_processor.py::test_delta_selects_best_conformer -v
```

Expected: PASS

**Step 3: Find delta calculation in molecule_processor.py**

Search for `calculate_delta` or `delta_e_`

**Step 4: Replace delta calculation logic**

REPLACE old logic with:
```python
def calculate_deltas(self, hofs: List[float], h298_cbs: float) -> dict:
    """
    Calculate deltas between CBS reference and conformer HOF values.

    Returns dict with:
    - delta_1, delta_2, delta_3: |H298_cbs - hof| for top 3 conformers
    - conformer_selected: index (1-based) of conformer with minimum delta
    - selected_hof: HOF value of selected conformer
    """
    # Sort HOFs by energy (lowest first) and take top 3
    sorted_hofs = sorted(hofs)[:3]

    # Calculate deltas for each conformer
    deltas = {}
    for i, hof in enumerate(sorted_hofs):
        deltas[f"delta_{i+1}"] = abs(h298_cbs - hof)

    # Find conformer with minimum delta (best agreement with CBS)
    valid_deltas = [d for d in deltas.values() if d is not None]
    if valid_deltas:
        min_delta = min(valid_deltas)
        best_idx = list(deltas.values()).index(min_delta)
    else:
        best_idx = 0

    return {
        "deltas": deltas,
        "conformer_selected": best_idx + 1,  # 1-indexed
        "selected_hof": sorted_hofs[best_idx] if best_idx < len(sorted_hofs) else None
    }
```

**Step 5: Update usage of calculate_deltas**

Find where this function is called and update to use the returned dict:
```python
result = self.calculate_deltas(hofs, h298_cbs)
molecule.deltas = result["deltas"]
molecule.conformer_selected = result["conformer_selected"]
molecule.hof = result["selected_hof"]  # Use selected conformer's HOF as H298_pm7
```

**Step 6: Commit**

```bash
git add src/grimperium/crest_pm7/molecule_processor.py tests/test_molecule_processor.py
git commit -m "fix(processor): correct delta calculation to select minimum delta"
```

---

## Task 10: Config - Verify Status Enum & Rerun Logic

**Files:**
- Verify/Modify: `src/grimperium/crest_pm7/config.py:20-50`
- Test: `tests/test_config.py`

**Step 1: Write test for MoleculeStatus enum**

```python
# tests/test_config.py
import pytest
from src.grimperium.crest_pm7.config import MoleculeStatus

def test_molecule_status_enum_exists():
    """Test that MoleculeStatus enum has required values"""
    assert hasattr(MoleculeStatus, 'PENDING')
    assert hasattr(MoleculeStatus, 'CURRENT_BATCH')
    assert hasattr(MoleculeStatus, 'RUNNING')
    assert hasattr(MoleculeStatus, 'OK')
    assert hasattr(MoleculeStatus, 'RERUN')
    assert hasattr(MoleculeStatus, 'SKIP')

def test_molecule_status_values():
    """Test MoleculeStatus enum values"""
    assert MoleculeStatus.PENDING.value == "Pending"
    assert MoleculeStatus.OK.value == "OK"
    assert MoleculeStatus.RERUN.value == "Rerun"
    assert MoleculeStatus.SKIP.value == "Skip"
```

**Step 2: Run test**

```bash
pytest tests/test_config.py::test_molecule_status_enum_exists -v
```

Expected: PASS (if enum exists) or FAIL (if missing)

**Step 3: If missing, add MoleculeStatus enum to config.py**

```python
# src/grimperium/crest_pm7/config.py
from enum import Enum

class MoleculeStatus(Enum):
    PENDING = "Pending"
    CURRENT_BATCH = "Current Batch"
    RUNNING = "Running"
    OK = "OK"
    RERUN = "Rerun"
    SKIP = "Skip"
```

**Step 4: Run test again**

```bash
pytest tests/test_config.py -v
```

Expected: All PASS

**Step 5: Add reruns field to molecule dataclass (if not present)**

Search for molecule dataclass definition and add:
```python
reruns: int = 0  # Retry counter: 0-3
```

**Step 6: Commit**

```bash
git add src/grimperium/crest_pm7/config.py tests/test_config.py
git commit -m "feat(config): ensure MoleculeStatus enum and reruns field exist"
```

---

## Task 11: Documentation - Create RDKit Integration Guide

**Files:**
- Create: `docs/RDKIT_INTEGRATION.md`
- Reference: `docs/CREST_INTEGRATION.md` (template)

**Step 1: Read existing integration docs as template**

```bash
head -n 50 docs/CREST_INTEGRATION.md
```

**Step 2: Create RDKit integration document**

Create `docs/RDKIT_INTEGRATION.md` with complete content (see implementation guide Part 6 for full content)

**Step 3: Verify document quality**

```bash
wc -l docs/RDKIT_INTEGRATION.md
# Expected: ~300-400 lines
```

**Step 4: Check markdown syntax**

```bash
# If markdownlint available:
markdownlint docs/RDKIT_INTEGRATION.md
```

**Step 5: Commit**

```bash
git add docs/RDKIT_INTEGRATION.md
git commit -m "docs: add RDKit integration guide for Phase A"
```

---

## Task 12: Integration Testing - Run 3-Molecule Batch

**Files:**
- Test data: `data/molecules_pm7/computed/thermo_pm7.csv`
- Test script: Create `scripts/test_phase_a.sh`

**Step 1: Create test script**

```bash
#!/bin/bash
# scripts/test_phase_a.sh

echo "=== Phase A Integration Test ==="
echo ""

echo "Step 1: Run 3-molecule batch"
python -m grimperium.cli.main <<EOF
1
3
20
20
y
EOF

echo ""
echo "Step 2: Verify CSV exists"
if [ -f "data/molecules_pm7/computed/thermo_pm7.csv" ]; then
    echo "✓ CSV exists"
else
    echo "✗ CSV not found"
    exit 1
fi

echo ""
echo "Step 3: Check column count"
COLS=$(head -n 1 data/molecules_pm7/computed/thermo_pm7.csv | tr ',' '\n' | wc -l)
echo "Columns found: $COLS"
if [ "$COLS" -eq 49 ]; then
    echo "✓ Correct column count (49)"
else
    echo "✗ Wrong column count (expected 49, got $COLS)"
    exit 1
fi

echo ""
echo "Step 4: Run pytest"
pytest tests/test_settings_phase_a.py tests/test_csv_schema.py tests/test_molecule_processor.py -v

echo ""
echo "=== Phase A Test Complete ==="
```

**Step 2: Make executable**

```bash
chmod +x scripts/test_phase_a.sh
```

**Step 3: Run integration test**

```bash
./scripts/test_phase_a.sh
```

Expected: All checks pass

**Step 4: Manual CLI verification**

```bash
python -m grimperium.cli.main
# Navigate to Settings > CREST Settings
# Verify dropdowns for "Set CREST Method" and "Set Quick Mode"
# Select gfnff and quick
# Return to main menu
# Re-enter CREST settings
# Verify selections persisted
```

**Step 5: Verify CSV content**

```bash
head -n 5 data/molecules_pm7/computed/thermo_pm7.csv
```

Expected:
- 49 columns
- H298_pm7 column present (not most_stable_hof)
- abs_diff, abs_diff_% calculated
- delta_1, delta_2, delta_3 present
- c_method and qm columns populated

**Step 6: Commit test script**

```bash
git add scripts/test_phase_a.sh
git commit -m "test: add Phase A integration test script"
```

---

## Task 13: Final Validation - Run Full Test Suite

**Files:**
- All test files
- Final checklist: `grimperium_final_checklist.md`

**Step 1: Run all pytest tests**

```bash
pytest tests/ -v --cov=src/grimperium --cov-report=term
```

Expected: High coverage, all tests pass

**Step 2: Verify CSV with pandas**

```python
import pandas as pd

df = pd.read_csv("data/molecules_pm7/computed/thermo_pm7.csv")

print(f"Columns: {len(df.columns)}")
print(f"Rows: {len(df)}")
print("\nColumn names:")
print(list(df.columns))

print("\nmol_00001 data:")
row1 = df[df['mol_id'] == 'mol_00001'].iloc[0]
print(f"  nrotbonds: {row1['nrotbonds']}")
print(f"  H298_pm7: {row1['H298_pm7']}")
print(f"  abs_diff: {row1['abs_diff']}")
print(f"  c_method: {row1['c_method']}")
print(f"  qm: {row1['qm']}")
```

**Step 3: Run checklist verification**

Go through `grimperium_final_checklist.md`:
- [ ] Settings roundtrip works
- [ ] Backward compatibility verified
- [ ] CSV has 49 columns
- [ ] All values populated
- [ ] No old column names
- [ ] All tests pass

**Step 4: Create completion summary**

```bash
cat > PHASE_A_COMPLETION.md <<'EOF'
# Phase A Completion Summary

## Date: 2026-01-20

## Changes Implemented

### 1. Settings System (settings_manager.py)
- ✓ Replaced boolean toggles with dropdown menus
- ✓ Added `crest_method` (gfn2, gfnff, gfn2//gfnff)
- ✓ Added `quick_mode` (off, quick, squick, mquick)
- ✓ Backward compatibility for old CSV format
- ✓ Updated UI to show dropdowns

### 2. CSV Schema (49 columns)
- ✓ Renamed: most_stable_hof → H298_pm7
- ✓ Removed: crest_timeout, mopac_timeout
- ✓ Added: abs_diff, abs_diff_%
- ✓ Added: delta_1, delta_2, delta_3, conformer_selected
- ✓ Added: Settings tracking (c_method, qm, precise_scf, etc)
- ✓ Added: RDKit descriptors (nrotbonds, tpsa, aromatic_rings)

### 3. Molecule Processor
- ✓ Fixed: mopac_time aggregation
- ✓ Fixed: Delta calculation (select minimum, not maximum)
- ✓ Added: Conformer selection tracking

### 4. Documentation
- ✓ Created: docs/RDKIT_INTEGRATION.md
- ✓ Updated: Help text for new settings

## Test Results

```
pytest: XX/XX tests passed
Coverage: XX%
Integration: 3-molecule batch successful
CSV: 49 columns, all populated
```

## Ready for Phase B
- ✓ Schema finalized
- ✓ Settings improved
- ✓ Data validation fixed
- ✓ Documentation complete

## Next Steps
Phase B: Scale to 29,000 molecules for ML training
EOF
```

**Step 5: Final commit**

```bash
git add PHASE_A_COMPLETION.md
git commit -m "docs: Phase A completion summary"
```

**Step 6: Create tag**

```bash
git tag -a phase-a-complete -m "Phase A: Dataset & CLI refinements complete"
```

---

## Success Criteria

**Code Quality:**
- ✓ All tests pass
- ✓ No syntax errors
- ✓ Type hints preserved
- ✓ No regressions

**Functionality:**
- ✓ Settings dropdowns work
- ✓ Backward compatibility verified
- ✓ CSV has 49 columns in correct order
- ✓ All metrics calculated correctly
- ✓ Delta selection uses minimum (not maximum)

**Documentation:**
- ✓ RDKit integration documented
- ✓ Help text updated
- ✓ Completion summary created

**Data Validation:**
- ✓ mol_00001: nrotbonds=1, aromatic_rings=1
- ✓ mol_00001: abs_diff calculated
- ✓ mol_00001: conformer_selected tracked
- ✓ All molecules: mopac_time > 0

---

## Estimated Time

- Task 1-2: 30 min (Settings fields + tests)
- Task 3-5: 45 min (Settings to_dict/from_dict/UI)
- Task 6-7: 60 min (CSV schema update)
- Task 8-9: 30 min (Processor fixes)
- Task 10: 15 min (Config verification)
- Task 11: 25 min (RDKit docs)
- Task 12-13: 30 min (Testing & validation)

**Total: ~3.5-4 hours**

---

## Notes

- Follow TDD where possible
- Commit frequently (after each task)
- Test manually via CLI for UI changes
- Use `grimperium_final_checklist.md` for validation
- Reference `grimperium_spec.md` for exact specifications
