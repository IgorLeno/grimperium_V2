# 📋 DOCUMENTO 2: GRIMPERIUM IMPLEMENTATION GUIDE

```markdown
# GRIMPERIUM V2 - IMPLEMENTATION GUIDE FOR CLAUDE

**Purpose:** Step-by-step surgical implementation instructions  
**Target Audience:** Claude AI (via Context7)  
**Date:** 2026-01-20  
**Status:** Ready for Implementation  

---

## PART 1: REPOSITORY STRUCTURE & FILE LOCATIONS

### Key Directories

```
grimperium_V2/
├── src/
│   └── grimperium/
│       ├── cli/
│       │   ├── __init__.py
│       │   ├── main.py                    ← Entry point
│       │   ├── settings_manager.py        ← 🔴 PRIMARY WORK HERE (70% of changes)
│       │   ├── batch_menu.py
│       │   └── ...
│       ├── crest_pm7/
│       │   ├── __init__.py
│       │   ├── config.py                  ← Status enum (verify)
│       │   ├── molecule_processor.py      ← 🟡 SECONDARY WORK (10% of changes)
│       │   ├── crest_handler.py
│       │   ├── mopac_handler.py
│       │   ├── batch_processor.py         ← 🟡 CSV WRITER LOCATION (TBD - 15% of changes)
│       │   └── ...
│       └── rdkit_handler.py               ← RDKit integration (verify it exists)
├── docs/
│       ├── CREST_INTEGRATION.md           ← Template for RDKit doc
│       ├── MOPAC_INTEGRATION.md           ← Template for RDKit doc
│       └── RDKIT_INTEGRATION.md           ← 🟢 CREATE NEW (25 min)
├── data/
│       └── molecules_pm7/
│           ├── computed/
│           │   └── thermo_pm7.csv         ← Test data (3 molecules)
│           └── ...
├── tests/
│       └── test_phase_a_validation.py    ← Validation tests
└── ...
```

### Files to Search For

**To find CSV writer:**
```bash
grep -r "to_csv\|\.csv\|thermo_pm7" src/grimperium/ --include="*.py"
grep -r "H298_pm7\|most_stable_hof" src/ --include="*.py"
grep -r "DataFrame\|pd.DataFrame" src/grimperium/crest_pm7/ --include="*.py"
```

**Expected matches:**
- `batch_processor.py` - Most likely location
- `database_manager.py` - Possible location
- `results_handler.py` - Possible location

---

## PART 2: SETTINGS SYSTEM OVERHAUL

### FILE: `src/grimperium/cli/settings_manager.py`

#### 2.1 Identify Current Structure

**Search for:**
```python
@dataclass
class CRESTSettings:
```

**You should find (~line 25-50):**
```python
@dataclass
class CRESTSettings:
    xtb: bool = False
    v3: bool = False
    gfnff: bool = False      # ← TO REMOVE
    quick: bool = False      # ← TO REMOVE
    nci: bool = False
    energy_window: float = 5.0
    rmsd_threshold: float = 1.5
    threads: int = 4
```

#### 2.2 STEP 1: Add New Fields to CRESTSettings

**BEFORE:**
```python
@dataclass
class CRESTSettings:
    xtb: bool = False
    v3: bool = False
    gfnff: bool = False      # DELETE THIS LINE
    quick: bool = False      # DELETE THIS LINE
    nci: bool = False
    energy_window: float = 5.0
    rmsd_threshold: float = 1.5
    threads: int = 4
```

**AFTER:**
```python
@dataclass
class CRESTSettings:
    # Valid values for dropdowns
    CREST_METHOD_OPTIONS = ["gfn2", "gfnff", "gfn2//gfnff"]
    QUICK_MODE_OPTIONS = ["off", "quick", "squick", "mquick"]
    
    # Settings
    xtb: bool = False
    v3: bool = False
    nci: bool = False
    crest_method: str = "gfn2"              # NEW: was gfnff boolean
    quick_mode: str = "off"                 # NEW: was quick boolean
    energy_window: float = 5.0
    rmsd_threshold: float = 1.5
    threads: int = 4
```

**Verification:**
```python
# Test defaults work
settings = CRESTSettings()
assert settings.crest_method == "gfn2"
assert settings.quick_mode == "off"
```

#### 2.3 STEP 2: Find & Update `to_dict()` Method

**Search for:** `def to_dict(self):`

**You should find (~line 100-150):**
```python
def to_dict(self) -> dict:
    return {
        "crest_xtb": str(self.crest.xtb),
        "crest_v3": str(self.crest.v3),
        "crest_gfnff": str(self.crest.gfnff),      # ← REMOVE THIS
        "crest_quick": str(self.crest.quick),      # ← REMOVE THIS
        "crest_nci": str(self.crest.nci),
        ...
    }
```

**CHANGE TO:**
```python
def to_dict(self) -> dict:
    return {
        "crest_xtb": str(self.crest.xtb),
        "crest_v3": str(self.crest.v3),
        "crest_method": self.crest.crest_method,        # NEW
        "crest_quick_mode": self.crest.quick_mode,      # NEW
        "crest_nci": str(self.crest.nci),
        ...
    }
```

**DO NOT include:**
- `crest_gfnff`
- `crest_quick`

#### 2.4 STEP 3: Update `from_dict()` Method

**Search for:** `def from_dict(self, data: dict):`

**You should find (~line 150-200):**
```python
def from_dict(self, data: dict) -> None:
    # Conversion mappings
    conversion_map = {
        ("crest_xtb", "xtb"): ("crest", "xtb"),
        ("crest_v3", "v3"): ("crest", "v3"),
        ("crest_gfnff", "gfnff"): ("crest", "gfnff"),  # ← HANDLE THIS
        ("crest_quick", "quick"): ("crest", "quick"),  # ← HANDLE THIS
        ...
    }
    
    for (old_key, alt_key), (section, attr) in conversion_map.items():
        if old_key in data:
            value = self._parse_bool(data[old_key])
            ...
```

**CHANGE TO (add backward compatibility):**
```python
def from_dict(self, data: dict) -> None:
    # Backward compatibility: old boolean toggles → new dropdown values
    
    # Handle old crest_gfnff toggle
    if "crest_gfnff" in data and "crest_method" not in data:
        if self._parse_bool(data["crest_gfnff"]):
            data["crest_method"] = "gfnff"
        else:
            data["crest_method"] = "gfn2"
    
    # Handle old crest_quick toggle
    if "crest_quick" in data and "crest_quick_mode" not in data:
        if self._parse_bool(data["crest_quick"]):
            data["crest_quick_mode"] = "quick"
        else:
            data["crest_quick_mode"] = "off"
    
    # Standard conversion mappings (keep existing ones)
    conversion_map = {
        ("crest_xtb", "xtb"): ("crest", "xtb"),
        ("crest_v3", "v3"): ("crest", "v3"),
        ("crest_nci", "nci"): ("crest", "nci"),
        ("crest_method",): ("crest", "crest_method"),         # NEW
        ("crest_quick_mode",): ("crest", "quick_mode"),       # NEW
        ...
    }
    
    for (old_key, *alt_keys), (section, attr) in conversion_map.items():
        if old_key in data:
            # For method and quick_mode, use string directly
            if attr in ["crest_method", "quick_mode"]:
                value = data[old_key]
                # Validate
                if attr == "crest_method":
                    assert value in CRESTSettings.CREST_METHOD_OPTIONS
                elif attr == "quick_mode":
                    assert value in CRESTSettings.QUICK_MODE_OPTIONS
            else:
                # For boolean settings, parse as before
                value = self._parse_bool(data[old_key])
            
            # Set value
            if section == "crest":
                setattr(self.crest, attr, value)
            elif section == "mopac":
                setattr(self.mopac, attr, value)
```

**Verification:**
```python
# Test backward compatibility
sm = SettingsManager()

# Old format
old_data = {
    "crest_gfnff": "true",
    "crest_quick": "true"
}
sm.from_dict(old_data)
assert sm.crest.crest_method == "gfnff"
assert sm.crest.quick_mode == "quick"

# New format
new_data = {
    "crest_method": "gfn2//gfnff",
    "crest_quick_mode": "squick"
}
sm2 = SettingsManager()
sm2.from_dict(new_data)
assert sm2.crest.crest_method == "gfn2//gfnff"
assert sm2.crest.quick_mode == "squick"
```

#### 2.5 STEP 4: Update HELP_TEXT Dictionary

**Search for:** `HELP_TEXT = {`

**You should find (~line 250-300):**
```python
HELP_TEXT = {
    "crest_gfnff": "Enable GFN-FF force field...",
    "crest_quick": "Enable quick mode...",
    ...
}
```

**CHANGE TO:**
```python
HELP_TEXT = {
    "crest_method": "Choose CREST quantum method: gfn2 (default, balanced), gfnff (faster), gfn2//gfnff (two-step refinement)",
    "crest_quick_mode": "Choose speed/accuracy tradeoff: off (full), quick (fast), squick (super-fast), mquick (fastest)",
    ...
}
```

#### 2.6 STEP 5: Update show_crest_summary() Method

**Search for:** `def show_crest_summary(self):`

**You should find (~line 350-400):**
```python
def show_crest_summary(self) -> None:
    table = Table(title="CREST Configuration")
    ...
    table.add_row("GFN-FF Force Field", "✓" if self.crest.gfnff else "✗")
    table.add_row("Quick Mode", "✓" if self.crest.quick else "✗")
    ...
    console.print(table)
```

**CHANGE TO:**
```python
def show_crest_summary(self) -> None:
    table = Table(title="CREST Configuration")
    ...
    table.add_row("CREST Method", self.crest.crest_method.upper())
    table.add_row("Quick Mode", self.crest.quick_mode)
    ...
    console.print(table)
```

#### 2.7 STEP 6: Replace display_crest_menu() Method

**Search for:** `def display_crest_menu(self):`

**You should find (~line 400-500):**
```python
def display_crest_menu(self) -> None:
    while True:
        console.print("\n[bold cyan]CREST Settings[/bold cyan]")
        
        # Current toggles (KEEP v3, nci as toggles)
        if questionary.confirm("Toggle xTB method?").ask():
            self.crest.xtb = not self.crest.xtb
        
        if questionary.confirm("Toggle v3 algorithm?").ask():
            self.crest.v3 = not self.crest.v3
        
        if questionary.confirm("Toggle GFN-FF force field?").ask():    # ← REMOVE
            self.crest.gfnff = not self.crest.gfnff
        
        if questionary.confirm("Toggle Quick Mode?").ask():             # ← REMOVE
            self.crest.quick = not self.crest.quick
        
        if questionary.confirm("Toggle NCI mode?").ask():
            self.crest.nci = not self.crest.nci
        
        if questionary.confirm("Back to main menu?").ask():
            break
```

**CHANGE TO:**
```python
def display_crest_menu(self) -> None:
    while True:
        console.print("\n[bold cyan]CREST Settings[/bold cyan]")
        
        # Show current settings
        console.print(f"Current CREST Method: {self.crest.crest_method}")
        console.print(f"Current Quick Mode: {self.crest.quick_mode}")
        
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
        
        elif choice == "Set Quick Mode":
            self.crest.quick_mode = questionary.select(
                "Choose quick mode:",
                choices=["off", "quick", "squick", "mquick"]
            ).ask()
        
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

**Verification:**
```bash
# Run CLI and navigate to CREST settings
# Expected behavior:
# - See "Set CREST Method" option
# - See "Set Quick Mode" option
# - Can select from dropdown (not toggle)
# - Selections persist
```

---

## PART 3: CSV WRITER LOCATION & SCHEMA UPDATE

### Step 1: Locate CSV Writer

**Run these commands to find:**

```bash
# Search 1: Find CSV output function
grep -r "\.to_csv\|to_csv(" src/grimperium/ --include="*.py" -n

# Search 2: Find thermo_pm7 file handling
grep -r "thermo_pm7" src/ --include="*.py" -n

# Search 3: Find DataFrame conversion
grep -r "pd\.DataFrame\|DataFrame\(" src/grimperium/crest_pm7/ --include="*.py" -n

# Search 4: Find most_stable_hof references
grep -r "most_stable_hof" src/ --include="*.py" -n
```

**Expected output locations:**
- `batch_processor.py` (most likely)
- `database_manager.py` (possible)
- `results_exporter.py` (possible)

**Document the exact location:**
```
CSV Writer File: ___________________________________
Function Name: ____________________________________
Line Number: ______________________________________
```

### Step 2: Identify Current Column Mapping

In the CSV writer, find the code that creates column mappings:

**You might see:**
```python
def to_dataframe(self):
    rows = []
    for molecule in self.molecules:
        row = {
            "mol_id": molecule.id,
            "status": molecule.status,
            "smiles": molecule.smiles,
            ...
            "most_stable_hof": molecule.hof,         # ← RENAME
            "crest_timeout": 20,                     # ← REMOVE
            "mopac_timeout": 20,                     # ← REMOVE
            ...
        }
        rows.append(row)
    return pd.DataFrame(rows)
```

### Step 3: Update Column Mapping

**REMOVE columns:**
```python
# DELETE these lines
"crest_timeout": ...,
"mopac_timeout": ...,
```

**RENAME columns:**
```python
# CHANGE this:
"most_stable_hof": molecule.hof,

# TO this:
"H298_pm7": molecule.hof,
```

**ADD new columns:**
```python
# Add after H298_pm7:
"abs_diff": abs(molecule.H298_cbs - molecule.hof),
"abs_diff_%": (abs(molecule.H298_cbs - molecule.hof) / abs(molecule.H298_cbs)) * 100,
"reruns": molecule.reruns if hasattr(molecule, 'reruns') else 0,
"nrotbonds": molecule.nrotbonds,
"tpsa": molecule.tpsa,
"aromatic_rings": molecule.aromatic_rings,

# CREST settings
"xtb": molecule.settings.crest.xtb,
"v3": molecule.settings.crest.v3,
"qm": molecule.settings.crest.quick_mode,        # NEW name from quick
"nci": molecule.settings.crest.nci,
"c_method": molecule.settings.crest.crest_method, # NEW name from gfnff
"energy_window": molecule.settings.crest.energy_window,
"rmsd_threshold": molecule.settings.crest.rmsd_threshold,
"threads": molecule.settings.crest.threads,

# MOPAC settings
"precise_scf": molecule.settings.mopac.precise_scf,
"scf_threshold": molecule.settings.mopac.scf_threshold,

# Delta calculations
"delta_1": molecule.deltas if len(molecule.deltas) > 0 else None,
"delta_2": molecule.deltas [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/923ec5e4-a41c-424c-8ee7-10a942385a05/WORKFLOW.md) if len(molecule.deltas) > 1 else None,
"delta_3": molecule.deltas [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/08dead72-cdf8-410c-97e3-0750ce889a3d/MOPAC_INTEGRATION.md) if len(molecule.deltas) > 2 else None,
"conformer_selected": molecule.conformer_selected,

# Other metrics
"mopac_time": molecule.mopac_time,  # Aggregated
```

**SET COLUMN ORDER** (exactly as specified):

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
df = df[column_order]  # Reorder columns
df.to_csv("path/to/thermo_pm7.csv", index=False)
```

**Verification:**
```python
# After writing CSV
df = pd.read_csv("path/to/thermo_pm7.csv")
assert len(df.columns) == 49, f"Expected 49 columns, got {len(df.columns)}"
assert list(df.columns) == column_order
assert "crest_timeout" not in df.columns
assert "mopac_timeout" not in df.columns
assert "H298_pm7" in df.columns
```

---

## PART 4: MOLECULE PROCESSOR FIXES

### FILE: `src/grimperium/crest_pm7/molecule_processor.py`

#### 4.1 Find mopac_time Aggregation Point

**Search for:** `mopac_time` or `mopac_execution_time`

**You should find (~line 200-250):**
```python
class ConformerResult:
    mopac_execution_time: float
    hof: float

class MoleculeResult:
    conformers: List[ConformerResult]
    
    def to_dict(self):
        return {
            "mopac_time": ???  # Currently empty or single value
        }
```

#### 4.2 Add Aggregation Logic

**CHANGE:**
```python
def to_dict(self):
    # Before: single or missing
    result = {
        "mopac_time": self.conformers.mopac_execution_time,
        ...
    }
```

**TO:**
```python
def to_dict(self):
    # Aggregate mopac_time across all successful conformers
    total_mopac_time = 0
    for conformer in self.conformers:
        if hasattr(conformer, 'mopac_execution_time') and conformer.mopac_execution_time > 0:
            total_mopac_time += conformer.mopac_execution_time
    
    result = {
        "mopac_time": total_mopac_time,  # Sum of all conformers
        ...
    }
```

**Verification:**
```python
# Test aggregation
result = processor.process(molecule)
result_dict = result.to_dict()
assert result_dict["mopac_time"] > 0, "mopac_time should be non-zero"
```

#### 4.3 Fix Delta Calculations

**Search for:** `delta_e_12` or `calculate_delta`

**You should find (~line 150-200):**
```python
def calculate_deltas(self, hofs: List[float], h298_cbs: float):
    # Current logic (might be wrong)
    deltas = {
        "delta_e_12": abs(h298_cbs - hofs) - abs(h298_cbs - hofs [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/923ec5e4-a41c-424c-8ee7-10a942385a05/WORKFLOW.md)),
        "delta_e_13": abs(h298_cbs - hofs) - abs(h298_cbs - hofs [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/08dead72-cdf8-410c-97e3-0750ce889a3d/MOPAC_INTEGRATION.md)),
        "delta_e_15": abs(h298_cbs - hofs) - abs(h298_cbs - hofs [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/05baed7e-48cc-45ea-9fa0-2771db313fc8/DATASETS.md)),
    }
```

**CHANGE TO (proper delta calculation):**
```python
def calculate_deltas(self, hofs: List[float], h298_cbs: float):
    # Correct logic: delta = |H298_cbs - hof|
    # Take top 3 conformers (lowest energy first)
    
    sorted_hofs = sorted(hofs)[:3]  # Top 3
    
    deltas = {
        "delta_1": abs(h298_cbs - sorted_hofs),
        "delta_2": abs(h298_cbs - sorted_hofs [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/923ec5e4-a41c-424c-8ee7-10a942385a05/WORKFLOW.md)) if len(sorted_hofs) > 1 else None,
        "delta_3": abs(h298_cbs - sorted_hofs [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/08dead72-cdf8-410c-97e3-0750ce889a3d/MOPAC_INTEGRATION.md)) if len(sorted_hofs) > 2 else None,
    }
    
    # Find conformer with minimum delta
    valid_deltas = [d for d in deltas.values() if d is not None]
    min_delta = min(valid_deltas)
    best_idx = list(deltas.values()).index(min_delta)
    
    return {
        "deltas": deltas,
        "conformer_selected": best_idx + 1,  # 1-indexed
        "selected_hof": sorted_hofs[best_idx]
    }
```

**Verification:**
```python
# Test delta calculation
hofs = [-3.68, -3.52, -3.45, -2.90, -2.80]
h298_cbs = 0.15234

result = processor.calculate_deltas(hofs, h298_cbs)

# Top 3 sorted HOFs
assert result["deltas"]["delta_1"] == abs(0.15234 - (-3.68))  # 3.83234
assert result["deltas"]["delta_2"] == abs(0.15234 - (-3.52))  # 3.67234
assert result["deltas"]["delta_3"] == abs(0.15234 - (-3.45))  # 3.60234

# Best is delta_3 (minimum)
assert result["conformer_selected"] == 3
assert result["selected_hof"] == -3.45
```

---

## PART 5: CONFIG & STATUS ENUM

### FILE: `src/grimperium/crest_pm7/config.py`

#### 5.1 Verify MoleculeStatus Enum

**Search for:** `class MoleculeStatus` or `MoleculeStatus = `

**You should find (~line 20-50):**
```python
from enum import Enum

class MoleculeStatus(Enum):
    PENDING = "Pending"
    RUNNING = "Running"
    OK = "OK"
    RERUN = "Rerun"
    SKIP = "Skip"
```

#### 5.2 If Missing, Add It

**If NOT found, add:**
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

#### 5.3 Update Rerun Logic

In batch processor, when molecule fails:

**Find:** `if not success:`

**Add logic:**
```python
if not success:
    molecule.reruns += 1
    
    if molecule.reruns >= 3:
        molecule.status = MoleculeStatus.SKIP
        molecule.error_message = f"Failed after {molecule.reruns} retries"
    else:
        molecule.status = MoleculeStatus.RERUN
        molecule.error_message = f"Retry attempt {molecule.reruns}"
```

---

## PART 6: CREATE RDKIT INTEGRATION DOC

### FILE: `docs/RDKIT_INTEGRATION.md` (NEW)

**Copy structure from:** `docs/CREST_INTEGRATION.md` or `docs/MOPAC_INTEGRATION.md`

**Create new file with sections:**

```markdown
# RDKit Integration

## Overview

RDKit (RDKit: Open-source cheminformatics software) is used in Phase A to compute molecular descriptors from SMILES strings before conformer generation.

## Key Capabilities

- Molecular descriptor computation
- SMILES parsing and validation
- Aromaticity perception
- Topological analysis

## Phase A Integration

RDKit is called at the very beginning of the pipeline:
1. Parse SMILES input
2. Compute 3 descriptors: nrotbonds, tpsa, aromatic_rings
3. Pass molecule to CREST for conformer generation

## Descriptors Computed

### 1. Number of Rotatable Bonds (nrotbonds)

**Definition:** Count of single bonds (excluding ring bonds and bonds to terminal atoms) that can rotate.

**RDKit Function:** `rdkit.Chem.Descriptors.NumRotatableBonds(mol)`

**Use Case:** Indicator of molecular flexibility. Molecules with more rotatable bonds tend to have larger conformer spaces.

**Typical Range:** 0-15 (most drug-like molecules)

**Example:**
- Ethane (CH3-CH3): 1
- Butane (CH3-CH2-CH2-CH3): 3
- Cyclohexane (rigid ring): 0

### 2. Topological Polar Surface Area (tpsa)

**Definition:** Sum of surfaces of polar atoms in a molecule.

**RDKit Function:** `rdkit.Chem.Descriptors.TPSA(mol)`

**Use Case:** Predicts drug-like properties and membrane permeability.

**Typical Range:** 0-200 Ų

**Example:**
- Benzene (C6H6): 0
- Water (H2O): 20
- Glucose: ~130

### 3. Aromatic Rings (aromatic_rings)

**Definition:** Count of aromatic rings in the molecule.

**RDKit Code:**
```python
aromatic_rings = sum(
    ring.IsAromatic() 
    for ring in Chem.GetSSSR(mol)
)
```

**Use Case:** Indicator of molecular rigidity and stability.

**Typical Range:** 0-5 (most molecules)

**Example:**
- Benzene: 1
- Naphthalene: 2
- Pyridine: 1

## Workflow

Input SMILES → RDKit Parse → Compute Descriptors → Output {nrotbonds, tpsa, aromatic_rings}

## Implementation

**File:** `src/grimperium/rdkit_handler.py`

**Key Function:**
```python
def compute_descriptors(smiles: str) -> Dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    return {
        "nrotbonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": Descriptors.TPSA(mol),
        "aromatic_rings": count_aromatic_rings(mol)
    }
```

## Performance

- **Execution Time:** <1ms per molecule
- **Scalability:** Linear with number of molecules
- **Bottleneck:** SMILES parsing (not descriptor computation)

## Limitations & Future Work

### Current Limitations
- Only 3 descriptors computed
- No 3D information (2D only)
- No molecular fingerprints

### Phase B Enhancements
- Additional RDKit descriptors (H-bond donors, rotatable bonds, ring count)
- Morgan fingerprints (2048-bit)
- MACCS keys fingerprints
- Molecular weight, lipophilicity (logP) predictions

## References

**RDKit Documentation:**
- Main: https://www.rdkit.org/
- Descriptors: https://www.rdkit.org/docs/GettingStartedInPython.html
- Descriptor List: https://www.rdkit.org/docs/source/rdkit.Chem.html#descriptors

**Scientific References:**
- Lipinski, C. A. (2000). Drug‐like properties and the causes of poor solubility and poor permeability. Journal of Pharmacological and Toxicological Methods, 44(1), 235-249.

## Status

✅ Phase A: Fully integrated  
🚀 Phase B: Expansion planned

***
**Last Updated:** 2026-01-20  
**Maintainer:** Grimperium Team
```

---

## PART 7: UPDATE HELP TEXT

### FILE: Various CLI files

**Update help for new settings:**

```
CREST Method (--crest-method):
  gfn2 (default)
    - Balanced speed and accuracy
    - Recommended for most molecules
    - ~1-2 minutes per molecule
    
  gfnff
    - Fast Gassian approximation
    - Lower accuracy but faster
    - ~30 seconds per molecule
    
  gfn2//gfnff
    - Two-step: GFN2 optimization followed by GFN-FF refinement
    - Best for large molecules
    - ~2-3 minutes per molecule

Quick Mode (--quick-mode):
  off (default)
    - Full optimization
    - Best results
    - Slower execution
    
  quick
    - Reduced iterations
    - 30-50% faster
    - Slightly lower accuracy
    
  squick
    - Super-quick mode
    - 50-70% faster
    - Noticeably lower accuracy
    
  mquick
    - Meta-quick (fastest)
    - 70-90% faster
    - Significant accuracy loss
```

---

## PART 8: TESTING CHECKLIST

After each change, verify:

### Settings Tests
```python
# Test 1: Default values
sm = SettingsManager()
assert sm.crest.crest_method == "gfn2"
assert sm.crest.quick_mode == "off"

# Test 2: to_dict
d = sm.to_dict()
assert "crest_method" in d
assert "crest_quick_mode" in d
assert "crest_gfnff" not in d
assert "crest_quick" not in d

# Test 3: from_dict (new format)
sm2 = SettingsManager()
sm2.from_dict({"crest_method": "gfnff", "crest_quick_mode": "quick"})
assert sm2.crest.crest_method == "gfnff"
assert sm2.crest.quick_mode == "quick"

# Test 4: from_dict (old format - backward compat)
sm3 = SettingsManager()
sm3.from_dict({"crest_gfnff": "true", "crest_quick": "false"})
assert sm3.crest.crest_method == "gfnff"
assert sm3.crest.quick_mode == "off"
```

### CSV Tests
```python
# Test 1: Column count
df = pd.read_csv("thermo_pm7.csv")
assert len(df.columns) == 49

# Test 2: Column names
expected = ["mol_id", "status", "smiles", ..., "assigned_mopac_timeout"]
assert list(df.columns) == expected

# Test 3: No old columns
assert "crest_timeout" not in df.columns
assert "mopac_timeout" not in df.columns
assert "most_stable_hof" not in df.columns

# Test 4: Metrics present
assert "abs_diff" in df.columns
assert "abs_diff_%" in df.columns

# Test 5: Deltas present
assert "delta_1" in df.columns
assert "delta_2" in df.columns
assert "delta_3" in df.columns
assert "conformer_selected" in df.columns

# Test 6: Settings present
assert "crest_method" in df.columns
assert "quick_mode" in df.columns
assert "precise_scf" in df.columns
```

---

## PART 9: QUICK REFERENCE - WHAT CHANGES WHERE

| File | Change | Lines | Type |
|------|--------|-------|------|
| settings_manager.py | Add c_method, quick_mode fields | 3 | Add |
| settings_manager.py | Remove gfnff, quick fields | 2 | Delete |
| settings_manager.py | Update to_dict | 2 | Modify |
| settings_manager.py | Update from_dict + compat | 20 | Modify |
| settings_manager.py | Update HELP_TEXT | 2 | Modify |
| settings_manager.py | Update show_crest_summary | 2 | Modify |
| settings_manager.py | Replace display_crest_menu | 40 | Replace |
| CSV Writer | Remove timeout columns | 2 | Delete |
| CSV Writer | Rename most_stable_hof | 1 | Rename |
| CSV Writer | Add metrics | 2 | Add |
| CSV Writer | Add deltas | 3 | Add |
| CSV Writer | Add settings cols | 12 | Add |
| CSV Writer | Reorder columns | 1 | Modify |
| molecule_processor.py | Add mopac_time aggregation | 5 | Add |
| molecule_processor.py | Fix delta calculation | 10 | Modify |
| config.py | Verify/Add MoleculeStatus enum | 10 | Verify |
| docs/RDKIT_INTEGRATION.md | Create new doc | 400 | New |

**Total Implementation Time:** ~2.5-3 hours

---

**Version:** 1.0  
**Status:** ✅ Ready  
**For:** Claude AI  
```




