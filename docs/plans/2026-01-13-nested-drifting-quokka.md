# CREST PM7 Batch Pipeline - Implementation Plan

## Overview

Implement a batch execution system for the CREST PM7 pipeline to process ~30,000 molecules from `thermo_cbs_clean.csv` with:
- **Sequential processing** (one molecule at a time)
- **Fixed timeouts** per batch (no dynamic prediction)
- **CSV-based tracking** with JSON detail files per molecule

---

## Architecture Summary

```
src/grimperium/crest_pm7/batch/
├── __init__.py
├── enums.py              # MoleculeStatus, BatchSortingStrategy, BatchFailurePolicy
├── models.py             # Pydantic models for CSV rows, JSON details, Batch
├── csv_manager.py        # BatchCSVManager
├── detail_manager.py     # ConformerDetailManager
├── execution_manager.py  # BatchExecutionManager
└── processor_adapter.py  # Adapter for fixed timeout processing

scripts/
└── init_batch_csv.py     # Initialize CSV from thermo_cbs_clean.csv

src/grimperium/cli/views/
└── batch_view.py         # CLI view for batch commands (or extend calc_view.py)

tests/unit/
└── test_batch/
    ├── test_csv_manager.py
    ├── test_detail_manager.py
    └── test_execution_manager.py
```

---

## Critical Files to Modify/Create

### New Files (9 files)
| File | Purpose |
|------|---------|
| `src/grimperium/crest_pm7/batch/__init__.py` | Package exports |
| `src/grimperium/crest_pm7/batch/enums.py` | New enums (3) |
| `src/grimperium/crest_pm7/batch/models.py` | Pydantic models (6) |
| `src/grimperium/crest_pm7/batch/csv_manager.py` | CSV tracking (~250 lines) |
| `src/grimperium/crest_pm7/batch/detail_manager.py` | JSON details (~100 lines) |
| `src/grimperium/crest_pm7/batch/execution_manager.py` | Batch orchestration (~200 lines) |
| `src/grimperium/crest_pm7/batch/processor_adapter.py` | Fixed timeout adapter (~80 lines) |
| `scripts/init_batch_csv.py` | CSV initialization (~150 lines) |
| `src/grimperium/cli/views/batch_view.py` | CLI interface (~200 lines) |

### Files to Modify (2 files)
| File | Change |
|------|--------|
| `src/grimperium/crest_pm7/__init__.py` | Export batch module |
| `src/grimperium/cli/app.py` | Register BatchView |

---

## Implementation Steps

### Step 1: Create Enums (`batch/enums.py`)

```python
class MoleculeStatus(str, Enum):
    PENDING = "Pending"
    TO_RUN = "To run"
    RUNNING = "Running"
    RERUN = "Rerun"
    SKIP = "Skip"
    OK = "OK"

class BatchSortingStrategy(str, Enum):
    RERUN_FIRST_THEN_EASY = "rerun_first_then_easy"
    RANDOM = "random"
    BY_NHEAVY = "by_nheavy"
    BY_NROTBONDS = "by_nrotbonds"

class BatchFailurePolicy(str, Enum):
    PARTIAL_OK = "partial_ok"
    ALL_OR_NOTHING = "all_or_nothing"
```

### Step 2: Create Pydantic Models (`batch/models.py`)

Models needed:
1. **BatchRowCSV** - 36 columns matching CSV schema
2. **ConformerDetail** - Per-conformer data for JSON
3. **Decision** - Audit trail entry
4. **ConformerDetailFile** - Complete JSON structure
5. **Batch** - Batch definition with timeout config
6. **BatchResult** - Execution summary

Key design: Reuse existing `PM7Result` fields but add batch-specific columns (status, retry_count, batch_id, etc.)

### Step 3: Create BatchCSVManager (`batch/csv_manager.py`)

Core methods:
```python
class BatchCSVManager:
    def __init__(self, csv_path: Path)
    def initialize_from_cbs_dataset(self, cbs_path, crest_timeout_min, mopac_timeout_min)
    def load_csv(self) -> pd.DataFrame
    def save_csv(self) -> None
    def select_batch(self, batch_id, batch_size, crest_timeout_min, mopac_timeout_min, strategy) -> Batch
    def mark_running(self, mol_id: str) -> None
    def mark_success(self, mol_id: str, result: PM7Result, crest_timeout_used, mopac_timeout_used) -> None
    def mark_rerun(self, mol_id: str, error_message: str) -> None
    def mark_skip(self, mol_id: str, reason: str) -> None
    def get_batch_stats(self, batch_id: str) -> dict
```

### Step 4: Create ConformerDetailManager (`batch/detail_manager.py`)

```python
class ConformerDetailManager:
    def __init__(self, detail_dir: Path)
    def create_detail_file(self, mol_id, batch_id, result: PM7Result) -> ConformerDetailFile
    def load_detail_file(self, mol_id: str) -> Optional[ConformerDetailFile]
    def update_decision(self, mol_id: str, decision: Decision) -> None
    def delete_detail_file(self, mol_id: str) -> None
```

---

## JSON Error Scenarios (CRÍTICO #3)

### Scenario 1: CREST Failure (No Conformers Generated)

```json
{
  "mol_id": "MOL00042",
  "batch_id": "batch_20260113_100000",
  "conformers": [],
  "decisions": [
    {
      "decision": "crest_failed",
      "reason": "CREST timeout after 600s",
      "confidence": "HIGH",
      "timestamp": "2026-01-13T10:10:00Z"
    }
  ],
  "issues": [
    "CREST_TIMEOUT: Exceeded 10 minute limit",
    "NO_CONFORMERS: CREST did not generate any conformers"
  ],
  "quality_grade": "FAILED"
}
```

**CSV Status**: `Rerun` (if retry_count < max_retries) or `Skip` (if exhausted)

### Scenario 2: CREST Success, MOPAC Partial Success

```json
{
  "mol_id": "MOL00123",
  "batch_id": "batch_20260113_100000",
  "conformers": [
    {
      "index": 0,
      "crest_status": "SUCCESS",
      "crest_geometry_file": "/tmp/crest_pm7/MOL00123/MOL00123_conf000.xyz",
      "mopac_status": "SUCCESS",
      "mopac_output_file": "/tmp/crest_pm7/MOL00123/MOL00123_conf000.out",
      "mopac_execution_time": 45.3,
      "mopac_timeout_used": 600.0,
      "energy_hof": -56.78,
      "hof_confidence": "HIGH",
      "hof_extraction_method": "FINAL_HEAT_OF_FORMATION",
      "is_successful": true
    },
    {
      "index": 1,
      "crest_status": "SUCCESS",
      "crest_geometry_file": "/tmp/crest_pm7/MOL00123/MOL00123_conf001.xyz",
      "mopac_status": "SCF_FAILED",
      "mopac_output_file": "/tmp/crest_pm7/MOL00123/MOL00123_conf001.out",
      "mopac_execution_time": 120.5,
      "mopac_timeout_used": 554.7,
      "mopac_error_message": "SCF DID NOT CONVERGE",
      "energy_hof": null,
      "hof_confidence": "LOW",
      "hof_extraction_method": null,
      "is_successful": false
    },
    {
      "index": 2,
      "crest_status": "SUCCESS",
      "crest_geometry_file": "/tmp/crest_pm7/MOL00123/MOL00123_conf002.xyz",
      "mopac_status": "TIMEOUT",
      "mopac_output_file": "/tmp/crest_pm7/MOL00123/MOL00123_conf002.out",
      "mopac_execution_time": 434.2,
      "mopac_timeout_used": 434.2,
      "mopac_error_message": "MOPAC timeout after 434.2s",
      "energy_hof": null,
      "hof_confidence": "LOW",
      "hof_extraction_method": null,
      "is_successful": false
    }
  ],
  "decisions": [
    {"decision": "selected_3_conformers", "reason": "nrotbonds=3 (medium flexibility)", "confidence": "HIGH", "timestamp": "..."},
    {"decision": "used_fixed_crest_timeout", "reason": "batch-level timeout: 10 min", "confidence": "HIGH", "timestamp": "..."},
    {"decision": "partial_success", "reason": "1/3 conformers succeeded", "confidence": "MEDIUM", "timestamp": "..."}
  ],
  "issues": [
    "MOPAC_SCF_FAILED: Conformer 1 SCF did not converge",
    "MOPAC_TIMEOUT: Conformer 2 exceeded timeout"
  ],
  "quality_grade": "B"
}
```

**CSV Status**: `OK` (at least 1 conformer succeeded)
**Quality Grade**: `B` (partial success)

### Scenario 3: CREST Success, All MOPAC Failed

```json
{
  "mol_id": "MOL00456",
  "batch_id": "batch_20260113_100000",
  "conformers": [
    {
      "index": 0,
      "crest_status": "SUCCESS",
      "crest_geometry_file": "/tmp/crest_pm7/MOL00456/MOL00456_conf000.xyz",
      "mopac_status": "GEOMETRY_ERROR",
      "mopac_error_message": "INITIAL GEOMETRY WAS BAD",
      "is_successful": false
    }
  ],
  "decisions": [
    {"decision": "selected_1_conformer", "reason": "rigid molecule (nrotbonds=0)", "confidence": "HIGH", "timestamp": "..."},
    {"decision": "all_conformers_failed", "reason": "0/1 conformers succeeded", "confidence": "HIGH", "timestamp": "..."}
  ],
  "issues": [
    "MOPAC_GEOMETRY_ERROR: Invalid initial geometry from CREST",
    "NO_HOF_EXTRACTED: No successful energy extraction"
  ],
  "quality_grade": "FAILED"
}
```

**CSV Status**: `Rerun` or `Skip`

### Scenario 4: Invalid SMILES (Pre-Processing Failure)

```json
{
  "mol_id": "MOL00789",
  "batch_id": "batch_20260113_100000",
  "conformers": [],
  "decisions": [
    {"decision": "smiles_parsing_failed", "reason": "RDKit could not parse SMILES", "confidence": "HIGH", "timestamp": "..."}
  ],
  "issues": [
    "INVALID_SMILES: Failed to parse 'C[C@H](invalid)C'",
    "PREPROCESSING_FAILED: Could not generate 3D coordinates"
  ],
  "quality_grade": "FAILED"
}
```

**CSV Status**: `Skip` (no point retrying invalid SMILES)

### Error Handling Logic in `create_detail_file()`

```python
def create_detail_file(
    self,
    mol_id: str,
    batch_id: str,
    result: PM7Result,
) -> ConformerDetailFile:
    """Create JSON detail file from PM7Result."""

    # Determine issues based on result state
    issues = list(result.issues)  # Copy from PM7Result

    # Add specific error categorizations
    if result.crest_status == CRESTStatus.FAILED:
        if "timeout" in (result.crest_error or "").lower():
            issues.append("CREST_TIMEOUT")
        else:
            issues.append(f"CREST_FAILED: {result.crest_error}")

    for conf in result.conformers:
        if conf.mopac_status == MOPACStatus.TIMEOUT:
            issues.append(f"MOPAC_TIMEOUT: Conformer {conf.index}")
        elif conf.mopac_status == MOPACStatus.SCF_FAILED:
            issues.append(f"MOPAC_SCF_FAILED: Conformer {conf.index}")
        elif conf.mopac_status == MOPACStatus.GEOMETRY_ERROR:
            issues.append(f"MOPAC_GEOMETRY_ERROR: Conformer {conf.index}")

    if not result.most_stable_hof:
        issues.append("NO_HOF_EXTRACTED")

    # Build detail file
    detail = ConformerDetailFile(
        mol_id=mol_id,
        batch_id=batch_id,
        conformers=[self._conformer_to_detail(c) for c in result.conformers],
        decisions=[Decision(**d) for d in self._build_decisions(result)],
        issues=issues,
        quality_grade=result.quality_grade.value,
    )

    # Save to file
    file_path = self.detail_dir / f"{mol_id}.json"
    file_path.write_text(detail.model_dump_json(indent=2))

    return detail
```

### Step 5: Create ProcessorAdapter (`batch/processor_adapter.py`)

#### Integration Analysis

Current flow in `MoleculeProcessor.process()`:
1. **CREST**: Uses `config.crest_timeout` directly (line 461 in `run_crest()`)
2. **MOPAC**: Uses `timeout_predictor.predict(nheavy, num_conformers)` to get total timeout
   - Total is then distributed across conformers: `per_conformer = remaining / remaining_conformers`

#### Solution: FixedTimeoutPredictor + Config Override

```python
from dataclasses import dataclass
from grimperium.crest_pm7.config import PM7Config, TimeoutConfidence
from grimperium.crest_pm7.molecule_processor import MoleculeProcessor, PM7Result
from grimperium.crest_pm7.timeout_predictor import TimeoutPredictor


class FixedTimeoutPredictor(TimeoutPredictor):
    """TimeoutPredictor that returns fixed values (no learning)."""

    def __init__(self, fixed_timeout_seconds: float) -> None:
        """Initialize with fixed timeout.

        Args:
            fixed_timeout_seconds: Fixed total timeout for all MOPAC conformers
        """
        super().__init__(recalibrate_interval=999999)  # Never recalibrate
        self._fixed_timeout = fixed_timeout_seconds

    def predict(
        self,
        nheavy: int,
        num_conformers: int,
    ) -> tuple[float, TimeoutConfidence]:
        """Always return fixed timeout with HIGH confidence."""
        return self._fixed_timeout, TimeoutConfidence.HIGH

    def add_observation(self, nheavy: int, execution_time: float) -> None:
        """No-op: don't learn from observations in batch mode."""
        pass


@dataclass
class FixedTimeoutProcessor:
    """Adapter for MoleculeProcessor with fixed timeouts."""

    base_config: PM7Config

    def process_with_fixed_timeout(
        self,
        mol_id: str,
        smiles: str,
        crest_timeout_seconds: float,
        mopac_timeout_seconds: float,
    ) -> PM7Result:
        """Process molecule with fixed timeouts.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            crest_timeout_seconds: Fixed CREST timeout (seconds)
            mopac_timeout_seconds: Fixed total MOPAC timeout (seconds)

        Returns:
            PM7Result with processing results
        """
        # Create config copy with fixed CREST timeout
        from copy import deepcopy
        config = deepcopy(self.base_config)
        config.crest_timeout = crest_timeout_seconds

        # Create processor with fixed MOPAC timeout predictor
        fixed_predictor = FixedTimeoutPredictor(mopac_timeout_seconds)
        processor = MoleculeProcessor(config, timeout_predictor=fixed_predictor)

        # Process and return
        return processor.process(mol_id, smiles)
```

**Key Points**:
- `FixedTimeoutPredictor` inherits from `TimeoutPredictor` for interface compatibility
- Override `predict()` to return fixed value
- Override `add_observation()` as no-op (don't learn in batch mode)
- Clone config to avoid modifying shared state

### Step 6: Create BatchExecutionManager (`batch/execution_manager.py`)

```python
class BatchExecutionManager:
    def __init__(self, csv_path, detail_dir, config: PM7Config)

    def create_batch(self, batch_id, batch_size, crest_timeout_min, mopac_timeout_min, strategy) -> Batch

    def execute_batch(self, batch: Batch) -> BatchResult:
        """Sequential execution loop"""
        for mol_id in batch.molecules:
            1. csv_manager.mark_running(mol_id)
            2. Get smiles from CSV
            3. processor_adapter.process_with_fixed_timeout(...)
            4. If success: csv_manager.mark_success(...)
            5. If failed: csv_manager.mark_rerun(...) or mark_skip(...)
            6. detail_manager.create_detail_file(...)
        return BatchResult

    def handle_interrupt(self, batch: Batch) -> None
    def display_batch_summary(self, result: BatchResult) -> None
```

### Step 7: Create Init Script (`scripts/init_batch_csv.py`)

```python
def init_batch_csv(
    cbs_csv_path: str,
    output_csv_path: str,
    crest_timeout_minutes: int = 10,
    mopac_timeout_minutes: int = 30,
) -> pd.DataFrame:
    """
    Generate crest_pm7_batch_status.csv from thermo_cbs_clean.csv

    - Map index to mol_id (e.g., MOL00001)
    - Calculate RDKit descriptors
    - Initialize status = Pending
    """
```

### Step 8: Create CLI BatchView (`cli/views/batch_view.py`)

```python
class BatchView(BaseView):
    name = "batch"
    title = "Batch Processing"

    def run_batch(self, batch_size, crest_timeout, mopac_timeout, strategy) -> None
    def show_status(self) -> None
    def resume_batch(self, batch_id: str) -> None
```

### Step 9: Register in CLI (`cli/app.py`)

Add BatchView to `_register_views()` method.

### Step 10: Add Tests

Unit tests for:
- `test_csv_manager.py`: CSV operations, status transitions
- `test_detail_manager.py`: JSON file operations
- `test_execution_manager.py`: Batch execution logic (mocked processors)

---

## CSV Schema (36 columns)

| Group | Columns |
|-------|---------|
| Identification (4) | mol_id, smiles, timestamp, phase |
| RDKit Descriptors (5) | nheavy, nrotbonds, tpsa, aromatic_rings, has_heteroatoms |
| Status Tracking (4) | status, retry_count, max_retries, last_error_message |
| CREST Execution (6) | crest_status, crest_conformers_generated, crest_time, crest_error, crest_output_directory, crest_ensemble_file |
| Timeout Config (2) | crest_timeout_max_minutes, mopac_timeout_max_minutes |
| MOPAC Execution (8) | num_conformers_selected, most_stable_hof, quality_grade, success, error_message, total_execution_time, actual_crest_timeout_used, actual_mopac_timeout_used |
| Delta-E (3) | delta_e_12, delta_e_13, delta_e_15 |
| Batch Tracking (4) | batch_id, batch_order, batch_sorting_strategy, batch_failure_policy |

---

## PM7Result → CSV Field Mapping (CRÍTICO #2)

### Complete Field-by-Field Mapping

```python
def pm7result_to_csv_row(
    result: PM7Result,
    batch_id: str,
    batch_order: int,
    crest_timeout_max_minutes: int,
    mopac_timeout_max_minutes: int,
) -> dict:
    """Map PM7Result to CSV row dict.

    Source: src/grimperium/crest_pm7/molecule_processor.py (PM7Result dataclass)
    """
    return {
        # === Group 1: Identification ===
        "mol_id": result.mol_id,                          # PM7Result.mol_id
        "smiles": result.smiles,                          # PM7Result.smiles
        "timestamp": result.timestamp.isoformat(),        # PM7Result.timestamp (datetime → ISO string)
        "phase": result.phase,                            # PM7Result.phase (str: "A", "B", "C", "PRODUCTION")

        # === Group 2: RDKit Descriptors ===
        "nheavy": result.nheavy,                          # PM7Result.nheavy (int)
        "nrotbonds": result.nrotbonds,                    # PM7Result.nrotbonds (int)
        "tpsa": round(result.tpsa, 2) if result.tpsa else None,  # PM7Result.tpsa (float, 2 decimals)
        "aromatic_rings": result.aromatic_rings,          # PM7Result.aromatic_rings (int)
        "has_heteroatoms": result.has_heteroatoms,        # PM7Result.has_heteroatoms (bool)

        # === Group 3: Status Tracking ===
        # NOTE: These are BATCH-MANAGED, not from PM7Result
        "status": "OK" if result.success else "Rerun",    # Determined by BatchCSVManager
        "retry_count": 0,                                 # Managed by BatchCSVManager
        "max_retries": 3,                                 # Config default
        "last_error_message": result.error_message,       # PM7Result.error_message (if failed)

        # === Group 4: CREST Execution ===
        "crest_status": result.crest_status.value,        # PM7Result.crest_status (CRESTStatus enum → str)
        "crest_conformers_generated": result.crest_conformers_generated,  # PM7Result.crest_conformers_generated (int)
        "crest_time": round(result.crest_time, 1) if result.crest_time else None,  # PM7Result.crest_time (float, 1 decimal)
        "crest_error": result.crest_error,                # PM7Result.crest_error (str or None)
        "crest_output_directory": None,                   # Extracted from CREST work_dir (runtime)
        "crest_ensemble_file": None,                      # Extracted from CREST output (runtime)

        # === Group 5: Timeout Configuration ===
        "crest_timeout_max_minutes": crest_timeout_max_minutes,   # From batch config
        "mopac_timeout_max_minutes": mopac_timeout_max_minutes,   # From batch config

        # === Group 6: MOPAC Execution ===
        "num_conformers_selected": result.num_conformers_selected,  # PM7Result.num_conformers_selected (int)
        "most_stable_hof": round(result.most_stable_hof, 2) if result.most_stable_hof else None,  # PM7Result.most_stable_hof (property, float)
        "quality_grade": result.quality_grade.value,      # PM7Result.quality_grade (QualityGrade enum → str)
        "success": result.success,                        # PM7Result.success (bool)
        "error_message": result.error_message,            # PM7Result.error_message (str or None)
        "total_execution_time": round(result.total_execution_time, 1) if result.total_execution_time else None,  # PM7Result.total_execution_time (float)
        "actual_crest_timeout_used": round(result.crest_time, 1) if result.crest_time else None,  # Same as crest_time
        "actual_mopac_timeout_used": _sum_mopac_times(result),  # Sum of all conformer MOPAC times

        # === Group 7: Delta-E (Energy Differences) ===
        # Source: conformer_selector.py:calculate_delta_e()
        # E1 = lowest energy, delta_e_12 = E2 - E1, etc.
        "delta_e_12": round(result.delta_e_12, 4) if result.delta_e_12 else None,  # PM7Result.delta_e_12
        "delta_e_13": round(result.delta_e_13, 4) if result.delta_e_13 else None,  # PM7Result.delta_e_13
        "delta_e_15": round(result.delta_e_15, 4) if result.delta_e_15 else None,  # PM7Result.delta_e_15

        # === Group 8/9: Batch Tracking ===
        "batch_id": batch_id,                             # From batch execution
        "batch_order": batch_order,                       # Position in batch (1-indexed)
        "batch_sorting_strategy": "rerun_first_then_easy",  # From batch config
        "batch_failure_policy": "partial_ok",             # From batch config
    }


def _sum_mopac_times(result: PM7Result) -> Optional[float]:
    """Sum MOPAC execution times across all conformers."""
    times = [c.mopac_execution_time for c in result.conformers if c.mopac_execution_time]
    return round(sum(times), 1) if times else None
```

### Delta-E Calculation Details

From `conformer_selector.py:calculate_delta_e()`:

```python
# Input: List of HOF energies from successful conformers
energies = [c.energy_hof for c in result.conformers if c.is_successful]

# Sorted lowest first
sorted_e = sorted(energies)  # e.g., [-56.78, -55.50, -54.20, -53.00, -51.50]

# E1 = sorted_e[0] = -56.78 (most stable)
# E2 = sorted_e[1] = -55.50
# delta_e_12 = E2 - E1 = -55.50 - (-56.78) = 1.28 kcal/mol

# delta_e_12: E2 - E1 (2nd most stable vs 1st)
# delta_e_13: E3 - E1 (3rd most stable vs 1st)
# delta_e_15: E5 - E1 (5th most stable vs 1st)
```

**Interpretation**:
- `delta_e_12 = 0` → Only 1 conformer succeeded
- `delta_e_12 > 3 kcal/mol` → Significant conformational energy difference
- `delta_e_15 = None` → Fewer than 5 successful conformers

---

## Key Design Decisions

### 1. Fixed Timeouts (NOT Dynamic)
- Timeouts set at batch creation time
- Same timeout for ALL molecules in a batch
- Stored in CSV columns: `crest_timeout_max_minutes`, `mopac_timeout_max_minutes`

### 2. Sequential Processing
- One molecule at a time (no parallelization)
- Simple loop in `execute_batch()`
- KeyboardInterrupt handling to save partial progress

### 3. Retry Logic
- `retry_count` starts at 0
- On failure: `retry_count += 1`
- If `retry_count >= max_retries`: status = Skip
- Else: status = Rerun

### 4. Integration with Existing Code
- Reuse `PM7Result`, `ConformerData` structures
- Reuse existing enums (`CRESTStatus`, `MOPACStatus`, etc.)
- Adapter pattern for fixed timeouts

---

## reset_batch() Implementation (IMPORTANTE #4)

The `reset_batch()` method is only called when `failure_policy == ALL_OR_NOTHING` and the batch completes with failures.

```python
def reset_batch(self, batch_id: str) -> int:
    """Reset all molecules in a batch for re-processing.

    Only applies when failure_policy == ALL_OR_NOTHING.
    Resets successful molecules back to Pending so entire batch can retry.

    Args:
        batch_id: The batch to reset

    Returns:
        Number of molecules reset
    """
    if self.df is None:
        self.load_csv()

    mask = self.df["batch_id"] == batch_id
    batch_rows = self.df[mask]

    if batch_rows.empty:
        return 0

    # Check if this batch uses all_or_nothing policy
    policy = batch_rows["batch_failure_policy"].iloc[0]
    if policy != BatchFailurePolicy.ALL_OR_NOTHING.value:
        LOG.warning(f"reset_batch called on {batch_id} but policy is {policy}")
        return 0

    reset_count = 0
    for idx in batch_rows.index:
        current_status = self.df.at[idx, "status"]

        if current_status == MoleculeStatus.OK.value:
            # Successful molecule → back to Pending (will retry)
            self.df.at[idx, "status"] = MoleculeStatus.PENDING.value
            self.df.at[idx, "batch_id"] = None
            self.df.at[idx, "batch_order"] = None
            reset_count += 1

        elif current_status in [MoleculeStatus.RERUN.value, MoleculeStatus.SKIP.value]:
            # Failed molecule → keep Rerun status but clear batch assignment
            self.df.at[idx, "batch_id"] = None
            self.df.at[idx, "batch_order"] = None
            reset_count += 1

        # Running/To_run molecules: should not exist after batch completion
        # but handle gracefully
        elif current_status in [MoleculeStatus.RUNNING.value, MoleculeStatus.TO_RUN.value]:
            self.df.at[idx, "status"] = MoleculeStatus.PENDING.value
            self.df.at[idx, "batch_id"] = None
            self.df.at[idx, "batch_order"] = None
            reset_count += 1

    self.save_csv()
    return reset_count
```

**When to call**:
```python
# In BatchExecutionManager.execute_batch():
if batch.failure_policy == BatchFailurePolicy.ALL_OR_NOTHING:
    if result.failed_count > 0:
        LOG.warning(f"Batch {batch.batch_id} has failures, resetting all molecules")
        self.csv_manager.reset_batch(batch.batch_id)
```

---

## Retry Timeout Strategy (IMPORTANTE #5)

### Decision: FIXED Timeouts (No Escalation)

Timeouts remain the same on retry:
- **Rationale**: If a molecule timed out at 10 min, it will likely timeout at 10 min again
- **Solution**: User should create a new batch with longer timeouts for problematic molecules

```python
def mark_rerun(self, mol_id: str, error_message: str) -> None:
    """Mark molecule for reprocessing with SAME timeout config."""
    idx = self._get_row_index(mol_id)

    current_retries = self.df.at[idx, "retry_count"]
    max_retries = self.df.at[idx, "max_retries"]

    self.df.at[idx, "retry_count"] = current_retries + 1
    self.df.at[idx, "last_error_message"] = error_message

    # Clear execution results (molecule will be re-processed from scratch)
    self._clear_execution_results(idx)

    if current_retries + 1 >= max_retries:
        self.df.at[idx, "status"] = MoleculeStatus.SKIP.value
        LOG.info(f"{mol_id}: Max retries ({max_retries}) reached → Skip")
    else:
        self.df.at[idx, "status"] = MoleculeStatus.RERUN.value
        LOG.info(f"{mol_id}: Retry {current_retries + 1}/{max_retries}")

    # Timeout config is NOT changed - uses same batch-level timeouts
    # crest_timeout_max_minutes and mopac_timeout_max_minutes remain unchanged

    self.save_csv()


def _clear_execution_results(self, idx: int) -> None:
    """Clear execution results for re-processing."""
    clear_fields = [
        "crest_status", "crest_conformers_generated", "crest_time",
        "crest_error", "crest_output_directory", "crest_ensemble_file",
        "num_conformers_selected", "most_stable_hof", "quality_grade",
        "success", "error_message", "total_execution_time",
        "actual_crest_timeout_used", "actual_mopac_timeout_used",
        "delta_e_12", "delta_e_13", "delta_e_15",
    ]
    for field in clear_fields:
        self.df.at[idx, field] = None
```

### Alternative Strategy (Future Enhancement)

If timeout escalation is needed later:
```python
# NOT IMPLEMENTED NOW - Future consideration
def mark_rerun_with_escalation(self, mol_id: str, error_message: str) -> None:
    """Mark rerun with timeout increase (FUTURE)."""
    if "timeout" in error_message.lower():
        current_crest = self.df.at[idx, "crest_timeout_max_minutes"]
        current_mopac = self.df.at[idx, "mopac_timeout_max_minutes"]
        # Increase by 50%
        self.df.at[idx, "crest_timeout_max_minutes"] = int(current_crest * 1.5)
        self.df.at[idx, "mopac_timeout_max_minutes"] = int(current_mopac * 1.5)
```

---

## mol_id Generation (IMPORTANTE #7)

### Format: `MOL{index:05d}`

```python
def generate_mol_id(index: int) -> str:
    """Generate mol_id from dataset index.

    Args:
        index: Row index from thermo_cbs_clean.csv (1-based from "Unnamed: 0" column)

    Returns:
        mol_id like "MOL00001", "MOL00123", "MOL30026"
    """
    return f"MOL{index:05d}"


# In init_batch_csv.py:
def init_batch_csv(cbs_csv_path: str, output_csv_path: str, ...):
    cbs_df = pd.read_csv(cbs_csv_path)

    rows = []
    for _, row in cbs_df.iterrows():
        # thermo_cbs_clean.csv has "Unnamed: 0" as the original index
        original_index = int(row["Unnamed: 0"])
        mol_id = generate_mol_id(original_index)
        smiles = row["smiles"]
        # ...
```

### Mapping Verification

| thermo_cbs_clean.csv index | mol_id |
|----------------------------|--------|
| 1 | MOL00001 |
| 2 | MOL00002 |
| 8 | MOL00008 |
| 12 | MOL00012 |
| 30026 | MOL30026 |

**Note**: The `Unnamed: 0` column in thermo_cbs_clean.csv is NOT sequential (1, 2, 3, 8, 9, 10, 12...). We preserve the original index for traceability back to the source dataset.

---

## Expanded Test Plan (IMPORTANTE #6)

### Unit Tests: `tests/unit/test_batch/`

#### `test_csv_manager.py`
```python
class TestBatchCSVManager:
    # === Initialization ===
    def test_initialize_from_cbs_dataset_creates_csv(self):
        """Verify CSV is created with correct 36 columns."""

    def test_initialize_generates_correct_mol_ids(self):
        """Verify MOL{index:05d} format for all rows."""

    def test_initialize_calculates_rdkit_descriptors(self):
        """Verify nheavy, nrotbonds, tpsa, aromatic_rings, has_heteroatoms."""

    def test_initialize_sets_default_status_pending(self):
        """All molecules start with status=Pending."""

    # === Batch Selection ===
    def test_select_batch_returns_correct_size(self):
        """Batch contains exactly batch_size molecules."""

    def test_select_batch_prioritizes_rerun_over_pending(self):
        """With rerun_first_then_easy strategy, Rerun comes first."""

    def test_select_batch_sorts_by_nheavy_ascending(self):
        """by_nheavy strategy sorts smallest molecules first."""

    def test_select_batch_random_is_actually_random(self):
        """Random strategy produces different orderings."""

    def test_select_batch_sets_timeout_config(self):
        """crest_timeout_max_minutes and mopac_timeout_max_minutes set correctly."""

    def test_select_batch_when_no_pending_or_rerun(self):
        """Returns empty batch when nothing to process."""

    # === Status Transitions ===
    def test_mark_running_updates_status(self):
        """Pending/To_run → Running."""

    def test_mark_success_updates_all_fields(self):
        """Running → OK with all result fields populated."""

    def test_mark_rerun_increments_retry_count(self):
        """Running → Rerun with retry_count++."""

    def test_mark_rerun_becomes_skip_at_max_retries(self):
        """Running → Skip when retry_count >= max_retries."""

    def test_mark_skip_is_terminal(self):
        """Skip status cannot be changed."""

    # === Edge Cases ===
    def test_select_batch_with_mixed_statuses(self):
        """Correctly handles mix of Pending, Rerun, OK, Skip."""

    def test_concurrent_batch_creation_rejected(self):
        """Cannot create new batch while Running molecules exist."""

    def test_empty_smiles_handling(self):
        """Gracefully handle empty or None SMILES."""

    def test_unicode_smiles_handling(self):
        """Handle SMILES with special characters."""


class TestBatchReset:
    def test_reset_batch_all_or_nothing_policy(self):
        """OK molecules reset to Pending, Rerun stays Rerun."""

    def test_reset_batch_partial_ok_policy_is_noop(self):
        """reset_batch does nothing for partial_ok policy."""

    def test_reset_batch_clears_batch_assignment(self):
        """batch_id and batch_order set to None after reset."""
```

#### `test_detail_manager.py`
```python
class TestConformerDetailManager:
    def test_create_detail_file_success(self):
        """JSON file created with correct structure."""

    def test_create_detail_file_crest_failure(self):
        """JSON has empty conformers list, issues populated."""

    def test_create_detail_file_partial_mopac_success(self):
        """JSON shows which conformers succeeded/failed."""

    def test_load_detail_file_returns_none_for_missing(self):
        """Returns None when file doesn't exist."""

    def test_update_decision_appends(self):
        """Adds decision to existing file without overwriting."""

    def test_delete_detail_file_removes_file(self):
        """File is deleted from disk."""

    # Edge Cases
    def test_create_detail_file_invalid_mol_id(self):
        """Handles mol_id with special characters in filename."""

    def test_load_detail_file_corrupted_json(self):
        """Gracefully handles malformed JSON."""
```

#### `test_execution_manager.py`
```python
class TestBatchExecutionManager:
    def test_execute_batch_sequential_order(self, mock_processor):
        """Molecules processed in batch_order sequence."""

    def test_execute_batch_all_success(self, mock_processor):
        """BatchResult.success_count equals total_count."""

    def test_execute_batch_mixed_results(self, mock_processor):
        """Correct counts for success, rerun, skip."""

    def test_execute_batch_keyboard_interrupt(self, mock_processor):
        """Partial results saved, Running reset to To_run."""

    def test_execute_batch_uses_fixed_timeouts(self, mock_processor):
        """Verifies FixedTimeoutProcessor is called with correct values."""

    # Edge Cases
    def test_execute_empty_batch(self):
        """Empty batch returns zero counts."""

    def test_execute_batch_processor_exception(self, mock_processor):
        """Exception in processor doesn't crash batch, molecule marked Rerun."""

    def test_execute_batch_disk_full(self, mock_processor):
        """Graceful handling when JSON write fails."""
```

#### `test_processor_adapter.py`
```python
class TestFixedTimeoutPredictor:
    def test_predict_always_returns_fixed_value(self):
        """predict() returns same value regardless of nheavy."""

    def test_add_observation_is_noop(self):
        """add_observation doesn't modify internal state."""

    def test_confidence_is_always_high(self):
        """TimeoutConfidence.HIGH returned."""


class TestFixedTimeoutProcessor:
    def test_config_crest_timeout_overridden(self, mock_processor):
        """Verifies config.crest_timeout is set to passed value."""

    def test_mopac_timeout_passed_to_predictor(self, mock_processor):
        """FixedTimeoutPredictor receives correct mopac_timeout."""

    def test_original_config_not_modified(self, original_config):
        """Base config unchanged after processing."""
```

### Integration Tests: `tests/integration/test_batch_pipeline.py`

```python
class TestBatchPipelineIntegration:
    @pytest.fixture
    def small_test_csv(self, tmp_path):
        """Create test CSV with 5 known molecules."""

    def test_full_batch_workflow(self, small_test_csv, mock_crest, mock_mopac):
        """End-to-end: init → select → execute → verify results."""

    def test_batch_resume_after_interrupt(self, small_test_csv):
        """Simulate interrupt, verify resume picks up correctly."""

    def test_multiple_batches_sequential(self, small_test_csv):
        """Run batch 1, then batch 2, verify no overlap."""

    def test_rerun_molecules_picked_up_in_next_batch(self, small_test_csv):
        """Molecules marked Rerun appear in subsequent batch."""
```

### Edge Case Tests

```python
class TestEdgeCases:
    def test_molecule_with_zero_rotatable_bonds(self):
        """Rigid molecule gets 1 conformer."""

    def test_molecule_with_many_rotatable_bonds(self):
        """Flexible molecule (nrotbonds > 4) gets max conformers."""

    def test_molecule_with_no_heavy_atoms(self):
        """Edge case: H2 molecule handling."""

    def test_very_long_smiles(self):
        """SMILES > 500 characters (macrocycle)."""

    def test_batch_with_single_molecule(self):
        """Batch size 1 works correctly."""

    def test_batch_size_larger_than_remaining(self):
        """Batch size 100 but only 5 Pending: returns 5."""

    def test_all_molecules_already_processed(self):
        """No Pending/Rerun left: empty batch returned."""

    def test_csv_with_duplicate_mol_ids(self):
        """Initialization should prevent/detect duplicates."""

    def test_timeout_zero_handling(self):
        """Timeout of 0 should be rejected with clear error."""

    def test_negative_timeout_handling(self):
        """Negative timeout should be rejected."""
```

### Test Fixtures (`tests/fixtures/batch_fixtures.py`)

```python
@pytest.fixture
def sample_pm7result_success():
    """PM7Result with 3 successful conformers."""

@pytest.fixture
def sample_pm7result_partial():
    """PM7Result with 1 success, 2 failures."""

@pytest.fixture
def sample_pm7result_crest_fail():
    """PM7Result with CREST failure (no conformers)."""

@pytest.fixture
def sample_pm7result_all_mopac_fail():
    """PM7Result with CREST success but all MOPAC failed."""

@pytest.fixture
def mock_molecule_processor():
    """Mocked MoleculeProcessor for controlled testing."""
```

### CLI Tests: `tests/integration/test_batch_cli.py`

```python
class TestBatchViewCLI:
    def test_batch_view_renders(self, cli_runner):
        """BatchView shows correctly in menu."""

    def test_run_batch_with_defaults(self, cli_runner):
        """Run batch with default options."""

    def test_run_batch_with_custom_timeouts(self, cli_runner):
        """Run batch with --crest-timeout 15 --mopac-timeout 45."""

    def test_show_status_displays_summary(self, cli_runner):
        """Status shows OK/Rerun/Skip/Pending counts."""

    def test_invalid_timeout_rejected(self, cli_runner):
        """Timeout <= 0 shows error message."""
```

### Verification Checklist

```bash
# 1. Unit tests (85%+ coverage required)
pytest tests/unit/test_batch/ -v --cov=src/grimperium/crest_pm7/batch --cov-report=term-missing
# Target: >85% coverage

# 2. Integration tests
pytest tests/integration/test_batch_pipeline.py -v

# 3. Linting/typing
ruff check src/grimperium/crest_pm7/batch/
mypy src/grimperium/crest_pm7/batch/ --strict
black src/grimperium/crest_pm7/batch/ --check

# 4. Pre-commit
pre-commit run --all-files

# 5. Manual smoke test
python scripts/init_batch_csv.py --input data/thermo_cbs_clean.csv --output data/test_batch.csv --limit 10
# Then run grimperium CLI and execute batch
```

---

## Dependencies

- **pandas**: CSV operations
- **pydantic**: Data validation models
- **rdkit**: Molecular descriptors (already in project)
- **rich**: CLI formatting (already in project)

---

## Estimated Implementation Order

1. `batch/enums.py` (15 min)
2. `batch/models.py` (45 min)
3. `batch/csv_manager.py` (90 min)
4. `batch/detail_manager.py` (30 min)
5. `batch/processor_adapter.py` (45 min)
6. `batch/execution_manager.py` (60 min)
7. `scripts/init_batch_csv.py` (30 min)
8. `cli/views/batch_view.py` (60 min)
9. Tests (90 min)
10. Integration & polish (60 min)

---

## Addition 1: Logging Configuration

**Location**: Add to `batch/execution_manager.py` and all batch managers

**What to add**:

```python
# At top of each file:
import logging

LOG = logging.getLogger(__name__)

# In BatchExecutionManager class:
def __init__(
    self,
    csv_path: Path,
    detail_dir: Path,
    config: PM7Config,
    logger: logging.Logger = None,
):
    self.csv_manager = BatchCSVManager(csv_path, logger=logger)
    self.detail_manager = ConformerDetailManager(detail_dir, logger=logger)
    self.pm7_calculator = pm7_calc
    self.logger = logger or LOG

# In execute_batch() method - add these LOG calls:
def execute_batch(self, batch: Batch) -> BatchResult:
    self.logger.info(f"🚀 Starting batch: {batch.batch_id}")
    self.logger.info(f"   Molecules: {len(batch.molecules)}")
    self.logger.info(f"   CREST timeout: {batch.crest_timeout_minutes} min")
    self.logger.info(f"   MOPAC timeout: {batch.mopac_timeout_minutes} min")

    for i, mol_id in enumerate(batch.molecules, 1):
        self.logger.info(f"[{i}/{len(batch.molecules)}] Processing {mol_id}")

        self.csv_manager.mark_running(mol_id)
        self.logger.debug(f"  Status: Running")

        try:
            result = self.processor_adapter.process_with_fixed_timeout(...)

            if result.success:
                self.csv_manager.mark_success(...)
                self.logger.info(f"  ✅ OK - HOF: {result.most_stable_hof:.2f} kcal/mol")
            else:
                self.csv_manager.mark_rerun(mol_id, result.error_message)
                self.logger.warning(f"  ⚠️  Rerun - {result.error_message}")

        except Exception as e:
            self.logger.error(f"  ❌ Error: {str(e)}", exc_info=True)
            self.csv_manager.mark_skip(mol_id, str(e))

        self.detail_manager.create_detail_file(mol_id, batch.batch_id, result)

    self.logger.info(f"✅ Batch completed: {batch.batch_id}")
    return result

# In handle_interrupt():
def handle_interrupt(self, batch: Batch) -> None:
    self.logger.warning(f"⚠️  Batch interrupted by user: {batch.batch_id}")
    self.logger.info("Saving current state...")
    # ... rest of logic
    self.logger.info("✅ State saved. Resume with: grimperium batch --resume {batch_id}")
```

**Expected output during execution**:
```
🚀 Starting batch: batch_20260113_100000
   Molecules: 10
   CREST timeout: 10 min
   MOPAC timeout: 30 min
[1/10] Processing MOL00001
  Status: Running
  ✅ OK - HOF: -56.78 kcal/mol
[2/10] Processing MOL00002
  Status: Running
  ⚠️  Rerun - CREST timeout after 600s
[3/10] Processing MOL00003
  Status: Running
  ❌ Error: Invalid SMILES geometry
✅ Batch completed: batch_20260113_100000
```

---

## Addition 2: display_batch_summary() Implementation

**Location**: Add to `batch/execution_manager.py` as complete method

**What to add**:

```python
def display_batch_summary(self, result: BatchResult) -> None:
    """Display formatted batch summary with rich formatting.

    Shows:
    - Batch ID and status
    - Success/Rerun/Skip counts and percentages
    - Execution timing (total, average, min, max)
    - Energy statistics (min/max HOF, range)
    - Issues summary
    """
    from rich.table import Table
    from rich.console import Console
    from rich.panel import Panel
    from datetime import timedelta

    console = Console()

    # === Header ===
    header = Panel(
        f"[bold green]✅ BATCH COMPLETED[/bold green]\n"
        f"Batch ID: [cyan]{result.batch_id}[/cyan]\n"
        f"Duration: {timedelta(seconds=int(result.total_time))}",
        title="Batch Summary",
        border_style="green",
    )
    console.print(header)

    # === Results Table ===
    results_table = Table(title="Processing Results")
    results_table.add_column("Status", style="cyan")
    results_table.add_column("Count", justify="right", style="magenta")
    results_table.add_column("Percentage", justify="right", style="magenta")

    total = result.total_count
    results_table.add_row(
        "[green]✅ OK[/green]",
        str(result.success_count),
        f"{100 * result.success_count / total:.1f}%",
    )
    results_table.add_row(
        "[yellow]⚠️  Rerun[/yellow]",
        str(result.rerun_count),
        f"{100 * result.rerun_count / total:.1f}%",
    )
    results_table.add_row(
        "[red]❌ Skip[/red]",
        str(result.skip_count),
        f"{100 * result.skip_count / total:.1f}%",
    )
    results_table.add_row(
        "[blue]⏳ Pending[/blue]",
        str(result.failed_count),
        f"{100 * result.failed_count / total:.1f}%",
    )
    console.print(results_table)

    # === Timing Table ===
    timing_table = Table(title="Execution Timing")
    timing_table.add_column("Metric", style="cyan")
    timing_table.add_column("Value", justify="right", style="magenta")

    avg_time = result.total_time / result.total_count if result.total_count > 0 else 0
    timing_table.add_row("Total Time", f"{result.total_time:.1f} seconds ({timedelta(seconds=int(result.total_time))})")
    timing_table.add_row("Average per Molecule", f"{avg_time:.1f} seconds")
    timing_table.add_row("Molecules/Hour", f"{3600 / avg_time:.1f}" if avg_time > 0 else "N/A")
    timing_table.add_row("Estimated Full Run (30k molecules)", f"{30026 * avg_time / 3600:.1f} hours" if avg_time > 0 else "N/A")
    console.print(timing_table)

    # === Energy Statistics (if available) ===
    if hasattr(result, 'min_hof') and result.min_hof is not None:
        energy_table = Table(title="Energy Statistics")
        energy_table.add_column("Statistic", style="cyan")
        energy_table.add_column("Value (kcal/mol)", justify="right", style="magenta")

        energy_table.add_row("Minimum HOF", f"{result.min_hof:.2f}")
        energy_table.add_row("Maximum HOF", f"{result.max_hof:.2f}")
        energy_table.add_row("Range", f"{result.max_hof - result.min_hof:.2f}")
        console.print(energy_table)

    # === Next Steps ===
    console.print("\n[bold]Next Steps:[/bold]")
    if result.rerun_count > 0:
        console.print(f"  • [yellow]{result.rerun_count} molecules marked for retry[/yellow]")
        console.print(f"    Create new batch with: [cyan]grimperium batch --batch-size {result.rerun_count}[/cyan]")

    if result.skip_count > 0:
        console.print(f"  • [red]{result.skip_count} molecules reached max retries[/red]")
        console.print(f"    Review in CSV: [cyan]data/crest_pm7_batch_status.csv[/cyan]")

    if result.success_count == result.total_count:
        console.print(f"  • [green]All {result.total_count} molecules processed successfully![/green]")
        console.print(f"    Results in: [cyan]data/conformer_details/[/cyan]")

    console.print()  # Blank line
```

**Expected output**:
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃   ✅ BATCH COMPLETED         ┃
┃ Batch ID: batch_20260113_100000
┃ Duration: 0:02:05             ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

Processing Results
┏━━━━━━━┳━━━━━━┳━━━━━━━━━━━━┓
┃ Status ┃ Count ┃ Percentage ┃
┡━━━━━━━╇━━━━━━╇━━━━━━━━━━━━┩
│ ✅ OK  │     8 │      80.0% │
│ ⚠️  Rerun │     1 │      10.0% │
│ ❌ Skip │     1 │      10.0% │
│ ⏳ Pending │     0 │       0.0% │
└────────┴───────┴────────────┘

Execution Timing
┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
┃ Metric               ┃ Value      ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
│ Total Time           │ 125.3s     │
│ Average per Molecule │ 12.5s      │
│ Molecules/Hour       │ 288.0      │
│ Estimated Full Run   │ 103.6 hrs  │
└──────────────────────┴────────────┘

Energy Statistics
┏━━━━━━━━━━━━━┳─────────────────┓
┃ Statistic   ┃ Value (kcal/mol)┃
┡━━━━━━━━━━━━━╇─────────────────┩
│ Minimum HOF │ -145.23         │
│ Maximum HOF │ -32.45          │
│ Range       │ 112.78          │
└─────────────┴─────────────────┘

Next Steps:
  • 1 molecules marked for retry
    Create new batch with: grimperium batch --batch-size 1
  • 1 molecules reached max retries
    Review in CSV: data/crest_pm7_batch_status.csv
```

---

## Addition 3: Parallelization Preparation (Future Enhancement)

### Current Status: Sequential Only
- **Processing model**: One molecule at a time
- **Reason**: Simplicity, easier debugging, predictable resource usage
- **Estimated 30k molecules**: ~200-300 hours (CPU-bound on CREST/MOPAC)

### Architecture is Ready for Parallelization

The current design allows easy upgrade to parallel processing WITHOUT code changes to managers:

**Thread-Safe Components**:
- ✅ **BatchCSVManager**: Uses pandas locks, atomic writes
- ✅ **ConformerDetailManager**: JSON files are per-molecule (no conflicts)
- ✅ **processor_adapter.process_with_fixed_timeout()**: Stateless function
- ✅ **PM7Config**: Immutable per-thread copy

**NOT Thread-Safe (Need Synchronization)**:
- ❌ **CSV DataFrame**: Currently shared, would need locks
  - **Solution**: Use pandas `threading.Lock()` or split CSV per-worker

**Future Implementation Path** (When Needed):

```python
# Phase 1: Add lock to CSV manager
class BatchCSVManager:
    def __init__(self, csv_path: Path):
        self._lock = threading.Lock()

    def mark_success(self, mol_id: str, ...):
        with self._lock:
            # Update CSV atomically
            self.df[...] = ...

# Phase 2: Wrap execute_batch() with ThreadPoolExecutor
from concurrent.futures import ThreadPoolExecutor

def execute_batch_parallel(self, batch: Batch, max_workers: int = 4) -> BatchResult:
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(self._process_molecule, mol_id, batch)
            for mol_id in batch.molecules
        ]
        results = [f.result() for f in futures]
    return self._aggregate_results(results)

def _process_molecule(self, mol_id: str, batch: Batch) -> dict:
    # Single molecule processing (thread-safe)
    self.csv_manager.mark_running(mol_id)
    result = self.processor_adapter.process_with_fixed_timeout(...)
    if result.success:
        self.csv_manager.mark_success(...)  # Locked
    return result

# Phase 3: Add CLI flag for parallel execution
@batch.command()
@click.option("--workers", type=int, default=1, help="Number of parallel workers")
def run_batch(..., workers: int):
    if workers > 1:
        result = batch_mgr.execute_batch_parallel(batch, max_workers=workers)
    else:
        result = batch_mgr.execute_batch(batch)
```

**When to Activate**:
- Estimated need: When 30k molecules becomes bottleneck (~6+ months of production)
- Triggers: Processing time > 200 hours, need faster iteration
- Cost-benefit: 4 workers → ~4x speedup (CREST/MOPAC are CPU-bound)

**Not Activated Now Because**:
- ✅ Sequential execution is simpler to debug
- ✅ Lower risk of race conditions
- ✅ Adequate for current Phase A (~10-20 molecules)
- ✅ Can activate incrementally when needed

---

## Summary of Additions

| Addition | Lines | Impact |
|----------|-------|--------|
| **Logging** | ~40 lines | Production-grade monitoring |
| **display_batch_summary()** | ~80 lines | User-friendly reporting |
| **Parallelization note** | ~60 lines | Future-proofing & clarity |
| **Total** | ~180 lines | Plan → 100% Complete ✅ |

---

## What We're NOT Adding (Intentionally)

- ❌ **Logging config file** (can use Python defaults)
- ❌ **Rich library styling** (already in project dependencies)
- ❌ **Actual parallelization** (future, not needed now)
- ❌ **Async/await** (overkill, threading sufficient)

---

## Final Checklist After These Additions

- [ ] Logging added to all batch managers
- [ ] display_batch_summary() produces rich-formatted output
- [ ] Parallelization note explains future path without blocking current implementation
- [ ] No breaking changes to existing architecture
- [ ] All 3 additions are backward-compatible
- [ ] Ready for Claude Code implementation

---

## Confirmed Decisions

| Decision | Choice |
|----------|--------|
| CLI Style | New `BatchView` (separate from CalcView) |
| mol_id Format | `MOL{index:05d}` (e.g., MOL00001, MOL00002) |
| Default CREST Timeout | 10 minutes |
| Default MOPAC Timeout | 30 minutes |
| JSON Detail Location | `data/conformer_details/` |

---

## Ready for Implementation

All design decisions confirmed. Implementation can proceed following the steps outlined above.
