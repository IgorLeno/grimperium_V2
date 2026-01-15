# Grimperium Phase A: Workflow Documentation

> **Version:** 1.0.0
> **Last Updated:** January 2026
> **Target:** CREST PM7 batch processing for ~30,000 molecules

---

## Table of Contents

1. [Overview](#overview)
2. [Single Molecule Workflow](#single-molecule-workflow)
3. [Batch Processing Architecture](#batch-processing-architecture)
4. [Data Flow and State Management](#data-flow-and-state-management)
5. [Timeout Management](#timeout-management)
6. [Parallelization Strategy](#parallelization-strategy)
7. [Monitoring and Status Tracking](#monitoring-and-status-tracking)
8. [Quality Assurance](#quality-assurance)
9. [Error Handling and Recovery](#error-handling-and-recovery)
10. [Performance Expectations](#performance-expectations)
11. [CLI Interface](#cli-interface)
12. [Output Artifacts](#output-artifacts)

---

## Overview

Grimperium Phase A processes molecules through a CREST + PM7 pipeline to generate
semi-empirical Heat of Formation (HOF) values. These values serve as baseline
predictions that will be corrected using delta-learning ML models trained on
CBS (Complete Basis Set) reference data.

### Pipeline Goal

```
SMILES → 3D Geometry → CREST Conformers → PM7 Optimization → HOF (kcal/mol)
```

### Key Components

| Component | Purpose |
|-----------|---------|
| `MoleculeProcessor` | Single molecule orchestration |
| `BatchExecutionManager` | Batch execution control |
| `BatchCSVManager` | State persistence and batch selection |
| `FixedTimeoutProcessor` | Timeout management adapter |
| `PM7Config` | Configuration parameters |

---

## Single Molecule Workflow

### Process Overview

The `MoleculeProcessor.process()` method orchestrates the complete pipeline for
a single molecule:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          SINGLE MOLECULE FLOW                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  1. INPUT                                                                   │
│     ├── mol_id (string identifier)                                         │
│     ├── SMILES (molecular structure)                                       │
│     └── input_xyz (optional, otherwise generated)                          │
│                                                                             │
│  2. DESCRIPTOR CALCULATION                                                  │
│     ├── nheavy (non-hydrogen atoms)                                        │
│     ├── nrotbonds (rotatable bonds)                                        │
│     ├── tpsa (topological polar surface area)                              │
│     ├── aromatic_rings (count)                                             │
│     └── has_heteroatoms (N, O, S, etc.)                                    │
│                                                                             │
│  3. 3D COORDINATE GENERATION (if no input_xyz)                             │
│     ├── RDKit MolFromSmiles()                                              │
│     ├── AddHs() for explicit hydrogens                                     │
│     ├── EmbedMolecule() with ETKDGv3                                       │
│     └── MMFF94 geometry optimization                                       │
│                                                                             │
│  4. CONFORMER COUNT DECISION                                                │
│     └── get_num_conformers(nrotbonds, config)                              │
│         ├── 0-2 bonds → 3 conformers                                       │
│         ├── 3-5 bonds → 5 conformers                                       │
│         ├── 6-9 bonds → 7 conformers                                       │
│         └── 10+ bonds → 10 conformers                                      │
│                                                                             │
│  5. TIMEOUT PREDICTION                                                      │
│     └── timeout_predictor.predict(nheavy, num_conformers)                  │
│         ├── Fixed timeout from config                                      │
│         └── Confidence level (HIGH/MEDIUM/LOW)                             │
│                                                                             │
│  6. CREST EXECUTION                                                         │
│     ├── run_crest(mol_id, input_xyz, config, smiles)                       │
│     ├── iMTD-GC conformer search                                           │
│     └── Returns: conformer XYZ files                                       │
│                                                                             │
│  7. PM7 OPTIMIZATION (per conformer)                                        │
│     ├── Dynamic timeout redistribution                                      │
│     ├── optimize_conformer() for each conformer                            │
│     ├── HOF extraction from .out file                                      │
│     └── Energy ranking                                                      │
│                                                                             │
│  8. ENERGY ANALYSIS                                                         │
│     ├── Select lowest HOF conformer                                         │
│     ├── Calculate delta_e_12, delta_e_13, delta_e_15                       │
│     └── Energy gaps for quality assessment                                  │
│                                                                             │
│  9. QUALITY GRADING                                                         │
│     └── _assign_quality_grade(result)                                       │
│         ├── A: Excellent (all conformers, low delta_e)                     │
│         ├── B: Good (most conformers successful)                           │
│         ├── C: Acceptable (>50% conformers)                                │
│         ├── D: Poor (single conformer only)                                │
│         └── F: Failed (no valid HOF)                                       │
│                                                                             │
│ 10. OUTPUT                                                                  │
│     └── PM7Result with all data                                            │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Stage Details

#### Stage 1: Input Validation

```python
# Input parameters
mol_id: str      # e.g., "mol_00001"
smiles: str      # e.g., "CCO" (ethanol)
input_xyz: Path  # Optional pre-computed geometry
```

#### Stage 2: Molecular Descriptors

Descriptors drive timeout prediction and quality assessment:

```python
def compute_molecular_descriptors(smiles: str) -> dict:
    """Calculate key molecular properties."""
    return {
        "nheavy": 10,           # Heavy atoms (non-H)
        "nrotbonds": 3,         # Rotatable bonds
        "tpsa": 45.2,           # Polar surface area
        "aromatic_rings": 1,    # Aromatic ring count
        "has_heteroatoms": True # N, O, S, etc.
    }
```

#### Stage 3: 3D Coordinate Generation

RDKit generates initial 3D geometry:

```python
def _generate_xyz_from_smiles(mol_id: str, smiles: str) -> Path | None:
    """Generate 3D coordinates using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # ETKDGv3 embedding with MMFF optimization
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    success = AllChem.EmbedMolecule(mol, params)

    if success == 0:
        AllChem.MMFFOptimizeMolecule(mol)
        # Write to XYZ file
        return xyz_path
    return None
```

#### Stage 4: Conformer Selection

Number of conformers based on molecular flexibility:

| Rotatable Bonds | Conformers | Rationale |
|-----------------|------------|-----------|
| 0-2 | 3 | Rigid molecules |
| 3-5 | 5 | Moderate flexibility |
| 6-9 | 7 | Flexible molecules |
| 10+ | 10 | Highly flexible |

#### Stage 5: Timeout Prediction

Fixed timeout strategy for Phase A batch processing:

```python
@dataclass
class FixedTimeoutPredictor:
    """Predict molecule processing timeout."""
    crest_timeout_seconds: float = 1800  # 30 min default
    mopac_timeout_seconds: float = 3600  # 60 min default
```

#### Stage 6: CREST Execution

See [CREST_INTEGRATION.md](CREST_INTEGRATION.md) for details.

#### Stage 7: PM7 Optimization

See [MOPAC_INTEGRATION.md](MOPAC_INTEGRATION.md) for details.

#### Stage 8: Quality Grading

```python
class QualityGrade(Enum):
    """Quality grades for PM7 results."""
    A = "A"  # Excellent: All conformers, low delta_e
    B = "B"  # Good: Most conformers successful
    C = "C"  # Acceptable: >50% conformers
    D = "D"  # Poor: Single conformer
    F = "F"  # Failed: No valid HOF
```

---

## Batch Processing Architecture

### Component Hierarchy

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         BATCH PROCESSING LAYERS                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Layer 4: CLI Interface                                                     │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  BatchCalculatorView (TUI)                                          │   │
│  │  - User configuration interface                                      │   │
│  │  - Progress display and monitoring                                   │   │
│  │  - Error reporting                                                   │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│  Layer 3: Execution Management                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  BatchExecutionManager                                              │   │
│  │  - Orchestrates batch execution                                      │   │
│  │  - Handles failure policies                                          │   │
│  │  - Coordinates status updates                                        │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│  Layer 2: Processing Adapters                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  FixedTimeoutProcessor                                              │   │
│  │  - Adapts MoleculeProcessor for batch mode                          │   │
│  │  - Manages timeout configuration                                     │   │
│  │  - Fixed timeout strategy for Phase A                                │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│  Layer 1: Core Processing                                                   │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  MoleculeProcessor                                                  │   │
│  │  - Single molecule processing                                        │   │
│  │  - CREST + PM7 execution                                             │   │
│  │  - Result aggregation                                                │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### BatchExecutionManager

Central coordinator for batch processing:

```python
class BatchExecutionManager:
    """Manage batch execution of molecules."""

    def __init__(
        self,
        csv_manager: BatchCSVManager,
        processor_adapter: FixedTimeoutProcessor,
        detail_dir: Path,
    ):
        self.csv_manager = csv_manager
        self.processor_adapter = processor_adapter
        self.detail_dir = detail_dir

    def execute_batch(
        self,
        batch: Batch,
        progress_callback: Callable | None = None,
    ) -> BatchResult:
        """Execute a batch of molecules."""
        # 1. Update timeout configuration
        self.processor_adapter.update_timeouts(
            crest_timeout_minutes=batch.crest_timeout_minutes,
            mopac_timeout_minutes=batch.mopac_timeout_minutes,
        )

        # 2. Process each molecule
        for mol in batch.molecules:
            self._process_molecule(mol, batch, result)

        # 3. Handle failure policy
        if batch.failure_policy == BatchFailurePolicy.ALL_OR_NOTHING:
            if result.failed_count > 0:
                self.csv_manager.reset_batch(batch.batch_id)

        return result
```

### BatchCSVManager

State persistence and batch selection:

```python
class BatchCSVManager:
    """Manage CSV state for batch processing."""

    # Status values
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    RERUN = "rerun"
    SKIP = "skip"

    def select_batch(
        self,
        batch_id: str,
        batch_size: int,
        sorting_strategy: SortingStrategy,
    ) -> list[dict]:
        """Select molecules for a new batch."""
        # 1. Filter pending molecules
        pending = self.df[self.df["status"] == self.PENDING]

        # 2. Apply sorting strategy
        sorted_df = self._apply_sorting_strategy(pending, sorting_strategy)

        # 3. Select batch_size molecules
        selected = sorted_df.head(batch_size)

        # 4. Mark as running
        for mol_id in selected["mol_id"]:
            self.mark_running(mol_id, batch_id)

        return selected.to_dict("records")
```

### Batch Data Model

```python
@dataclass
class Batch:
    """Represents a batch of molecules to process."""
    batch_id: str
    molecules: list[BatchMolecule]
    crest_timeout_minutes: float = 30.0
    mopac_timeout_minutes: float = 60.0
    failure_policy: BatchFailurePolicy = BatchFailurePolicy.CONTINUE

@dataclass
class BatchMolecule:
    """A molecule in a batch."""
    mol_id: str
    smiles: str
    batch_order: int

@dataclass
class BatchResult:
    """Results from batch execution."""
    batch_id: str
    total_count: int
    success_count: int = 0
    failed_count: int = 0
    rerun_count: int = 0
    skip_count: int = 0
    total_time: float = 0.0
    min_hof: float | None = None
    max_hof: float | None = None
```

---

## Data Flow and State Management

### CSV State Machine

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          STATUS STATE MACHINE                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│                              ┌──────────┐                                   │
│                              │ PENDING  │◄─────────────────────────┐        │
│                              └────┬─────┘                          │        │
│                                   │                                │        │
│                          select_batch()                            │        │
│                                   │                            reset_batch()│
│                                   ▼                                │        │
│                              ┌──────────┐                          │        │
│                              │ RUNNING  │──────────────────────────┤        │
│                              └────┬─────┘                          │        │
│                                   │                                │        │
│                    ┌──────────────┼──────────────┐                 │        │
│                    │              │              │                 │        │
│              mark_success()  mark_failed()  mark_rerun()           │        │
│                    │              │              │                 │        │
│                    ▼              ▼              ▼                 │        │
│              ┌──────────┐  ┌──────────┐  ┌──────────┐              │        │
│              │ SUCCESS  │  │  FAILED  │  │  RERUN   │──────────────┘        │
│              └──────────┘  └──────────┘  └──────────┘                       │
│                                                                             │
│                          ┌──────────┐                                       │
│                          │   SKIP   │  (manual exclusion)                   │
│                          └──────────┘                                       │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### CSV Schema

```csv
mol_id,smiles,status,batch_id,batch_order,timestamp,success,h298_pm7,quality_grade,conformers_found,error_message
mol_00001,CCO,success,batch_001,1,2026-01-15T10:30:00,-45.23,A,5,
mol_00002,CC(=O)O,pending,,,,,,,,
mol_00003,c1ccccc1,failed,batch_001,2,2026-01-15T10:35:00,,F,0,CREST timeout
```

### Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| `mol_id` | str | Unique molecule identifier |
| `smiles` | str | SMILES string |
| `status` | str | Current status (pending/running/success/failed/rerun/skip) |
| `batch_id` | str | Batch identifier if assigned |
| `batch_order` | int | Order within batch |
| `timestamp` | str | ISO timestamp of last update |
| `success` | bool | Whether processing succeeded |
| `h298_pm7` | float | PM7 HOF result (kcal/mol) |
| `quality_grade` | str | Quality grade (A-F) |
| `conformers_found` | int | Number of CREST conformers |
| `error_message` | str | Error details if failed |

---

## Timeout Management

### Fixed Timeout Strategy (Phase A)

For reliable batch processing of 30k molecules, Phase A uses fixed timeouts:

```python
@dataclass
class FixedTimeoutProcessor:
    """Fixed timeout processor for batch mode."""

    def __init__(
        self,
        config: PM7Config,
        crest_timeout_minutes: float = 30.0,
        mopac_timeout_minutes: float = 60.0,
    ):
        self.config = config
        self._timeout_predictor = FixedTimeoutPredictor(
            crest_timeout_seconds=crest_timeout_minutes * 60,
            mopac_timeout_seconds=mopac_timeout_minutes * 60,
        )
        # CRITICAL: Update config.crest_timeout for subprocess
        self.config.crest_timeout = crest_timeout_minutes * 60
```

### Timeout Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          TIMEOUT FLOW                                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  User Configuration (CLI)                                                   │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  crest_timeout_minutes = 30                                         │   │
│  │  mopac_timeout_minutes = 60                                         │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│  FixedTimeoutProcessor.update_timeouts()                                    │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  self._timeout_predictor.crest_timeout_seconds = 30 * 60 = 1800     │   │
│  │  self._timeout_predictor.mopac_timeout_seconds = 60 * 60 = 3600     │   │
│  │  self.config.crest_timeout = 30 * 60 = 1800  # CRITICAL!            │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│  MoleculeProcessor.process()                                                │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  timeout = self.timeout_predictor.predict(nheavy, conformers)       │   │
│  │  → Returns: (1800, Confidence.HIGH)                                 │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│  run_crest() / optimize_conformer()                                         │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  subprocess.run(..., timeout=config.crest_timeout)                  │   │
│  │  subprocess.run(..., timeout=per_conformer_timeout)                 │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Recommended Timeouts

| Molecule Size | CREST Timeout | MOPAC Timeout | Total Budget |
|---------------|---------------|---------------|--------------|
| Small (<20 heavy) | 15 min | 30 min | ~45 min |
| Medium (20-40 heavy) | 30 min | 60 min | ~90 min |
| Large (40+ heavy) | 60 min | 120 min | ~180 min |

---

## Parallelization Strategy

### Phase A: Sequential Processing

For Phase A initial implementation, molecules are processed sequentially within
each batch. This ensures:

1. **Simplicity**: Easier debugging and error tracking
2. **Reproducibility**: Deterministic execution order
3. **Resource Control**: No memory/CPU contention

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     PHASE A: SEQUENTIAL PROCESSING                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Batch 001: [mol_1, mol_2, mol_3, ..., mol_100]                            │
│                                                                             │
│  Processing Timeline:                                                       │
│  ────────────────────────────────────────────────────────────────────►      │
│                                                                             │
│  │ mol_1  │ mol_2  │ mol_3  │ ... │ mol_100 │                              │
│  │  CREST │  CREST │  CREST │     │  CREST  │                              │
│  │  PM7   │  PM7   │  PM7   │     │  PM7    │                              │
│                                                                             │
│  Estimated time for 100 molecules:                                          │
│  - Average: 5-10 min/molecule                                               │
│  - Total: 8-17 hours/batch                                                  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Future: Parallel Processing (Phase B+)

Future phases may implement parallel processing:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     FUTURE: PARALLEL PROCESSING                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Worker Pool (4 workers)                                                    │
│                                                                             │
│  Worker 1: │ mol_1  │ mol_5  │ mol_9   │ ...                               │
│  Worker 2: │ mol_2  │ mol_6  │ mol_10  │ ...                               │
│  Worker 3: │ mol_3  │ mol_7  │ mol_11  │ ...                               │
│  Worker 4: │ mol_4  │ mol_8  │ mol_12  │ ...                               │
│                                                                             │
│  Speedup: ~4x (limited by I/O and memory)                                   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Monitoring and Status Tracking

### Progress Callback

```python
def progress_callback(mol_id: str, current: int, total: int) -> None:
    """Called after each molecule is processed."""
    percentage = (current / total) * 100
    print(f"[{current}/{total}] {percentage:.1f}% - Processing {mol_id}")
```

### Status Summary

```python
def get_status_summary(csv_manager: BatchCSVManager) -> dict:
    """Get current processing status."""
    counts = csv_manager.get_status_counts()
    return {
        "pending": counts.get("pending", 0),
        "running": counts.get("running", 0),
        "success": counts.get("success", 0),
        "failed": counts.get("failed", 0),
        "rerun": counts.get("rerun", 0),
        "skip": counts.get("skip", 0),
        "total": sum(counts.values()),
        "completion_rate": counts.get("success", 0) / max(1, sum(counts.values())),
    }
```

### Batch Result Statistics

```python
@dataclass
class BatchResult:
    """Comprehensive batch results."""
    batch_id: str
    total_count: int
    success_count: int = 0
    failed_count: int = 0
    rerun_count: int = 0
    skip_count: int = 0
    total_time: float = 0.0

    # HOF statistics
    min_hof: float | None = None
    min_hof_mol_id: str | None = None
    max_hof: float | None = None
    max_hof_mol_id: str | None = None

    # Timestamps
    timestamp_start: datetime | None = None
    timestamp_end: datetime | None = None

    @property
    def success_rate(self) -> float:
        """Calculate success rate."""
        return self.success_count / max(1, self.total_count)

    @property
    def avg_time_per_molecule(self) -> float:
        """Average processing time per molecule."""
        return self.total_time / max(1, self.total_count)
```

---

## Quality Assurance

### Quality Grade System

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          QUALITY GRADING SYSTEM                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Grade A: EXCELLENT                                                         │
│  ├── All conformers optimized successfully                                  │
│  ├── HOF extracted with high confidence                                     │
│  ├── Energy gaps within expected ranges                                     │
│  └── Example: delta_e_12 < 2 kcal/mol                                       │
│                                                                             │
│  Grade B: GOOD                                                              │
│  ├── >80% conformers successful                                             │
│  ├── HOF extracted reliably                                                 │
│  └── Minor issues (e.g., one conformer failed)                              │
│                                                                             │
│  Grade C: ACCEPTABLE                                                        │
│  ├── >50% conformers successful                                             │
│  ├── HOF available but with caveats                                         │
│  └── Consider for training with caution                                     │
│                                                                             │
│  Grade D: POOR                                                              │
│  ├── Only 1-2 conformers successful                                         │
│  ├── HOF may not represent global minimum                                   │
│  └── Flag for manual review                                                 │
│                                                                             │
│  Grade F: FAILED                                                            │
│  ├── No valid HOF obtained                                                  │
│  ├── CREST or MOPAC completely failed                                       │
│  └── Requires investigation or exclusion                                    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Expected Grade Distribution (Target)

| Grade | Expected % | Description |
|-------|------------|-------------|
| A | 70-80% | Production-ready data |
| B | 10-15% | Usable with minor caveats |
| C | 5-10% | Review before use |
| D | 2-5% | Likely problematic |
| F | <2% | Processing failure |

---

## Error Handling and Recovery

### Error Categories

```python
class ErrorCategory(Enum):
    """Categories of processing errors."""
    SMILES_INVALID = "invalid_smiles"
    GEOMETRY_FAILED = "geometry_generation_failed"
    CREST_TIMEOUT = "crest_timeout"
    CREST_ERROR = "crest_execution_error"
    MOPAC_TIMEOUT = "mopac_timeout"
    MOPAC_ERROR = "mopac_execution_error"
    HOF_EXTRACTION = "hof_extraction_failed"
    UNEXPECTED = "unexpected_error"
```

### Recovery Strategies

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          ERROR RECOVERY STRATEGIES                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Error Type              │ Recovery Strategy                                │
│  ────────────────────────┼──────────────────────────────────────────────   │
│  SMILES_INVALID          │ Mark skip, report for dataset cleaning          │
│  GEOMETRY_FAILED         │ Mark rerun, try alternative embedding           │
│  CREST_TIMEOUT           │ Mark rerun with higher timeout                  │
│  CREST_ERROR             │ Analyze output, may need parameter adjustment   │
│  MOPAC_TIMEOUT           │ Mark rerun with higher timeout                  │
│  MOPAC_ERROR             │ Check SCF convergence, adjust keywords          │
│  HOF_EXTRACTION          │ Rerun with verbose output                       │
│  UNEXPECTED              │ Log full traceback, mark failed                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Failure Policies

```python
class BatchFailurePolicy(Enum):
    """How to handle failures within a batch."""

    CONTINUE = "continue"
    # Process all molecules, mark failures individually
    # Use for: Initial processing, maximizing throughput

    ALL_OR_NOTHING = "all_or_nothing"
    # Reset entire batch if any failure occurs
    # Use for: Critical batches requiring 100% success

    STOP_ON_FAILURE = "stop_on_failure"
    # Stop batch on first failure
    # Use for: Debugging, testing new configurations
```

---

## Performance Expectations

### Per-Molecule Timing

| Stage | Small Mol | Medium Mol | Large Mol |
|-------|-----------|------------|-----------|
| Descriptor calc | <1s | <1s | <2s |
| 3D generation | 1-5s | 5-15s | 15-60s |
| CREST | 1-10 min | 5-30 min | 20-60 min |
| PM7 (all conf) | 1-5 min | 5-20 min | 15-60 min |
| **Total** | **3-20 min** | **15-60 min** | **45-120 min** |

### Batch Processing Estimates

For 30,000 molecules with CHON-only dataset:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     BATCH PROCESSING ESTIMATES                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Configuration:                                                             │
│  - Dataset: CBS Reference (CHON-only), 29,568 molecules                    │
│  - Batch size: 100 molecules                                                │
│  - CREST timeout: 30 min                                                    │
│  - MOPAC timeout: 60 min                                                    │
│                                                                             │
│  Optimistic Estimate (small molecules, no failures):                        │
│  - Average: 5 min/molecule                                                  │
│  - Per batch: 8.3 hours                                                     │
│  - Total: ~296 batches × 8.3 hours = 102 days                              │
│                                                                             │
│  Realistic Estimate (mixed sizes, ~5% failures):                            │
│  - Average: 10 min/molecule                                                 │
│  - Per batch: 16.7 hours                                                    │
│  - Total: ~296 batches × 16.7 hours = 205 days                             │
│                                                                             │
│  Conservative Estimate (larger molecules, ~10% failures):                   │
│  - Average: 20 min/molecule                                                 │
│  - Per batch: 33 hours                                                      │
│  - Total: ~296 batches × 33 hours = 407 days                               │
│                                                                             │
│  Parallel Processing Speedup:                                               │
│  - 4 workers: ~25-100 days                                                  │
│  - 8 workers: ~13-50 days                                                   │
│  - Cloud scaling: ~1-2 weeks                                                │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Success Rate Targets

| Metric | Target | Acceptable | Concern |
|--------|--------|------------|---------|
| Overall Success | >98% | >95% | <90% |
| Grade A+B | >85% | >75% | <60% |
| Timeout Rate | <5% | <10% | >15% |
| HOF Extraction | >99% | >97% | <95% |

---

## CLI Interface

### Interactive Configuration

The CLI provides an interactive TUI for batch configuration:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     CREST PM7 BATCH CALCULATOR                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Dataset Configuration                                                      │
│  ├── Dataset: CBS Reference (CHON-only) [29,568 molecules]                 │
│  ├── CSV Status: 1,234 pending, 567 success, 23 failed                     │
│  └── Sorting: By heavy atoms (ascending)                                    │
│                                                                             │
│  Batch Settings                                                             │
│  ├── Batch ID: batch_20260115_001                                          │
│  ├── Batch Size: 100 molecules                                              │
│  ├── CREST Timeout: 30 minutes                                              │
│  └── MOPAC Timeout: 60 minutes                                              │
│                                                                             │
│  Actions                                                                    │
│  ├── [Start Batch] - Begin processing                                       │
│  ├── [Preview Batch] - Show selected molecules                              │
│  ├── [View Status] - Show current status                                    │
│  └── [Back] - Return to main menu                                           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### CLI Commands

```bash
# Start interactive TUI
grimperium

# Quick batch start (future)
grimperium batch start --size 100 --timeout 30

# Check status (future)
grimperium batch status

# Export results (future)
grimperium batch export --format json
```

---

## Output Artifacts

### Directory Structure

```
data/molecules_pm7/
├── computed/
│   ├── phase_a_results.json      # Aggregated results
│   └── batch_001/
│       ├── batch_summary.json    # Batch statistics
│       ├── mol_00001.json        # Detailed result
│       ├── mol_00002.json
│       └── ...
├── crest_output/
│   ├── mol_00001/
│   │   ├── crest_best.xyz        # Best conformer
│   │   ├── crest_conformers.xyz  # All conformers
│   │   └── crest.output          # CREST log
│   └── ...
├── mopac_output/
│   ├── mol_00001/
│   │   ├── conf_0.mop            # MOPAC input
│   │   ├── conf_0.out            # MOPAC output
│   │   └── conf_0.arc            # Archive file
│   └── ...
└── status.csv                     # Master status file
```

### Result JSON Format

```json
{
  "results": [
    {
      "mol_id": "mol_00001",
      "smiles": "CCO",
      "success": true,
      "h298_pm7": -45.23,
      "quality_grade": "A",
      "conformers_found": 5,
      "conformers_successful": 5,
      "delta_e_12": 0.5,
      "delta_e_13": 1.2,
      "execution_time": 324.5,
      "timestamp": "2026-01-15T10:30:00Z"
    }
  ],
  "metadata": {
    "batch_id": "batch_001",
    "total_molecules": 100,
    "success_count": 98,
    "failed_count": 2,
    "total_time": 32450.0,
    "timestamp": "2026-01-15T10:30:00Z"
  }
}
```

---

## Appendix: Quick Reference

### Key Files

| File | Purpose |
|------|---------|
| `molecule_processor.py` | Single molecule orchestration |
| `execution_manager.py` | Batch execution control |
| `csv_manager.py` | State persistence |
| `processor_adapter.py` | Timeout management |
| `crest_runner.py` | CREST execution |
| `mopac_optimizer.py` | PM7 optimization |

### Configuration Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `batch_size` | 100 | 10-500 | Molecules per batch |
| `crest_timeout` | 30 min | 15-120 min | CREST time limit |
| `mopac_timeout` | 60 min | 30-180 min | MOPAC time limit |
| `sorting_strategy` | `by_heavy_atoms` | - | Batch selection order |
| `failure_policy` | `continue` | - | Error handling mode |

### Status Transitions

```
pending → running → success
                 → failed
                 → rerun → pending
                 → skip
```

---

## See Also

- [CREST_INTEGRATION.md](CREST_INTEGRATION.md) - CREST configuration and usage
- [MOPAC_INTEGRATION.md](MOPAC_INTEGRATION.md) - MOPAC/PM7 details
- [README.md](../README.md) - Project overview
