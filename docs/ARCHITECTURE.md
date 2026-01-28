# Grimperium Architecture

**Document Version:** 2.2.0
**Last Updated:** 2026-01-28
**Next Review:** 2026-06-19 (or after 100K molecules processed)
**Version:** v2.2
**Status:** Production-Ready (Single-Process)

## Table of Contents

1. [System Overview](#system-overview)
2. [Known Limitations](#known-limitations)
3. [Data Flow](#data-flow)
4. [Future Considerations](#future-considerations)
5. [Decision Log](#decision-log)

---

## System Overview

### High-Level Architecture

```
┌─────────────────────────────────────────────────────┐
│                   GRIMPERIUM v2.2                   │
├─────────────────────────────────────────────────────┤
│                                                     │
│  CLI Layer (main.py)                                │
│  ↓                                                  │
│  BatchOrchestrator (coordinator)                    │
│  ├─ BatchDataManager (load CSV)                     │
│  ├─ BatchScheduler (prioritize)                     │
│  ├─ CalculationExecutor (CREST + MOPAC)             │
│  ├─ MoleculePersister (atomic writes)               │
│  └─ BatchReporter (summaries)                       │
│  ↓                                                  │
│  CSV Backend (thermo_pm7.csv)                       │
│  ↓                                                  │
│  ~/.grimperium/config.toml (settings)               │
│  ~/.grimperium/validation_errors.log (audit)        │
│  ~/.grimperium/batch_summary.log (history)          │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### Technology Stack

| Component | Technology | Rationale |
|-----------|-----------|-----------|
| **Data Format** | CSV | Simple, auditable, versionable in git |
| **Settings** | TOML | Human-readable, standardized |
| **In-Memory** | Pandas DataFrame → Molecule dataclasses | Type-safe, clear semantics |
| **Persistence** | Atomic file writes (temp + replace) | Crash-safe, no corruption |
| **Process Model** | Single-threaded, sequential | Simple, deterministic, safe |
| **Logging** | Python logging | Structured, filterable |

### CLI Batch Progress Tracking

The CLI batch view renders a CSV-driven, 5-stage progress bar using a
daemon CSV poller and a thread-safe Queue. The tracker reads only the CSV
and never touches processing code directly.

Stages (CSV transitions):
1. `status`: `Pending`/`Selected` → `Running` (RDKit parameters)
2. `crest_status`: `NOT_ATTEMPTED` → `XTB_PREOPT` (pre-optimization; xTB if enabled)
3. `crest_status`: `XTB_PREOPT`/`NOT_ATTEMPTED` → `CREST_SEARCH` (CREST search)
4. `mopac_status`: `NOT_ATTEMPTED` → `RUNNING` (MOPAC PM7 calculation)
5. `status`: `Running` → `OK` (final calculations)

Completion counters update when `status` reaches `OK`, `Rerun`, or `Skip`.
Final result values still use `crest_status` = `SUCCESS`/`FAILED` and
`mopac_status` = `OK`/`FAILED`.
The UI renders a 30-character bar (6 chars per stage) with a live spinner.

### Key Design Decisions

#### Decision 1: CSV as Source of Truth
- **Why:** Simple, auditable, versionable, human-readable
- **Tradeoff:** No concurrent access, no transactions
- **Constraint:** Must remain single-process

#### Decision 2: Atomic 1-Write-Per-Molecule
- **Why:** Ensures CSV never corrupted, even if process crashes
- **How:** Write to temp file, atomic rename (os.replace)
- **Guarantee:** Read after crash = always valid state

#### Decision 3: Type-Safe Dataclass Hierarchy
- **Why:** Catch errors at parse time, not runtime
- **Structure:** MoleculeIdentity + Properties + Results + Meta
- **Benefit:** IDE hints, mypy validation, clear schema

#### Decision 4: Settings as Typed Objects
- **Why:** Avoid dict access bugs (typos, wrong types)
- **Persistence:** TOML file, auto-loaded
- **Thread-Safe:** Immutable per-session (no concurrent updates)

#### Decision 5: Visibility First (Permissive Mode)
- **Why:** Recovery is good, but visibility is better
- **How:** Always print validation warnings + persistent logs
- **Monitoring:** grep logs for error patterns

---

## Known Limitations

### Limitation 1: CSV Backend is NOT Transactional

**Current State:**
```
✅ Atomic writes (1 write/molecule)
✅ Single-process safety
✅ Crash-safe persistence
❌ No file-level locking
❌ No rollback capability
❌ No concurrent access
```

**Why This Works Now:**
```
Single-process execution
  → Only one writer at a time
  
Linear processing order
  → Deterministic, no race conditions
  
Atomic writes (temp + replace)
  → Even if crash mid-write, previous state intact
```

**Why This Will Break:**
```python
# DO NOT DO THIS (will corrupt data):
import multiprocessing

def process_mol(mol):
    mol.results.delta_1 = calculate(mol)
    persister.save_molecule(mol)

pool = multiprocessing.Pool(4)
pool.map(process_mol, molecules)  # ← BREAKS: no file locking
# Result: Molecules overwrite each other randomly
```

**What Happens:**
```
Process 1 reads CSV
Process 2 reads CSV (same version)
Process 1 updates mol_001, writes CSV
Process 2 updates mol_002, writes CSV  # ← Overwrites Process 1's work!
Result: mol_001 lost
```

**Solution Strategy (for future, NOT NOW):**

**Option A: SQLite (Recommended for ≤1M molecules)**
```python
# Future (v3.0+):
from grimperium.backends import SQLitePersister

persister = SQLitePersister("grimperium.db")
# Benefits: ACID, indexes, query support, concurrent reads
```

**Option B: PostgreSQL (For distributed)**
```python
# Future (v4.0+):
from grimperium.backends import PostgreSQLPersister

persister = PostgreSQLPersister("postgresql://...")
# Benefits: Full transactional, multi-server, scale unlimited
```

**Option C: Event Sourcing (For audit trail)**
```python
# Future (v3.5+):
from grimperium.backends import EventSourcingPersister

persister = EventSourcingPersister("grimperium.events")
# Benefits: Full history, temporal queries, no overwrites
```

**For Now (v2.2):**
- ✅ Use CSV with single-process
- ✅ Document this limitation (you're reading it)
- ❌ DO NOT attempt multiprocessing
- ❌ DO NOT run multiple instances on same data
- 🔄 Migrate to SQLite if molecules exceed 100K

---

### Limitation 2: Permissive Mode Can Normalize Errors

**Problem:**
```
Day 1:  1 invalid row (perceptible)
Day 5:  15 invalid rows (normalizing)
Day 20: 200 invalid rows (invisible in logs)

Without monitoring, error pattern goes undetected
→ Future data may have systematic issues
```

**Mitigation (v2.2):**
```python
# Always printed to console + persisted to file
BatchOrchestrator prints:
1. Validation warnings (if strict=False)
2. Batch summary (always)
3. Errors (always)

Files created:
~/.grimperium/validation_errors.log  # All validation errors
~/.grimperium/batch_summary.log      # Run summaries
```

**Monitoring Strategy:**

```bash
# Check for increasing error rates:
grep "Skipped" ~/.grimperium/batch_summary.log | tail -10

# Find error patterns:
grep "timestamp_added" ~/.grimperium/validation_errors.log | wc -l

# Alert if errors exceed threshold:
if [ $(wc -l < ~/.grimperium/validation_errors.log) -gt 100 ]; then
  echo "WARNING: Many validation errors. Investigate CSV quality."
fi
```

**Human Review Checklist:**
- [ ] After first run: Review validation_errors.log
- [ ] Weekly: Check if error count stable
- [ ] If errors > 5% of total: Stop and investigate
- [ ] Monthly: Analyze batch_summary.log trends

---

### Limitation 3: BatchOrchestrator is Bottleneck for Parallelism

**Current Concentrates:**
```
BatchOrchestrator
├─ Data loading
├─ Scheduling
├─ Execution (calls CalculationExecutor)
├─ Status tracking
├─ Persistence (atomic writes)
└─ Summaries
```

**Why Problematic for Parallelism:**
```python
# Today: Sequential, all in orchestrator
for mol in scheduled:
    _process_molecule(mol)  # Happens in orchestrator

# Tomorrow (multiprocessing):
# If multiple processes call _process_molecule:
# - Status updates race
# - Persistence conflicts
# - Summary counts wrong
```

**Future Refactor Strategy (NOT NOW):**

**Phase 1: Extract concerns (v3.0)**
```python
# Separate responsibility:

class RerunManager:
     """Handle retry logic only"""
     def mark_retry(self, mol_id):
         self.retries[mol_id] += 1
    """Handle retry logic only"""
    def mark_retry(self, mol_id):
        self.retries[mol_id] += 1

class PersistenceQueue:
    """Thread-safe write queue"""
    def enqueue(self, mol):
        self.queue.put(mol)  # Thread-safe
    
    def worker(self):
        while True:
            mol = self.queue.get()
            persister.save_molecule(mol)  # Atomic

class BatchOrchestrator:
    """Coordinator only"""
    def run(self):
        # Just coordinate, don't execute
        for mol in scheduled:
            executor.submit(process_mol, mol)
            persistence_queue.enqueue(mol)
```

**Phase 2: Async I/O (v3.5)**
```python
# Use asyncio instead of threads:
async def run_batch():
    tasks = [
        executor.run_crest_async(mol)
        for mol in scheduled
    ]
    results = await asyncio.gather(*tasks)
```

**Phase 3: Distributed (v4.0)**
```python
# Move to task queue (Celery):
@task
def process_molecule(mol_id):
    # Runs on worker process
    mol = load_molecule(mol_id)
    calculate(mol)
    save_molecule(mol)

# Enqueue tasks
for mol_id in scheduled_ids:
    process_molecule.delay(mol_id)
```

**For Now (v2.2):**
- ✅ BatchOrchestrator runs sequentially
- ✅ Simple, safe, predictable
- ❌ DO NOT attempt parallelism without refactor
- 📝 See Phase 1-3 above for future

---

## Data Flow

### Typical Batch Processing Flow

```
1. User runs: grimperium batch process

2. BatchOrchestrator.run()
   ├─ Step 1: Load CSV
   │  ├─ CSVDataLoader.load_dataframe()
   │  │  ├─ Parse CSV
   │  │  ├─ Validate columns
   │  │  └─ Validate rows (strict or permissive)
   │  │
   │  └─ BatchDataManager.load_batch()
   │     └─ For each CSV row:
   │        ├─ Molecule.from_csv_dict(row)
   │        └─ MoleculeValueConverter handles:
   │           ├─ Empty values → (None, "empty")
   │           ├─ Invalid values → (None, "invalid") [log warning]
   │           ├─ Zero values → (0.0, None) [not falsy!]
   │           └─ NA markers ("nan", "none") → (None, "invalid")
   │
   ├─ Step 2: Print validation report (if errors)
   │  └─ Console: [YELLOW] ⚠️ X rows skipped
   │     File: ~/.grimperium/validation_errors.log
   │
   ├─ Step 3: Schedule
   │  └─ BatchScheduler.schedule(molecules, max_reruns)
   │     ├─ FAILED (reruns < max) → Priority 1
   │     └─ PENDING → Priority 2
   │
   ├─ Step 4: For each scheduled molecule
   │  ├─ Set status = RUNNING
   │  ├─ CalculationExecutor.run_crest(mol)
   │  │  ├─ Run external CREST process (timeout)
   │  │  ├─ Parse output
   │  │  └─ Update mol.results.crest_*
   │  │
   │  ├─ CalculationExecutor.run_mopac_top_3(mol)
   │  │  ├─ Run MOPAC on 3 conformers
   │  │  ├─ Parse output
   │  │  └─ Update mol.results.delta_1/2/3
   │  │
   │  └─ ✅ MoleculePersister.save_molecule(mol) [ATOMIC]
   │     ├─ Read current CSV
   │     ├─ Update mol's row
   │     ├─ Write to temp file
   │     └─ Atomic rename (os.replace)
   │        → CSV never corrupted even if crash
   │
   ├─ Step 5: Print batch summary
   │  └─ Console: 📊 {total, complete, errors, elapsed}
   │     File: ~/.grimperium/batch_summary.log (append)
   │
   └─ Step 6: Return summary dict
      └─ CLI displays results

3. Next run: Steps 1-6 again
   → New molecules: status=PENDING
   → Failed molecules: status=FAILED, reruns++
   → Complete molecules: skipped (status=COMPLETE)
```

### Error Recovery Flow

```
Scenario: MOPAC times out for mol_001

Step 1: CalculationExecutor.run_mopac_top_3(mol_001)
        ↓
        raises CalculationError("MOPAC timeout")

Step 2: BatchOrchestrator._process_molecule catches error
        ├─ mol_001.meta.reruns += 1
        ├─ mol_001.meta.status = FAILED (if reruns < max_reruns)
        │  or SKIPPED (if reruns >= max_reruns)
        ├─ mol_001.results.error_message = "MOPAC timeout"

Step 3: MoleculePersister.save_molecule(mol_001) [ATOMIC]
        → CSV updated with FAILED status + retry count

Next Run:
  BatchScheduler.schedule() sees:
  ├─ mol_001 status=FAILED, reruns=1, max_reruns=3
  └─ Eligible for rerun → scheduled

Retry:
  CalculationExecutor tries again
  → If succeeds: status=COMPLETE
  → If fails again: reruns=2, FAILED again
  → After 3 failures: status=SKIPPED (no more retries)
```

### CSV Format Evolution

**v2.0: Initial columns**
```csv
mol_id,status,reruns,smiles,multiplicity,charge,nheavy,timestamp_added,H298_cbs,batch_id
```

**v2.1: Added optional columns (backward compatible)**
```csv
...,timestamp_started,timestamp_completed,crest_status,crest_conformers_generated,
crest_time,mopac_status,mopac_time,delta_1,delta_2,delta_3,most_stable_hof,error_message
```

---

## Future Considerations

### If You Need to Parallelize

**Prerequisites (MANDATORY before parallelizing):**

1. **Add File Locking (MANDATORY)**
   ```python
   from filelock import FileLock
   
   with FileLock(self.csv_path + ".lock"):
       self.persister.save_molecule(mol)
   ```

   class RerunTracker:
   ```python
   # Separate thread-safe retry tracking
   class ReruntimeTracker:
       def __init__(self):
           self.lock = threading.Lock()
           self.retries = {}
       
       def increment_rerun(self, mol_id):
           with self.lock:
               self.retries[mol_id] = self.retries.get(mol_id, 0) + 1
   ```

3. **Make Summary Updates Atomic (MANDATORY)**
   ```python
   # Use Queue for summary updates
   self.summary_queue = queue.Queue()
   
   # Worker thread
   def summary_worker():
       while True:
           event = self.summary_queue.get()
           self.summary[event['key']] += event['value']
   ```



### If You Need to Scale Beyond CSV

**Recommended migration path:**

1. **100K molecules:** Optimize CSV
   - Index by mol_id
   - Cache hot rows
   - Maybe: JSONL (append-only) for new writes

2. **1M molecules:** SQLite
   ```python
   # Local SQLite DB (ACID, indexes, queries)
   # No network, runs on same machine
   db = sqlite3.connect("grimperium.db")
   # Full transactional support
   # Can use with multiprocessing
   ```

3. **10M molecules:** PostgreSQL
   ```python
   # Remote database
   # Full enterprise features
   # Distributed writes
   ```

4. **100M+ molecules:** Distributed system
   - Kafka for event log
   - Spark for analytics
   - Time-series DB for metrics

### If You Need Real-Time Monitoring

```python
# Add metrics collection:
class Metrics:
    def __init__(self):
        self.total_processed = 0
        self.total_errors = 0
        self.avg_time_per_mol = 0

# Export to monitoring system:
# - Prometheus (/metrics endpoint)
# - CloudWatch
# - DataDog
```

### If You Need Interactive Dashboard

```python
# Add REST API:
from fastapi import FastAPI

app = FastAPI()

@app.get("/api/batch/{batch_id}/status")
def get_batch_status(batch_id: str):
    return {"status": "running", "progress": 42}

# Serve web dashboard
# Real-time updates via WebSocket
```

---

## Decision Log

### Decision: Why CSV Over SQLite (for now)?

**Considered:** SQLite, PostgreSQL, MongoDB, JSONL

**Why CSV won:**
```
✅ Human-readable (git-friendly)
✅ Simple (no server, no setup)
✅ Auditable (version control friendly)
✅ Works offline
✅ Excel-compatible (users can open/edit)
❌ No transactions (documented)
❌ No concurrency (documented)
```

**When to switch:**
- Molecules > 100K
- Need concurrent access
- Need transactional safety
- Need querying (SQLite wins)

---

### Decision: Why Atomic Writes?

**Problem:** CSV corruption if crash mid-write
```
Process writes:
  Line 1: mol_001,...   ✓ written
  Line 2: mol_002,...   ✓ written
  CRASH
  Line 3: mol_003,...   ✗ not written (corruption!)
```

**Solution:** Write to temp, atomic rename
```
Process writes to:
  ~/.temp_xyz.csv
  ├─ Line 1: mol_001
  ├─ Line 2: mol_002
  ├─ Line 3: mol_003
  └─ Sync to disk
  
os.replace(temp, actual)  ← Atomic at OS level

If crash during replace:
  → Old file still valid
  → New file half-written, orphaned
```

**Guarantee:** Read after crash = always valid state

---

### Decision: Why Settings as TOML?

**Considered:** JSON, YAML, INI, Python files

**Why TOML:**
```
✅ Human-readable
✅ Typed values
✅ Standard format
✅ Comments allowed
✅ Validated at parse time
❌ Fewer editors than JSON
```

**Example:**
```toml
[calculation]
max_reruns = 3
crest_timeout_minutes = 30

[database]
data_dir = "data"
csv_file = "thermo_pm7.csv"
```

---

### Decision: Why Typed Dataclasses?

**Compared to:** dict access, Pydantic models

**Why dataclasses:**
```
✅ IDE hints (autocomplete)
✅ mypy validation (type checking)
✅ Fast (no validation overhead)
✅ Lightweight
✅ Standard library (Python 3.7+)
❌ Less validation than Pydantic
```

**Tradeoff:** Simplicity over features

---

## Appendix: Command Reference

### Run Batch Processing

```bash
# Standard (development, strict validation)
grimperium batch process

# Production (skip invalid rows, continue)
grimperium batch process --permissive

# Dry-run (no calculations, just load + schedule)
grimperium batch process --dry-run

# With custom settings
grimperium batch process --max-reruns 5 --crest-timeout 60
```

### View Logs

```bash
# Last 10 validation errors
tail -10 ~/.grimperium/validation_errors.log

# All summaries
cat ~/.grimperium/batch_summary.log

# Watch logs in real-time
tail -f ~/.grimperium/batch_summary.log
```

### Reset Configuration

```bash
# Show current settings
grimperium settings show

# Reset to defaults
rm ~/.grimperium/config.toml
```

### Database Operations

```bash
# Check CSV integrity
grimperium validate --strict

# Check CSV integrity (permissive)
grimperium validate --permissive

# Count molecules
grimperium info --count

# Export to other format
grimperium export --format json > molecules.json
```

---

## Questions & Answers

**Q: What if CSV gets corrupted?**
A: With atomic writes (v2.2), crash-safe. If corrupted manually: restore from git history.

**Q: Can I run two instances?**
A: NO. Not supported. You'll get random overwrites.

**Q: Can I use multiprocessing?**
A: NO. Not safe. Will corrupt CSV. Refactor Phase 1 required.

**Q: Can I edit CSV manually?**
A: YES, but: 1) Stop grimperium, 2) Edit, 3) Run grimperium validate --strict to check, 4) Resume.

**Q: How do I back up?**
A: CSV is your source. Commit to git. Or: `cp data/thermo_pm7.csv data/thermo_pm7.csv.backup.$(date +%s)`

**Q: What if molecule has delta_1=0.0?**
A: Valid! v2.2 handles correctly (not converted to None).

**Q: What about zero values in results?**
A: Also valid (energy can be zero). MoleculeValueConverter preserves.

---

**Document Version:** 2.2.0
**Last Updated:** 2026-01-19
**Next Review:** 2026-06-19 (or after 100K molecules processed)
