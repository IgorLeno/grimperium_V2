# Known Issues — Grimperium

## CBS_SUSPECT: Molecules with anomalous H298_cbs in source dataset

**Status:** Known, documented, NOT corrected (upstream data)
**Discovered:** March 2026 (Batch 12)
**How found:** PM7 Baseline — MARE = 757 kcal/mol vs P50 = 5.14 kcal/mol
(mean >> median by 147x is a classic outlier signal)

### Description

13 molecules in `thermo_pm7.csv` carry `H298_cbs` values in the range of
−17,000 to −145,000 kcal/mol. A CHON molecule with 7–11 heavy atoms has a
physically plausible H298_cbs range of roughly −300 to +50 kcal/mol.

Most probable cause: CBS total electronic energy in Hartrees stored without
conversion to kcal/mol (e.g., −232 Ha × 627.5 kcal/mol/Ha ≈ −145,620 kcal/mol).
Origin is upstream of this repository and is not traceable from the current import
history.

### Quantified Impact (without filter)

| Metric | Corrupted | Clean (13 rows removed) |
|---|---|---|
| MARE | 757 kcal/mol | 6.22 kcal/mol |
| Bias | +746 kcal/mol | −5.00 kcal/mol |
| R² | −0.0068 | 0.9845 |

### Treatment

`cbs_quality_flag = "SUSPECT"` marks these 13 rows in the CSV.

| File | Symbol | Filter |
|---|---|---|
| `databases_view.py` | `_handle_pm7_baseline` | `_filter_suspect_rows()` |
| `database_analyzer.py` | `_target_delta_stats` | `cbs_quality_flag == "OK"` |
| `database_analyzer.py` | `_top_delta_outliers` | `cbs_quality_flag == "OK"` |
| `ml/data_loader.py` | `load_ml_data` | `cbs_quality_flag == "OK"` |

Rows are **retained in the CSV for traceability** but excluded from all
analysis and training.

### Reproduce

```python
import pandas as pd
df = pd.read_csv("data/thermo_pm7.csv")
suspects = df[df["cbs_quality_flag"] == "SUSPECT"][
    ["mol_id", "smiles", "nheavy", "H298_cbs"]
]
print(suspects)  # 13 rows
```

## Sync delivery / journal recovery (resolved 2026-07-12)

Previously, `/report/success|failure` applied results outside `ResultLedger`.
A lost HTTP ACK left the worker offline queue intact; a later `/sync_results`
retry could double-apply (especially `reruns`). Startup recovery also treated
any non-RUNNING status as proof of apply.

**Current contract:** one transactional path (`SyncResultApplicationService`),
per-item sync outcomes (including `conflict` / `stale_attempt` without aborting
the batch HTTP response), `attempt_id` leases on claim (cleared on AoN reset /
reclaim), immutable `result_id` reservation from first prepare, at most one
terminal result per attempt, PREPARED resume gated on matching lease and
`OperationKind`, legacy payloads rejected on terminal OK/Skip, `force_skip`
without rerun increment, exact journal proof before commit, durable idempotent
worker dead-letter for `conflict`/`stale_attempt` before queue confirm, and
heartbeat ownership/`attempt_id` validation with worker-side lease-loss
signaling. `BatchView` scientific outputs live under `runs/<run_id>/` with safe
finalization and selection compensation.

### Remaining known limitations (transactional path)

- CREST/MOPAC subprocesses are not cooperatively cancelled when a lease is lost
  mid-calculation; the worker gates publication after the pipeline returns.
- Corrupt/truncated JSONL lines in dead-letter or journal may still abort load
  for that line (skipped for dead-letter; journal still fails hard on bad JSON).
