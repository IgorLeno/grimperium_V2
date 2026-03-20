---
name: grimperium-csv
description: Referência do schema CSV atual e das regras de mudança do backend de batch
user-invocable: false
---

# CSV Schema Reference

## Source Of Truth

- `src/grimperium/crest_pm7/batch/csv_manager.py`
- `tests/test_csv_schema.py`

If prose and code disagree, trust those files.

## Current Shape

`BatchCSVManager.get_schema()` currently defines a 58-column schema grouped as:

- identity: `mol_id`, `status`, `smiles`
- molecular properties: `multiplicity`, `charge`, `nheavy`, `H298_cbs`,
  `H298_pm7`, `target_delta_kcalmol`, `abs_diff_%`, `cbs_quality_flag`
- batch info: `batch_id`, `timestamp`, `total_time`, `reruns`
- RDKit descriptors
- CREST settings and runtime fields
- MOPAC status, selected conformer rank, descriptors, runtime
- batch configuration fields

## Current Guardrails

- `target_delta_kcalmol` is the signed learning target.
- `abs_diff_%` is the current percent-error field.
- The selected PM7 conformer rank is tracked as `k_selected_pm7`.
- Result-reset behavior is defined by `RESULT_COLUMNS` in `csv_manager.py`.

## Removed Historical Fields

Do not treat these as current schema:

- `delta_1`, `delta_2`, `delta_3`
- `conformer_selected`
- `most_stable_hof`
- `reference_hof`
- legacy timeout and retry columns not present in the current schema tests

## When Editing The Schema

1. Update `BatchCSVManager.get_schema()` and related grouped constants.
2. Update `RESULT_COLUMNS` if reset behavior changes.
3. Update the schema tests.
4. Update any docs or changelog entry that references the CSV contract.
