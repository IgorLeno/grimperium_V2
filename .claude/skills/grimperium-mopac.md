---
name: grimperium-mopac
description: Referência especializada para o parser MOPAC atual baseado em arquivos .out
user-invocable: false
---

# MOPAC `.out` Parsing Reference

## Source Of Truth

- `src/grimperium/crest_pm7/mopac_descriptors.py`
- `src/grimperium/crest_pm7/mopac_optimizer.py`
- `tests/unit/test_mopac_descriptors.py`
- `tests/integration/test_mopac_real_output.py`

## Current State

The active parser reads human-readable MOPAC `.out` files. Historical `.aux`
guidance is obsolete for current Grimperium behavior.

Important implications:

- Do not add `.aux`-specific parsing helpers as if they were active.
- HOMO/LUMO values are now read from summary output formats instead of inferred
  from eigenvalue indexing.
- Keyword normalization in the optimizer removes the need to depend on `AUX`
  output for descriptor extraction.

## Descriptor Coverage

The current parser extracts:

- dipole moment
- ionization potential
- HOMO, LUMO, and gap
- COSMO area and volume
- gradient norm
- SCF cycles
- point group
- runtime

Heat of formation is handled elsewhere in the pipeline.

## Supported HOMO/LUMO Formats

The parser currently supports multiple `.out` summary variants, including:

- `HOMO LUMO ENERGIES (EV) = ...`
- `HOMO (SOMO) / LUMO ENERGIES (EV) = ...`
- separate `HOMO ENERGY (EV)` and `LUMO ENERGY (EV)` lines

When changing parsing behavior, update both unit and integration coverage.
