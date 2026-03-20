# Grimperium Context Snapshot

This is a non-authoritative orientation note. The active local behavior contract
is [`CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md).

## Current Repo Shape

- Scientific Python project for thermodynamic-property prediction with
  delta-learning.
- Core package layout:
  `src/grimperium/{cli,crest_pm7,core,data,ml,models,utils}`.
- Active datasets in `data/`:
  `thermo_cbs_chon.csv` and `thermo_pm7.csv`.

## Pipeline Notes

- The batch CSV backend is defined by `BatchCSVManager`.
- Current MOPAC descriptor parsing is based on `.out` files.
- CSV delta-related fields currently include `target_delta_kcalmol`,
  `abs_diff_%`, and `cbs_quality_flag`.

## Tests

- Test suites live in `tests/{unit,integration,experiments,cli,ml}`.
- On 2026-03-20, `pytest --collect-only -q` collected 691 tests in this repo.

Use this file for quick orientation only. Validate any important claim against
the actual code and tests.
