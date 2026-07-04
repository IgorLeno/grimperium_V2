# PR7 Package Boundary Review

**Status:** Decision recorded — extraction deferred to a dedicated task.  
**Context:** Required DoD item from `docs/plans/a-estrutura-atual-confirma-smooth-prism.md` (Specs §11, §819).  
**Author:** IgorLeno  
**Date:** 2026-07-03

---

## 1 — Purpose

Record the dependency matrix and location decision for the "Results"
capability: the code that renders analysis reports and error-analysis
charts for trained ML models. This review answers: **where should that
code live, and what can it import?**

---

## 2 — What is the "Results" capability?

The primary artefact is `src/grimperium/cli/views/results_view.py`.
It is the single file responsible for all post-training analysis output:
error analysis, HTML report generation, chart production, and model
metadata display.

### Current surface

| Location | Role |
|---|---|
| `cli/views/results_view.py` | Results view controller and all rendering logic |
| `ml/error_analysis.py` | Statistical error analysis (lazy-imported) |
| `ml/html_report.py` | HTML report generation (lazy-imported) |
| `ml/charts.py` | Chart generation (lazy-imported) |
| `ml/persistence.py` | Model metadata loader (top-level import) |

`results_view.py` is a CLI view: it inherits `BaseView`, uses Rich
panels, and lives inside `cli/views/`.  Its four `ml/` dependencies are
all imported lazily inside methods (except `persistence`, which is
imported at module level on line 21).

---

## 3 — Dependency matrix

The table below lists what a dedicated `results` module/package **may
and may not** import, derived from the clean-boundary rules in the
architectural plan.

| Category | Examples | Allowed? |
|---|---|---|
| Standard library | `json`, `pathlib`, `typing` | ✅ |
| Scientific / numeric | `numpy`, `pandas`, `matplotlib` | ✅ |
| Rich / CLI rendering | `rich.*` | ✅ (while in CLI) |
| ML layer | `ml/persistence`, `ml/error_analysis`, `ml/charts`, `ml/html_report` | ✅ read-only consumer |
| Canonical contract (PR1) | `calculation/contracts/`, `grimperium.DictStrAny` | ✅ |
| CREST / MOPAC / PM7 | `crest_pm7/*` | ❌ |
| Subprocess / xTB / workers | any subprocess, worker, server | ❌ |
| Core orchestration | `core/delta_learning`, `core/batch_orchestrator` | ❌ |
| Circular back-imports | `cli/controllers`, `cli/views/models_view` | ❌ (except through controller protocol) |

The Results capability is a **read-only display consumer**: it reads
trained model artefacts, pre-computed CSVs, and `ml/` analytics — it
does not write science, does not run pipelines, and does not own the
training loop.

---

## 4 — Location options

### Option A: `packages/grimperium-results/` (separate installable)

- Separate `pyproject.toml`, own version, own CI entry.
- Strict import boundary enforced by package isolation.
- **Cost:** Creates the `packages/` tree, which does not currently
  exist. Adds Poetry workspace configuration, import-path changes, CI
  matrix expansion, and publishing overhead. This is a non-trivial
  infrastructure investment.
- **Benefit:** Hard boundary — `crest_pm7` cannot accidentally creep in
  at import time.
- **Verdict:** Disproportionate to the current codebase size. Revisit
  when the project nears a public release with distinct install targets.

### Option B: `src/grimperium/results/` (sub-package, same wheel)

- New sub-package inside the existing wheel.
- Boundary enforced by convention + linting (e.g. `ruff` import rules or
  a custom check).
- Move `results_view.py` → `grimperium/results/view.py` (or
  `results/cli_view.py`).
- Move `ml/error_analysis.py`, `ml/html_report.py`, `ml/charts.py` →
  `grimperium/results/` once they are used exclusively there.
- **Cost:** Rename and re-export exercise; all `cli/views/` imports of
  `results_view` must update. Low risk, one focused PR.
- **Benefit:** Clearer physical location for results logic; easier to
  extract later if needed.
- **Verdict:** Preferred. Lower cost, no new infrastructure, boundary
  is explicit in the directory name.

---

## 5 — Decision

**Adopt Option B (`src/grimperium/results/`).**

Rationale:
1. The `packages/` workspace adds non-zero packaging overhead with no
   installation-time benefit at the current project scale.
2. `src/grimperium/results/` is the minimal correct signal: the directory
   name documents the boundary, and a `ruff` boundary rule can enforce it.
3. `results_view.py` already has clean lazy imports for `ml/*` — the
   structural move is low-risk.
4. If a separate installable is ever needed (e.g. a headless analytics
   server), the sub-package structure makes the extraction straightforward.

**Extraction is deferred** to a dedicated task outside this plan.  No
code moves in this PR.

---

## 6 — Import cycle risk assessment

Current cycle risk: **none observed**.

- `results_view.py` lazily imports `ml/*` — no cycle through `cli`.
- `ml/` does not import from `cli/`.
- Moving to `grimperium/results/` does not introduce new cycles as long
  as `results/` imports from `ml/` only (never the reverse).

One import to resolve at extraction time: `ml/persistence.py` is
imported at module level in `results_view.py` (line 21). That import
would move cleanly into `grimperium/results/`.

---

## 7 — Action items at extraction time

When a dedicated extraction PR is opened:

- [ ] Create `src/grimperium/results/__init__.py`
- [ ] Move `cli/views/results_view.py` → `results/cli_view.py` (or keep
  the `cli_view` name to signal it is still a CLI-facing class)
- [ ] Update `cli/views/__init__.py` re-export
- [ ] Consider co-locating `ml/error_analysis`, `ml/html_report`,
  `ml/charts` under `results/` if they have no other consumers in `ml/`
- [ ] Add a `ruff` per-file import restriction rule to block
  `crest_pm7.*` from `grimperium.results.*`
- [ ] Verify `mypy --strict` passes end-to-end
- [ ] Update `CLAUDE.md` architecture section

---

## 8 — Analysis-only boundary: ResultsView, result_evaluator, regression tests

The Results capability is **analysis/display only**. Two guarantees underpin
the boundary and must survive any future extraction:

- **`ResultsView` does not train and does not dispatch predictions.** It reads
  trained-model artefacts and pre-computed CSVs and renders them (error
  analysis, charts, HTML report). It never calls the training loop
  (`ml/trainer.py`), never runs the delta-learning fit, and never invokes a
  runner/pipeline to produce new science. Prediction dispatch belongs to the
  calculation layer (`calculation/runners/`, `cli/calc_pipeline.py`), not here.
- **`crest_pm7/result_evaluator.py` is operational-only.** After PR7a it holds
  only operational success/quality evaluation used by the pipeline; the
  scientific baseline-evaluation logic lives in the regression harness at
  `tests/regression/baseline_evaluation.py` (exercised by
  `tests/regression/test_baseline_evaluation.py`). The Results capability must
  not import `result_evaluator` to recompute science — that separation is what
  keeps analysis reproducible and free of pipeline side effects.

Regression coverage that pins this boundary:

- `tests/regression/test_baseline_evaluation.py` — scientific baseline metrics
  live in the test harness, not in the display/analysis path.
- `tests/unit/test_results_view_*` and `tests/cli/test_results_view.py` —
  exercise the view as a read-only consumer.

If the capability moves to `src/grimperium/results/`, these guarantees become
enforceable with a `ruff` per-file import rule (no `crest_pm7.*`, no
`ml/trainer`, no `calculation/runners`).

---

## 9 — References

- Architectural plan: `docs/plans/a-estrutura-atual-confirma-smooth-prism.md` §11, §819
- Current view: `src/grimperium/cli/views/results_view.py`
- Operational evaluator: `src/grimperium/crest_pm7/result_evaluator.py`
- Regression baseline: `tests/regression/baseline_evaluation.py`
- PR7 gap closure: `docs/plans/serene-tumbling-bengio.md` Batch 4
