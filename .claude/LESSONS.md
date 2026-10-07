# Grimperium Lessons

[2026-06-24] Context: adding Method B model metadata defaults during PR5 CLI integration
Mistake: importing `grimperium.calculation.methods.feature_schema` at module import time from `grimperium.ml.persistence` created a cycle through `feature_schema -> ml.features -> ml.__init__ -> ml.persistence`.
Rule: modules imported by `grimperium.ml.__init__` must not import calculation method registries or feature-schema catalog modules at module load time; use local imports inside functions when metadata defaults need those catalogs.

[2026-06-24] Context: adding PR6A batch output contract tests for kcal/mol to kJ/mol conversion
Mistake: hardcoding a manually calculated expected kJ/mol string produced an incorrect test expectation.
Rule: tests for unit conversion must derive expected values from the shared conversion constant or compare numerically with `pytest.approx`, not from hand-calculated string literals.

[2026-06-25] Context: adding PR6C execution-manager tests for method metadata written to batch operational state
Mistake: asserting the whole `extra_fields` dict made the test reject valid operational row context that was written alongside the required method metadata.
Rule: tests for state-manager row updates must assert required fields and invariants, not exact whole-row payloads, unless the complete row contract is the behavior under test.

[2026-06-25] Context: reconciling PR6D legacy `thermo_pm7.csv` schema after PR6C moved operational state to `BatchStateManager`
Mistake: keeping worker assignment columns in `BatchCSVManager` let the legacy CSV grow past the 61-column scientific contract.
Rule: after migrating a responsibility to a split manager, remove the old schema columns and writer APIs from the legacy owner, and add boundary tests that assert the moved fields exist only in the new owner.

[2026-07-03] Context: adding watchdog startup-recovery regression tests on Python 3.14
Mistake: testing the watchdog by sleeping briefly and cancelling the infinite task made the test nondeterministic and hit the known `asyncio.to_thread` hang in this checkout.
Rule: tests for watchdog startup recovery must call a one-shot recovery helper and patch `asyncio.to_thread` to run synchronously; do not rely on sleep-plus-cancel timing for infinite async loops.

[2026-10-04] Context: adding stored CREST outputs as test fixtures
Mistake: `git add` of the fixture directory silently skipped `crest_conformers.xyz` because `.gitignore` ignores CREST output names everywhere; the commit passed locally but would fail on a clean checkout.
Rule: after staging fixtures copied from tool output, compare `git status`/`git show --stat` against the files on disk, and run the new test from a `git archive HEAD` export before reporting; unignore fixtures with an explicit `!tests/fixtures/...` rule, never `git add -f`.

[2026-10-05] Context: conformer-validation pilot, arms B/C on an imine (cbs_01102)
Mistake: the CONFPASS port passed 13 golden cases that were all C/H/O, then crashed on the first N-H imine because current RDKit keeps a stereo-defining hydrogen that RDKit 2020.09.5 drops; the original's NH2 hydrogen choice also turned out to follow `PYTHONHASHSEED`.
Rule: golden files for a port validated on an older library must cover every heteroatom class present in the target dataset, and the reference run must pin `PYTHONHASHSEED` whenever the original iterates over sets.

[2026-10-05] Context: retrying failed pairs of the conformer-validation pilot
Mistake: `run --limit 8 --retry-errors` planned 7 new molecules (new CREST runs) because `--limit` counts molecules with pending arms, not the first N of the sample.
Rule: before any resume or retry of `scripts/conformer_validation run`, read the `--dry-run` plan (pairs and CREST count) and narrow with `--arms`/`--limit` until it matches exactly the intended pairs.

[2026-10-07] Context: building the experimental ΔHf table from the ATcT main table
Mistake: trusting the image `alt` SMILES as the structure selected the benzene cation (233 kcal/mol) as benzene, because ATcT ions often carry a neutral SMILES; the error only showed up as a 208 kcal/mol CBS outlier.
Rule: when a source's display SMILES defines identity, cross-check charge (and formula) against the source's own metadata before screening, and inspect the largest outliers of any reference comparison before reporting statistics.

[2026-10-07] Context: NIST WebBook fetcher, search pages parsed by regex after a 10-formula pilot
Mistake: the pilot passed, but the full run hit result names with markup (`<sup>`) and single-species pages without section links; species were skipped silently as "truncated"/"unknown".
Rule: a scraper must count every page kind it classifies and treat any non-expected kind as a failure to inspect; re-run the parser over the whole cache (offline) after the fetch, before building on it.
