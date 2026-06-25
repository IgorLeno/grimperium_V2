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
