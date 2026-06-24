# Grimperium Lessons

[2026-06-24] Context: adding Method B model metadata defaults during PR5 CLI integration
Mistake: importing `grimperium.calculation.methods.feature_schema` at module import time from `grimperium.ml.persistence` created a cycle through `feature_schema -> ml.features -> ml.__init__ -> ml.persistence`.
Rule: modules imported by `grimperium.ml.__init__` must not import calculation method registries or feature-schema catalog modules at module load time; use local imports inside functions when metadata defaults need those catalogs.
