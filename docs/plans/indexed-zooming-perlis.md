# Fix Dipole Parsing + Force PRECISE Keyword

## Context

Two critical bugs in the CREST-PM7 pipeline:
1. **Dipole not parsed** — `mopac_dipole_debye` returns `NaN` in 100% of cases. MOPAC2016 outputs a table format (`SUM X Y Z TOTAL`), but the regex only matches a summary format (`DIPOLE = X.XXX DEBYE`) that this version doesn't produce.
2. **PRECISE keyword missing** — Despite `PM7Config.mopac_precise_scf=True` being the default, the `.mop` files lack `PRECISE`. The config propagation path allows `False` to override. For production quality, PRECISE must always be forced ON.

## Changes

### Fix 1: Dipole parsing fallback — `src/grimperium/crest_pm7/mopac_descriptors.py`

**Lines 73-83** — Replace current dipole regex block with dual-format parser:
1. Try summary format first: `DIPOLE = X.XXX DEBYE` (backward compat)
2. Fallback to table format: extract TOTAL column from `SUM` line
3. Store result if either succeeds

### Fix 2: Force PRECISE keyword — `src/grimperium/crest_pm7/mopac_optimizer.py`

**Lines 83-88** — Replace conditional PRECISE logic:
1. Always set `precise_scf = True` (unconditional)
2. Log WARNING if `config.mopac_precise_scf=False` (user awareness)
3. Always append `PRECISE` + `SCFCRT=` to keywords

### Tests (new files)

| File | Tests |
|------|-------|
| `tests/unit/test_mopac_descriptors.py` | Dipole table format, summary format, fallback priority |
| `tests/unit/test_mopac_optimizer.py` | PRECISE always present, warning logged when config=False |
| `tests/integration/test_mopac_real_output.py` | Full descriptor extraction from real `.out` content |

## Verification

```bash
pytest tests/unit/test_mopac_descriptors.py tests/unit/test_mopac_optimizer.py tests/integration/test_mopac_real_output.py -v
ruff check src/grimperium/crest_pm7/mopac_descriptors.py src/grimperium/crest_pm7/mopac_optimizer.py
```
