# Fix: LibreOffice CSV Float Display + Add HOMO/LUMO Test Coverage

## Context

HOMO/LUMO/GAP values appear as integers (-8696, -81, 8615) when opening the CSV
in LibreOffice Calc on pt_BR locale. The **code is correct** — raw CSV contains
proper decimals (-8.696, -0.081, 8.615). The issue is that pt_BR uses dot as
thousands separator, so `8.696` (dot + 3 digits) is read as "eight thousand six
hundred ninety-six".

This plan fixes the display ambiguity and adds missing test coverage.

## Changes

### 1. Add `float_format='%.6f'` to CSV writing

**File:** `src/grimperium/core/csv_utils.py:59`

```python
# BEFORE
df.to_csv(f, index=False)

# AFTER
df.to_csv(f, index=False, float_format="%.6f")
```

This ensures all floats have 6 decimal places (e.g., `-8.696000`), which
LibreOffice cannot misinterpret as thousands separators. The overhead is ~4
bytes per float — negligible for 29k rows.

### 2. Add HOMO/LUMO unit tests

**File:** `tests/unit/test_mopac_descriptors.py`

Add a new `TestHomoLumoParsing` class with 3 tests:

1. **Standard format** — parse `-8.696 -0.081` from standard MOPAC output,
   validate HOMO, LUMO, and GAP values
2. **Various decimal places** — parse values with 7 decimal places to ensure
   no truncation
3. **Near-zero values** — parse `-0.001 0.002` to test edge case handling

Uses existing `_parse_out_file` import (already present in test file line 9).

## Files Modified

| File | Change |
|------|--------|
| `src/grimperium/core/csv_utils.py` | Add `float_format="%.6f"` to `to_csv()` call (line 59) |
| `tests/unit/test_mopac_descriptors.py` | Add `TestHomoLumoParsing` class with 3 tests |

## Verification

```bash
# 1. Run new + existing tests
pytest tests/unit/test_mopac_descriptors.py -v

# 2. Verify CSV output format (re-run pipeline or manual check)
python3 -c "
import pandas as pd
df = pd.read_csv('data/thermo_pm7.csv')
df.to_csv('/tmp/test_format.csv', index=False, float_format='%.6f')
" && head -3 /tmp/test_format.csv | cut -d',' -f46,47,48
# Expected: -8.696000,-0.081000,8.615000

# 3. Full test suite
pytest tests/ -v --cov=src/grimperium

# 4. Linting
ruff check src/grimperium/core/csv_utils.py
```
