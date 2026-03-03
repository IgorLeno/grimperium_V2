# Plan: Update PM7Config Default Values for Higher Precision

## Context

The CREST-PM7 pipeline currently uses conservative default settings designed to be fast rather than accurate. For production-quality scientific results, three defaults need to change:

- **`xtb_preopt`**: Pre-optimizing geometry with xTB before CREST reduces bad conformations and improves convergence.
- **`crest_opt_level`**: "tight" (`2`) is more rigorous than "normal" (`1`), especially for flexible molecules.
- **`mopac_precise_scf`**: The MOPAC `PRECISE` keyword tightens SCF convergence, improving reliability of HOF, HOMO/LUMO, and dipole values.

Expected trade-off: ~20–40% longer calculation times per molecule, justified by higher data quality.

---

## Files to Modify

### 1. `src/grimperium/crest_pm7/config.py`

Change three field defaults in `PM7Config`:

| Line | Field | Old Default | New Default |
|------|-------|-------------|-------------|
| 146 | `crest_opt_level` | `1` | `2` |
| 147 | `crest_optlev_label` | `"normal"` | `"tight"` |
| 158 | `mopac_precise_scf` | `False` | `True` |
| 162 | `xtb_preopt` | `False` | `True` |

**Note on `crest_optlev_label`**: `__post_init__` (lines ~245–253) already derives and overwrites this field from `crest_opt_level` via the map `{0: "loose", 1: "normal", 2: "tight"}`. Updating both the field default and `crest_opt_level` keeps them consistent at declaration time and avoids any transient inconsistency warnings.

### 2. `tests/unit/test_settings_manager.py`

**Line 143** asserts on the old default and will fail:
```python
# Before
assert settings_dict["crest_xtb_preopt"] is False
# After
assert settings_dict["crest_xtb_preopt"] is True
```

No other test assertions on these specific defaults were found. The assertion at line 206 (`manager.mopac.precise is False`) tests `MOPACSettings.precise`, a separate class unrelated to `PM7Config.mopac_precise_scf`.

---

## Implementation Steps

1. Edit `src/grimperium/crest_pm7/config.py`:
   - `crest_opt_level: int = 2`
   - `crest_optlev_label: str = "tight"`
   - `mopac_precise_scf: bool = True`
   - `xtb_preopt: bool = True`

2. Edit `tests/unit/test_settings_manager.py` line 143:
   - `assert settings_dict["crest_xtb_preopt"] is True`

---

## Verification

```bash
# Type checking
mypy --strict src/grimperium/crest_pm7/config.py

# Linting
ruff check src/grimperium/crest_pm7/config.py

# Full test suite
pytest tests/ -v --tb=short
```

Expected outcome: all tests pass with the updated assertion. The `__post_init__` validation will confirm `crest_opt_level=2` maps correctly to `crest_optlev_label="tight"`.
