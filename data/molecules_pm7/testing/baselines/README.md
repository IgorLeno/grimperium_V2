# CREST-PM7 Pipeline Baselines

This directory contains baseline expectations for validating CREST-PM7 pipeline results.

## File Structure

### `phase_a_baseline.json`
**Purpose**: Phase A validation baseline - 3 simple molecules for smoke testing

**Scope**: Contains methane, ethane, and benzene to validate perfect pipeline execution with minimal conformational complexity.

**Success Criteria**: Requires 100% success rate for all metrics:
- `min_success_rate`: 1.0 (100% overall pipeline success)
- `min_hof_extraction_rate`: 1.0 (100% HOF extraction success)
- `min_baseline_pass_rate`: 1.0 (100% HOF match within tolerance)
- `min_grade_ab_rate`: 1.0 (100% grade A or B)

**Molecules**: 3 molecules (methane, ethane, benzene)

### `extended_baseline.json`
**Purpose**: Extended test set for Phase B/C and edge cases

**Scope**: Contains 17 additional molecules beyond Phase A, including:
- **Phase B** (14 molecules): Standard organic molecules with moderate conformational complexity
- **Phase C** (3 molecules): Challenging edge cases (azulene, cubane, invalid structure)

**Success Criteria by Phase**:
- **Phase B**: 95% success rate, 90% baseline pass rate (allows for occasional convergence issues)
- **Phase C**: 85% success rate, 80% baseline pass rate (edge cases with lower thresholds)

**Molecules**: 20 molecules total (17 + 3 from Phase A = 20 unique structures tracked)

### `phase_a_expected.json` (legacy)
**Status**: Kept for backward compatibility, but new code should use `phase_a_baseline.json`

## Metric Definitions

- **min_success_rate**: Overall pipeline success (CREST+MOPAC complete without critical errors)
- **min_hof_extraction_rate**: HOF value successfully extracted from MOPAC output
- **min_baseline_pass_rate**: Extracted HOF matches expected baseline within tolerance (±2.5 kcal/mol)
- **min_grade_ab_rate**: Results receive quality grade A or B (not C or FAILED)

## Usage

### Phase A Testing
```bash
python scripts/phase_a_quick_test.py
```
Automatically loads `phase_a_baseline.json`

### Manual Baseline Validation
```bash
python scripts/utils/baseline_validator.py results.json data/molecules_pm7/testing/baselines/phase_a_baseline.json
```

## Tolerance

All baselines use ±2.5 kcal/mol tolerance for HOF matching, defined as `tolerance_kcal_mol: 2.5` in each JSON file.

## Source

All HOF values are from MOPAC PM7 literature values and validated experimental data.

