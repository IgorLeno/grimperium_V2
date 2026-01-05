# Grimperium Dataset Analysis

## Overview

**Source:** `thermo_cbs_opt.csv`
**Total Molecules:** 52,837
**CBS Method:** CBS-QB3 (Complete Basis Set composite method)

## Schema

| Column | Type | Description | Range |
|--------|------|-------------|-------|
| smiles | string | SMILES molecular structure | - |
| multiplicity | int | Spin multiplicity | 1-2 |
| charge | int | Total molecular charge | -1, 0, +1, +2 |
| nheavy | int | Number of heavy atoms | 1-22 |
| H298_cbs | float | CBS enthalpy at 298K (kcal/mol) | -325,407 to +164,949 |
| H298_b3 | float | B3LYP enthalpy at 298K (kcal/mol) | -263,359 to +174,323 |

## Statistics

### H298_cbs (Target Variable)

- **Mean:** -320.37 kcal/mol
- **Std:** 7,230.27 kcal/mol
- **Min:** -325,407.00 kcal/mol
- **Max:** 164,949.00 kcal/mol
- **Median:** -25.83 kcal/mol

### Molecular Properties

**Heavy Atom Distribution:**
- Most molecules: 9-11 heavy atoms (81%)
- Range: 1-22 heavy atoms

**Charge Distribution:**
- Neutral (0): 99.95%
- Charged (-1): <0.01% (2 molecules)
- Charged (+1): 0.04% (21 molecules)
- Charged (+2): <0.01% (2 molecules)

**Multiplicity Distribution:**
- Singlet (1): 98.65%
- Doublet (2): 1.35%

## Data Quality

- **Completeness:** 100% (no missing values)
- **Duplicates:** 446 duplicate SMILES found (may be conformers or computational duplicates)
- **Invalid values:** 0
- **SMILES validity:** All valid (basic check)

## Usage Notes

1. **Train/Test Split:** Use 80/20 stratified by nheavy
2. **Delta Learning:** Use H298_b3 as semiempirical proxy (PM7 unavailable)
3. **Feature Engineering:** Extract Morgan fingerprints from SMILES
4. **Validation:** Reserve 10% for final validation (70/10/20 split)
