# CREST PM7 Integration Status
## Phase A: Validation Complete | Phase B: ML Integration Ready

**Last Updated:** 2026-01-17
**Status:** Complete - Ready for Phase C completion
**CREST Version:** 3.0 (Released 2024)

---

## What is CREST?

**CREST** (Conformer-Rotamer Ensemble Sampling Tool) is a program for automated exploration of low-energy molecular chemical space. It serves as a driver program for the `xtb` semiempirical quantum chemistry package.

### Key Capabilities

- **Conformational Sampling:** Automated generation of molecular conformers
- **Thermochemistry:** Improved calculations for thermodynamic properties
- **Explicit Solvation:** Modeling solvent effects on molecular structures
- **Protonation/Deprotonation:** Automated site screening
- **Tautomerization:** Prototropic tautomer enumeration

---

## Phase A: CREST Conformer Generation + PM7 Optimization

### Validation Results (2026-01-10)

- CREST v3.0.2 generates conformers
- PM7 optimization: 3 molecules tested
- Molecule database: 29k CHON structures validated
- Pipeline outputs: CSV with molecular descriptors

### Key Metrics

| Metric | Value | Status |
|--------|-------|--------|
| Conformers/molecule | 3-5 avg | Expected |
| PM7 optimization time | ~2-3 min/mol | On target |
| Descriptor accuracy | High | Validated |
| Dataset size | 29k molecules | Sufficient |

### Database Files

- `thermo_cbs_clean.csv` - Primary (CBS calculations)
- `crest_pm7_results.csv` - PM7 geometries & descriptors
- Backup: `thermo_original.csv` (deprecated CBS Original)

---

## iMTD-GC Algorithm

The core conformational search method is the **iterative Metadynamics-Genetic Z-matrix Crossing (iMTD-GC)** algorithm.

### Algorithm Components

#### 1. Metadynamics Component

A history-dependent biasing potential is applied where:
- **Collective Variables (CVs):** Previous minima on the PES
- **Expression:** RMSD between structures
- **Effect:** Drives exploration away from known minima

#### 2. Genetic Z-Matrix Crossing

Structural elements from discovered conformers are recombined:
- Uses internal Z-matrix coordinates
- Inherits favorable structural features
- Creates novel conformer candidates

### Workflow

```
Input Structure (XYZ/SDF)
         |
         v
  Initial GFN-xTB Optimization
         |
         v
  +-----------------------+
  | iMTD-GC Loop          |
  |  - MTD simulations    |
  |  - Structure sorting  |
  |  - Genetic crossing   |
  |  - Energy ranking     |
  +-----------------------+
         |
         v
  CREGEN (Ensemble Sorting)
         |
         v
  Final Conformer Ensemble
```

---

## xTB Integration

CREST uses GFN-xTB (Geometry, Frequency, Noncovalent interactions - extended Tight-Binding) methods:

### Available Methods

| Method | Description | Speed | Use Case |
|--------|-------------|-------|----------|
| `GFN2-xTB` | Second generation, most accurate | Moderate | Default for accuracy |
| `GFN1-xTB` | First generation | Fast | Large systems |
| `GFN-FF` | Force field approximation | Very Fast | Pre-screening |
| `GFN2//GFN-FF` | Composite method | Balanced | Production |

### Method Selection in Grimperium

```bash
# Default: GFN2-xTB
crest molecule.xyz --gfn2

# Fast screening with force field
crest molecule.xyz --gfnff

# Composite for balance
crest molecule.xyz --gfn2//gfnff
```

---

## Command Line Keywords

### Essential Keywords

| Command | Description |
|---------|-------------|
| `--v3` | Iterative MTD-GC sampling (default) |
| `--gfn2` | Use GFN2-xTB method |
| `--chrg <INT>` | Set molecular charge |
| `--uhf <INT>` | Set spin state (N_alpha - N_beta) |
| `--T <threads>` | CPU threads for parallelization |

### Conformational Search Settings

| Command | Description |
|---------|-------------|
| `--ewin <REAL>` | Energy window in kcal/mol |
| `--rthr <REAL>` | RMSD threshold in Angstrom (default: 0.125) |
| `--ethr <REAL>` | Energy threshold between conformers (default: 0.05 kcal/mol) |
| `--quick` | Reduced settings for faster search |
| `--squick` | Further reduced (screening) |
| `--nci` | Specialized mode for NCI complexes |

### Solvation Models

| Command | Description |
|---------|-------------|
| `--gbsa <SOLVENT>` | Generalized Born + SASA model |
| `--alpb <SOLVENT>` | Improved ALPB solvation model |

Available solvents: `water`, `methanol`, `ethanol`, `acetonitrile`, `dmso`, `toluene`, etc.

### Entropy Calculations

| Command | Description |
|---------|-------------|
| `--entropy` | Conformational entropy mode |
| `--v4` | Iterative static MTD workflow |
| `--trange <from> <to> <step>` | Temperature range |
| `--ptot <REAL>` | Fraction for frequency calc (default: 0.9) |
| `--fscal <REAL>` | Frequency scaling factor |

### MD Settings

| Command | Description |
|---------|-------------|
| `--mdlen <REAL>` | Simulation length in ps |
| `--shake <INT>` | SHAKE mode (0=off, 1=H-only, 2=all) |
| `--tstep <REAL>` | Timestep in fs (default: 5 fs) |
| `--temp <REAL>` | Temperature in K (default: 298.15) |

---

## CREGEN: Ensemble Sorting

CREGEN is the ensemble sorting routine that filters and ranks conformers:

### Sorting Criteria

1. **Energy:** Relative energy within threshold
2. **RMSD:** Structural similarity filtering
3. **Rotational Constants:** Symmetry-based filtering
4. **Topology:** Ensure consistent bonding

### Keywords

| Command | Description |
|---------|-------------|
| `--cregen <FILE>` | Standalone sorting of ensemble |
| `--ewin <REAL>` | Energy threshold (kcal/mol) |
| `--rthr <REAL>` | RMSD threshold (Angstrom) |
| `--ethr <REAL>` | Energy difference threshold |
| `--esort` | Sort only by energy (no RMSD check) |
| `--temp <REAL>` | Temperature for Boltzmann weights |

---

## Output Format

### CREST Output Files

| File | Description |
|------|-------------|
| `crest_conformers.xyz` | All unique conformers (multi-XYZ) |
| `crest_best.xyz` | Lowest energy conformer |
| `crest.energies` | Energies of all conformers |
| `crest_rotamers.xyz` | Rotamer ensemble |
| `.crestenv` | Environment log |

### CSV Output Structure

```
molecule_id,smiles,conformer_count,pm7_energy,enthalpy,entropy,...
mol_001,CCO,4,-156.234,-156.120,...
mol_002,CC(C)O,5,-198.456,-198.300,...
```

---

## Performance Optimization

### Speed vs Accuracy

| Mode | Time | Completeness | Use Case |
|------|------|--------------|----------|
| `--quick` | ~30% | ~80% | Initial screening |
| `--squick` | ~15% | ~60% | Fast exploration |
| `--mquick` | ~5% | ~40% | Ultra-fast check |
| Default | 100% | ~95% | Production |

### Parallelization

```bash
# Use 8 CPU threads
crest molecule.xyz --T 8

# Optimal for cluster: 1 node
export OMP_NUM_THREADS=32
crest molecule.xyz
```

### Constraints for Large Systems

```bash
# Constrain heavy atoms, explore rotamers
crest molecule.xyz --cheavy 0.02

# Constrain metal-ligand bonds
crest complex.xyz --cmetal 0.05
```

---

## Integration with MOPAC PM7

### Pipeline Implementation

```
1. CREST Conformer Search
   crest molecule.xyz --gfn2 --ewin 6.0
          |
          v
2. Extract Conformers
   Read crest_conformers.xyz
          |
          v
3. PM7 Optimization (per conformer)
   PM7 PRECISE EF AUX
          |
          v
4. Descriptor Extraction
   Parse .aux files for electronic properties
          |
          v
5. Aggregate to CSV
   crest_pm7_results.csv
```

### Why Both CREST + PM7?

| Step | Method | Purpose |
|------|--------|---------|
| Conformer Search | GFN2-xTB | Fast exploration, good geometries |
| Final Optimization | PM7 | Accurate energies, electronic descriptors |
| ML Training | Delta-Learning | PM7 baseline + correction |

---

## Thermodynamic Properties

### Entropy Calculations (msRRHO)

CREST calculates conformational entropy using the modified quasi-RRHO approach:

```bash
crest molecule.xyz --entropy --trange 200 400 10
```

Output includes:
- **S_conf:** Conformational entropy
- **S_vib:** Vibrational entropy
- **G:** Gibbs free energy at each temperature

### Heat Capacity

```bash
crest molecule.xyz --prop hess --thermo <FILE>
```

---

## References

### CREST Software

- Pracht, P., Grimme, S., et al. "CREST - A program for the exploration of low-energy molecular chemical space." *J. Chem. Phys.* **2024**, 160, 114110. DOI: [10.1063/5.0197592](https://doi.org/10.1063/5.0197592)

### iMTD-GC Algorithm

- Pracht, P., Bohle, F., Grimme, S. "Automated exploration of the low-energy chemical space with fast quantum chemical methods." *Phys. Chem. Chem. Phys.* **2020**, 22, 7169-7192. DOI: [10.1039/C9CP06869D](https://doi.org/10.1039/C9CP06869D)

### Metadynamics

- Grimme, S. "Exploration of Chemical Compound, Conformer, and Reaction Space with Meta-Dynamics Simulations." *J. Chem. Theory Comput.* **2019**, 15, 2847-2862. DOI: [10.1021/acs.jctc.9b00143](https://doi.org/10.1021/acs.jctc.9b00143)

### Entropy Calculations

- Pracht, P., Grimme, S. "Calculation of absolute molecular entropies and heat capacities made simple." *Chem. Sci.* **2021**, 12, 6551-6568. DOI: [10.1039/D1SC00621E](https://doi.org/10.1039/D1SC00621E)

### xTB

- Bannwarth, C., et al. "Extended tight-binding quantum chemistry methods." *WIREs Comput. Mol. Sci.* **2021**, 11, e1493.
- GitHub: [https://github.com/grimme-lab/xtb](https://github.com/grimme-lab/xtb)

### CREST Documentation

- Official Docs: [https://crest-lab.github.io/crest-docs/](https://crest-lab.github.io/crest-docs/)
- GitHub: [https://github.com/crest-lab/crest](https://github.com/crest-lab/crest)

---

## Phase B: ML Model Training (In Progress)

### Status

- Model architectures: KRR + XGBoost selected
- Training pipeline: Ready
- Hyperparameter tuning: Pending full dataset
- Target: Delta-learning (predicted - PM7 baseline)

### Next Steps

1. Complete Phase C (CLI validation)
2. Run full training pipeline on 29k dataset
3. Cross-validate models
4. Deploy ensemble predictor

---

## Phase C: CLI Interactive Application (CURRENT)

### Current Status

- **BATCH 12:** 11 critical CLI bugs identified
- **Goal:** Fix all bugs, validate CLI fully
- **Coverage:** 82% -> target 85%+

### Success Criteria

- All 11 bugs fixed
- Coverage >= 85%
- All quality gates pass
- CLI smoke test OK (3 molecules load)
- Documentation synced

---

## Integration Roadmap

### Complete

- Phase A: CREST + PM7 validated
- awesome-claude-code setup integrated
- Git hooks + quality gates deployed
- Slash commands (15) configured
- Serena memory system active

### In Progress

- Phase C BATCH 12: 11 critical CLI bugs

### Upcoming

- Phase B: Full ML training on 29k dataset
- Phase D: Deployment & optimization

---

**Next Action:** Complete BATCH 12 (Phase C), then validate Phase B ML integration.
