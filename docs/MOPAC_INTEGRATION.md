# MOPAC PM7 Integration in Grimperium

## Overview

**Program:** OpenMOPAC 2025
**Method:** PM7 (Parametrized Model 7)
**Role:** Single-point energy and property calculation
**Speed:** 0.5-5 seconds per molecule (typical)
**Accuracy:** Semiempirical (error corrected via delta-learning)

MOPAC is a semiempirical quantum chemistry package that has been continuously
developed since the 1970s. PM7 is the latest major parametrization, providing
good thermochemical accuracy at a fraction of DFT computational cost.

---

## What is PM7?

### Definition

**PM7** = Parametrized Model 7

PM7 is a semiempirical quantum chemistry method that:
- Uses minimal basis sets
- Employs parameterized integrals
- Includes D2 dispersion correction
- Is trained on experimental and high-level computational data

### Method Hierarchy

```
         Accuracy                    Speed
            ↑                          ↑
            │  CBS/CCSD(T)            │  Force Fields (MM)
            │  └─ "Gold standard"     │  └─ Seconds (large systems)
            │                         │
            │  DFT/B3LYP              │  Semiempirical (PM7)
            │  └─ ~2-5 kcal/mol       │  └─ 1-5 seconds
            │                         │
            │  Semiempirical (PM7)    │  DFT/B3LYP
            │  └─ ~5-10 kcal/mol      │  └─ Minutes to hours
            │                         │
            │  Force Fields (MM)      │  CBS/CCSD(T)
            │  └─ ~10-20 kcal/mol     │  └─ Hours to days
            ↓                         ↓
```

### PM7 Parametrization

PM7 is parametrized using:
- Experimental heats of formation
- Experimental geometries
- DFT reference data
- High-level calculations

**Supported Elements:**
H, Li, Be, B, C, N, O, F, Na, Mg, Al, Si, P, S, Cl, K, Ca, Sc, Ti, V, Cr, Mn,
Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd,
Ag, Cd, In, Sn, Sb, Te, I, Cs, Ba, La, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg,
Tl, Pb, Bi, Po, At

### Advantages

1. **Very Fast:** Seconds per molecule vs. minutes/hours for DFT
2. **Good Thermochemistry:** ~5-10 kcal/mol MAE for heats of formation
3. **Dispersion Included:** D2 correction built into parametrization
4. **Robust SCF:** Generally converges well
5. **Extensive Element Coverage:** Most of the periodic table

### Limitations

1. **Fixed Parametrization:** Cannot be improved for specific chemistry
2. **Semiempirical Approximations:** Less accurate than DFT
3. **Basis Set Limitations:** Minimal basis cannot capture all effects
4. **Rare Elements:** Some lanthanides/actinides poorly parameterized

---

## Input Format

MOPAC uses a specific text input format with the `.mop` extension.

### Input File Structure

```
Line 1: KEYWORDS (space-separated)
Line 2: Title/description
Line 3: Comment (can be blank)
Line 4+: Geometry (Z-matrix or Cartesian)
```

### Example Input File

```mop
PM7 PRECISE SCFCRT=1.0D-5
Ethanol molecule - Heat of Formation calculation
Single point from CREST best conformer

 C     0.0000000000    0.0000000000    0.0000000000
 C     1.5200000000    0.0000000000    0.0000000000
 O     2.0500000000    1.2800000000    0.0000000000
 H    -0.3600000000    1.0100000000    0.0000000000
 H    -0.3600000000   -0.5100000000    0.8700000000
 H    -0.3600000000   -0.5100000000   -0.8700000000
 H     1.8800000000   -0.5100000000    0.8700000000
 H     1.8800000000   -0.5100000000   -0.8700000000
 H     2.9400000000    1.2800000000    0.0000000000
```

### Geometry Formats

**Cartesian Coordinates (recommended):**
```
ELEMENT   X         Y         Z
C         0.000     0.000     0.000
C         1.520     0.000     0.000
```

**Z-Matrix (internal coordinates):**
```
C
C   1    1.52
O   2    1.41    1   109.5
```

**Note:** Grimperium uses Cartesian coordinates from CREST output for
simplicity and direct coordinate transfer.

---

## Critical Keywords

### Essential Keywords for Grimperium

| Keyword | Default | Effect | Grimperium |
|---------|---------|--------|-----------|
| `PM7` | N/A | Use PM7 method | **Always** |
| `1SCF` | OFF | Single-point (no optimization) | **Always** |
| `PRECISE` | OFF | 100x tighter SCF criteria | **Optional** |
| `SCFCRT=<E>` | 1.0D-4 | SCF convergence threshold | **Default** |
| `ITRY=<N>` | 1000 | Max SCF iterations | **Default** |
| `PULAY` | OFF | SCF acceleration | **If needed** |

### Single-Point Keywords (Phase A)

For Grimperium Phase A, we perform single-point energy calculations:

```
PM7 1SCF
```

This gives:
- Heat of formation (ΔHf)
- HOMO/LUMO energies
- Dipole moment
- Mulliken charges

**No geometry optimization** - we use the CREST-optimized structure.

### SCF Convergence Keywords

**Standard (default):**
```
PM7 1SCF
```

**Tighter convergence:**
```
PM7 1SCF PRECISE
# or
PM7 1SCF SCFCRT=1.0D-6
```

**For difficult cases:**
```
PM7 1SCF PULAY ITRY=2000
```

### Output Control Keywords

```
# Minimal output (faster)
PM7 1SCF AUX(0)

# Include thermodynamic data
PM7 1SCF THERMO

# Include bond orders
PM7 1SCF BONDS

# Include electrostatic potential
PM7 1SCF ESP
```

---

## MOPAC Output Structure

### Output File (.out)

```
output.out
├── Header/Version information
├── Input echo
├── SCF convergence details
│   ├── Energy at each iteration
│   └── Convergence status
├── Final Results
│   ├── HEAT OF FORMATION
│   ├── TOTAL ENERGY
│   ├── ELECTRONIC ENERGY
│   ├── CORE-CORE REPULSION
│   ├── GRADIENT (if optimizing)
│   ├── DIPOLE
│   └── NO. OF FILLED LEVELS
├── Orbital Energies
│   ├── HOMO
│   ├── LUMO
│   └── Gap
├── Mulliken Charges
└── Timing information
```

### Key Output Sections

**Heat of Formation:**
```
          FINAL HEAT OF FORMATION =       -52.26784 KCAL/MOL
```

**HOMO-LUMO:**
```
          HOMO LUMO ENERGIES (EV) =      -10.492  1.234
```

**Dipole Moment:**
```
          DIPOLE           X         Y         Z       TOTAL
          POINT-CHG.     1.234    -0.567     0.123     1.367
          HYBRID         0.456     0.234    -0.089     0.521
          SUM            1.690    -0.333     0.034     1.723
```

### Auxiliary Files

When using `AUX` keyword:
- `.aux` - Machine-readable auxiliary output
- `.arc` - Archive file with trajectory
- `.den` - Density matrix

---

## Error Handling

### Exit Codes

**Exit Code 0:** Success
- SCF converged
- All requested calculations completed

**Non-zero Exit:** Error
- Check `.out` file for details
- Common issues documented below

### Common Errors

#### "SCF DID NOT CONVERGE"

**Symptoms:**
```
 SCF CALCULATIONS DID NOT CONVERGE
 TOTAL ENERGY, E, IS UNRELIABLE
```

**Causes:**
1. Bad input geometry
2. Unusual electronic structure
3. Near-degenerate orbitals
4. High-spin states

**Solutions:**
```mop
# Add PULAY acceleration
PM7 1SCF PULAY

# Increase iterations
PM7 1SCF ITRY=2000

# Relax convergence (not recommended)
PM7 1SCF SCFCRT=1.0D-3
```

#### "GEOMETRY IS INCORRECT"

**Symptoms:**
```
 IMPOSSIBLE NUCLEAR DISTANCE
```

**Causes:**
- Atoms too close (< 0.5 Å)
- Missing atoms
- Overlapping coordinates

**Solutions:**
- Check input geometry
- Validate coordinate system
- Use better initial structure

#### "CALCULATION ABNORMALLY TERMINATED"

**Symptoms:**
```
 ** CALCULATION ABANDONED **
```

**Causes:**
- Memory exhausted
- Numerical overflow
- File I/O error
- License issue

**Solutions:**
- Check system memory
- Verify input validity
- Check disk space

---

## Grimperium Integration

### PM7 in Phase A Workflow

```
CREST Best Conformer (.xyz)
          ↓
    Convert to MOPAC format
          ↓
    MOPAC PM7 1SCF
          ↓
    Parse Output
          ├── Heat of Formation
          ├── HOMO/LUMO
          ├── Dipole
          └── SCF iterations
          ↓
    PM7Result dataclass
```

### PM7Result Dataclass

```python
@dataclass
class PM7Result:
    """Result from PM7 single-point calculation."""

    mol_id: str
    smiles: str
    timestamp: datetime

    # Energy values
    homo: float | None           # eV
    lumo: float | None           # eV
    homo_lumo: float | None      # Gap in eV
    hof: float | None            # Heat of formation (kcal/mol)

    # Properties
    dipole: float | None         # Dipole moment (Debye)
    total_energy: float | None   # Total energy (eV)

    # Convergence info
    scf_iterations: int | None   # Iterations to converge
    converged: bool = True       # SCF converged?

    # Status
    success: bool = True
    error_message: str | None = None
    quality_grade: QualityGrade = QualityGrade.A
```

### MOPAC Optimizer Module

Location: `src/grimperium/crest_pm7/mopac_optimizer.py`

```python
def optimize_conformer(
    mol_id: str,
    xyz_file: Path,
    config: PM7Config,
    timeout: float,
    nheavy: int,
    conf_index: int = 0,
) -> MOPACResult:
    """
    Run PM7 single-point on conformer.

    Args:
        mol_id: Molecule identifier
        xyz_file: Input XYZ file from CREST
        config: PM7 configuration
        timeout: Timeout in seconds
        nheavy: Number of heavy atoms
        conf_index: Conformer index

    Returns:
        MOPACResult with energies and properties
    """
    ...
```

### Energy Extraction

Location: `src/grimperium/crest_pm7/energy_extractor.py`

Multiple extraction patterns for robustness:

```python
# Pattern 1: Standard output (most reliable)
FINAL HEAT OF FORMATION =      -52.26784 KCAL/MOL

# Pattern 2: Thermodynamic summary
HEAT OF FORMATION =      -52.26784 KCAL

# Pattern 3: Compact format (older MOPAC)
HOF =      -52.26784

# Pattern 4: Error case
TOTAL ENERGY, E, IS UNRELIABLE
```

---

## Performance Expectations

### Timing by Molecule Size

**Small molecules (1-10 heavy atoms):**
- PM7 time: 0.1 - 0.5 seconds
- SCF iterations: 3-10
- Example: Methane, ethanol, benzene

**Medium molecules (10-30 heavy atoms):**
- PM7 time: 0.5 - 2 seconds
- SCF iterations: 5-15
- Example: Aspirin, caffeine, nicotine

**Large molecules (30-100 heavy atoms):**
- PM7 time: 2 - 10 seconds
- SCF iterations: 10-30
- Example: Peptides, small proteins

**Very large (100+ heavy atoms):**
- PM7 time: 10 - 60 seconds
- SCF iterations: 20-100
- May require PULAY acceleration

### Convergence Statistics

For CHON-only dataset (29,568 molecules):
- **Average SCF iterations:** 5-10
- **Average time:** 0.5 - 2 seconds
- **SCF failures:** <0.1%
- **Overall success:** >99%

---

## Delta-Learning Context

### Why PM7 + Delta-Learning?

PM7 alone has ~5-10 kcal/mol error for heats of formation.
Delta-learning reduces this significantly.

**Strategy:**
```
PM7_value = PM7(molecule)           # Fast, ~5-10 kcal/mol error
CBS_value = CBS_reference(molecule) # Slow, gold standard

delta = CBS_value - PM7_value       # What PM7 missed

# Train ML model to predict delta from molecular features
model.fit(features, delta_values)

# For new molecules:
corrected = PM7(new_mol) + model.predict(features(new_mol))
```

**Result:** CBS-quality predictions at PM7 speed.

### Phase A Role

Phase A generates PM7 values for ~30k molecules.

These become training data for Phase B:
1. PM7 values (Y_pm7)
2. CBS reference values (Y_cbs) - from Chemperium dataset
3. Molecular features (X) - from RDKit descriptors

Phase B trains: delta = Y_cbs - Y_pm7

---

## Troubleshooting

### Problem: SCF Not Converging

**Symptoms:**
```
 SCF CALCULATIONS DID NOT CONVERGE
```

**Step 1: Check geometry**
```python
# Verify structure is reasonable
from rdkit import Chem
mol = Chem.MolFromMolFile('molecule.sdf')
if mol is None:
    print("Invalid structure")
```

**Step 2: Add PULAY**
```mop
PM7 1SCF PULAY
```

**Step 3: Relax criteria**
```mop
PM7 1SCF SCFCRT=1.0D-3  # Not recommended for production
```

### Problem: Negative HOMO-LUMO Gap

**This can happen legitimately:**
- Near-degenerate orbitals
- High-spin states
- Radical species

**Action:**
- Log it as unusual
- Continue processing
- Delta-learning will handle

### Problem: Very Slow Convergence

**Normal PM7:** 3-10 iterations
**Slow:** 20-50 iterations
**Very slow:** 100+ iterations

**Solutions:**
1. Check for unusual elements
2. Verify geometry quality
3. Add PULAY acceleration
4. Consider skipping (error log)

### Problem: Energy NaN or Inf

**Causes:**
- Bad geometry (overlapping atoms)
- Numerical overflow
- Memory corruption

**Solutions:**
- Validate input coordinates
- Check for atom collisions
- Restart calculation

---

## Best Practices

### 1. Input Preparation

```python
def prepare_mopac_input(xyz_file: Path, mol_id: str) -> str:
    """Generate MOPAC input from XYZ file."""
    coords = read_xyz(xyz_file)

    # Validate coordinates
    for i, (elem1, x1, y1, z1) in enumerate(coords):
        for j, (elem2, x2, y2, z2) in enumerate(coords[i+1:], i+1):
            dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            if dist < 0.5:
                raise ValueError(f"Atoms {i} and {j} too close: {dist} Å")

    # Generate input
    lines = [
        "PM7 1SCF",
        f"Molecule {mol_id}",
        "",  # Blank comment line
    ]
    for elem, x, y, z in coords:
        lines.append(f" {elem:2s} {x:15.10f} {y:15.10f} {z:15.10f}")

    return "\n".join(lines)
```

### 2. Output Parsing

```python
def parse_mopac_output(output_file: Path) -> dict:
    """Parse MOPAC output file."""
    results = {
        "converged": False,
        "hof": None,
        "homo": None,
        "lumo": None,
    }

    with open(output_file) as f:
        content = f.read()

    # Check convergence
    if "SCF CALCULATIONS DID NOT CONVERGE" in content:
        results["converged"] = False
        return results

    results["converged"] = True

    # Extract HOF (multiple patterns for robustness)
    patterns = [
        r"FINAL HEAT OF FORMATION\s*=\s*([-\d.]+)",
        r"HEAT OF FORMATION\s*=\s*([-\d.]+)",
        r"HOF\s*=\s*([-\d.]+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, content)
        if match:
            results["hof"] = float(match.group(1))
            break

    # Extract HOMO/LUMO
    homo_lumo = re.search(
        r"HOMO LUMO ENERGIES.*?=\s*([-\d.]+)\s+([-\d.]+)",
        content
    )
    if homo_lumo:
        results["homo"] = float(homo_lumo.group(1))
        results["lumo"] = float(homo_lumo.group(2))

    return results
```

### 3. Error Recovery

```python
def run_mopac_with_retry(mol_id: str, xyz_file: Path, config: PM7Config) -> MOPACResult:
    """Run MOPAC with retry logic."""

    # First attempt: standard
    result = run_mopac(mol_id, xyz_file, "PM7 1SCF")
    if result.converged:
        return result

    # Second attempt: with PULAY
    result = run_mopac(mol_id, xyz_file, "PM7 1SCF PULAY")
    if result.converged:
        return result

    # Third attempt: relaxed criteria
    result = run_mopac(mol_id, xyz_file, "PM7 1SCF PULAY ITRY=2000")
    if result.converged:
        result.quality_grade = QualityGrade.B  # Mark as lower confidence
        return result

    # Failed all attempts
    return MOPACResult(
        mol_id=mol_id,
        success=False,
        error_message="SCF failed after 3 attempts",
    )
```

---

## References

### Official Documentation

- [OpenMOPAC Website](https://openmopac.net/)
- [PM7 Method Description](https://openmopac.net/Manual/pm7.html)
- [MOPAC Keywords](https://openmopac.net/Manual/allkeys.html)
- [Error Messages](https://openmopac.net/Manual/error_messages.html)

### Key Publications

1. Stewart, J.J.P. *J. Mol. Model.* **2013**, 19, 1-32. "Optimization of
   parameters for semiempirical methods VI: more modifications to the NDDO
   approximations and re-optimization of parameters"

2. Stewart, J.J.P. *J. Mol. Model.* **2007**, 13, 1173-1213. "Optimization of
   parameters for semiempirical methods V: Modification of NDDO approximations
   and application to 70 elements"

3. Ryde, U.; Grimme, S. *Chem. Rev.* **2021**, 121, 12459-12521. "Quantum
   Mechanical Methods for Biomolecular Simulations"

### Related Grimperium Documentation

- [CREST Integration](./CREST_INTEGRATION.md) - Conformational sampling
- [Workflow Documentation](./WORKFLOW.md) - Full pipeline description

---

## Appendix: MOPAC Keyword Reference

### Method Keywords

```
PM7      - PM7 method (recommended)
PM6      - PM6 method (older)
AM1      - Austin Model 1
RM1      - Recife Model 1
```

### SCF Keywords

```
1SCF           - Single-point (no optimization)
SCFCRT=<E>     - SCF convergence criterion
ITRY=<N>       - Maximum SCF iterations
PULAY          - PULAY DIIS acceleration
CAMP           - Camp-King damping
```

### Output Keywords

```
AUX            - Write auxiliary file
AUX(0)         - Minimal auxiliary output
THERMO         - Thermodynamic analysis
BONDS          - Bond order analysis
ESP            - Electrostatic potential
GRAPH          - Generate graph data
```

### Optimization Keywords (not used in Phase A)

```
BFGS           - BFGS optimizer
EF             - Eigenvector following
TS             - Transition state search
GRADIENTS      - Calculate gradients
```

### Precision Keywords

```
PRECISE        - 100x tighter criteria
GNORM=<E>      - Gradient norm criterion
RELSCF=<E>     - Relative SCF convergence
```

---

**Last Updated:** 2026-01-15
**Version:** Phase A
**Grimperium Version:** 0.2.0
