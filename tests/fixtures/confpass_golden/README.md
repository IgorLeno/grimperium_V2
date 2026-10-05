# CONFPASS golden files

Reference output of the original CONFPASS PART 1 (Goodman-lab/CONFPASS,
MIT, commit `1b5efb69585ea1f51bedccfed1d9d07133c18a53`) on 14 CREST
ensembles. The ported backend must reproduce `golden.json` on the SDF
files in `inputs/`.

## Files

- `inputs/*.sdf` — one file per case, written by the production adapter
  (`build_confpass_candidates`) from a stored `runs/cbt_*/crest_conformers.xyz`
  (`amidino_alanine`: `runs/conformer_validation/cbs_01102`)
  parsed with `parse_crest_ensemble`. Records keep CREST's order.
- `inputs.json` — per case: source run, SMILES perceived from the geometry
  (RDKit `DetermineBonds`, charge 0), atom and conformer counts, energy
  inversions in CREST's order, SHA-256 of the source XYZ and of the SDF.
- `golden.json` — CONFPASS commit and source hashes, environment, and per
  case: `pipe_x_as` / `pipe_x` / `pipe_as` priority lists (0-based SDF
  positions; the CLI prints them 1-based), the kept dihedral columns
  (sorted), the number of clustering rows, the clusters at
  `n = round(0.8 · conformers)`, or the exception CONFPASS raised.

## Regenerating

```bash
poetry run python scripts/confpass_golden/build_inputs.py --runs runs --output tests/fixtures/confpass_golden
```

```bash
PYTHONHASHSEED=0 ~/Estagio/confpass-original/.venv/bin/python scripts/confpass_golden/run_original.py --confpass-src ~/Estagio/confpass-original/src/confpass --fixtures tests/fixtures/confpass_golden
```

`runs/` is not versioned, so the SDF inputs are the source of truth here;
`inputs.json` keeps the hash of each original XYZ.

Isolated environment (outside the project): Python 3.8.20 via `uv`,
numpy 1.20.0, pandas 1.3.4, scikit-learn 1.0.1, natsort 8.0.2,
rdkit-pypi 2020.09.5. CONFPASS was developed with RDKit 2019.09.3 and
numpy 1.21.2; 2019.09 has no PyPI wheel, and rdkit-pypi 2020.09.5 pins
numpy 1.20.0. Only the PART 1 modules are imported; the PART 2 pickle
model was not downloaded.

## Observed behaviour the port must match or deliberately change

- **Energies are never read.** Priority `pipe_as` takes the first member of
  each cluster, i.e. file order is treated as energy order. Four inputs
  carry small inversions in CREST's order (largest 0.028 kcal/mol,
  triacetin); the golden files reflect CREST's order as written.
- **No variable dihedral ⇒ failure.** `methyl_acetate`: the only rotatable
  dihedral (ester C–O) is constant across the ensemble and removed, the
  clustering frame is empty, and every method raises (`ValueError` for
  `pipe_x_as`/`pipe_x`, `IndexError` for `pipe_as`). The port needs an
  explicit fallback instead.
- **Dihedral column order depends on `PYTHONHASHSEED`** (built through
  `set` iteration). Over 20 seeds the priority lists and clusters of the
  first 13 cases were identical; only column order changed, so the columns
  are stored sorted, and those cases are byte-identical across seeds.
- **Equivalent terminal hydrogens follow the seed.** `amidino_alanine`: the
  NH2 dihedral end is taken from `set(neighbours) - {partner}`, so CONFPASS
  picks H 11 or H 12 depending on `PYTHONHASHSEED` (H 11 in 13 of 20 seeds;
  the priorities change from position 13 on). The port always takes the
  lowest atom index; the golden file is generated with `PYTHONHASHSEED=0`,
  which picks the same hydrogen.
- **Stereo-defining hydrogen.** `amidino_alanine` also has an imine N-H.
  Current RDKit keeps it when reading the molblock without hydrogens
  (C=N stereo perceived from 3D); RDKit 2020.09.5 drops it. The port drops
  it too, otherwise the heavy-atom naming runs out of names.
