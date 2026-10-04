# CREST atom-order fixtures

Three stored CREST runs, copied verbatim from `runs/` (not versioned), that
pin the assumption `SmilesTopology` rests on: CREST writes every conformer
in the atom order of its input, and that input keeps the RDKit `AddHs`
order of the SMILES it was embedded from.

## Files

Per run directory:

- `rdkit_initial.xyz` — RDKit embedding of the SMILES (explicit hydrogens).
- `input.xyz` — the xTB-optimised structure handed to CREST.
- `crest_conformers.xyz` — the final CREST ensemble.

`cases.json` holds, per run, a SMILES written in the source atom order
(non-canonical RDKit SMILES of the heavy atoms of `rdkit_initial.xyz`, bonds
perceived with `rdDetermineBonds`, charge 0) and the canonical SMILES of the
same molecule. The runs never stored their source SMILES, so the
source-order SMILES is reconstructed from the geometry.

SHA-256 of the copied files:

| run | file | sha256 |
|---|---|---|
| cbt_003 | rdkit_initial.xyz | `6f812c33b4869efabbd4dad79c76cf26c19263116d4ff5c1c7be7f646ce66fa4` |
| cbt_003 | input.xyz | `0582c666d025dde895a2dae472634179f661d51646f5036c086a399d05fb8efb` |
| cbt_003 | crest_conformers.xyz | `ff98af499ff25dea038f64f6b7446a62f68366b438bff56512a37588d4178609` |
| cbt_013 | rdkit_initial.xyz | `6cecac9e685b6664e4f59134cf50a3a4b28d6aaa5f53b7d536ad13ca7fb6d3d1` |
| cbt_013 | input.xyz | `1de3013ca791afc00944f4e4339425e09ce5b62c9b73cdac9a54cba011dc7ae8` |
| cbt_013 | crest_conformers.xyz | `b415542116f91dddc7d22bb8cf5dcc6a05a2dc5c45aa84d29352f156a7336285` |
| cbt_020 | rdkit_initial.xyz | `4f89502fc2ba8ee5a81c4ebbb4c8c3da951b0f9b1e0e778a25fdd7318ceccf15` |
| cbt_020 | input.xyz | `15beb96da871f560be05ffd96cf3f84692721b925d752a7738f44a7f9b3a704e` |
| cbt_020 | crest_conformers.xyz | `0423aec7545621eae9c0ffe1ffae3d5af6736606a4c68737d157149ac9b745f5` |

## Survey of all stored runs (2026-10-04)

All 37 `runs/cbt_*` directories were inspected. 35 have a complete
`crest_conformers.xyz` (1 to 13 778 conformers each, about 105 000 in
total); `cbt_032` only has a truncated `crest_ensemble.xyz`, which
`parse_crest_ensemble` refuses, and `cbt_033` has no ensemble. In all 35:

- `rdkit_initial.xyz`, `input.xyz`, `crest_input_copy.xyz` and every
  conformer have the same element at every position;
- every conformer has the same perceived connectivity as `rdkit_initial.xyz`,
  and `require_matching_order` accepts each one with that topology;
- `rdkit_initial.xyz` has the `AddHs` layout (heavy atoms first, then
  hydrogens grouped by ascending heavy parent).

A canonical SMILES re-parsed does **not** keep the source order in general:
it was refused for glycerol (`cbt_013`) and methyl acetate (`cbt_020`) and
accepted for the other 33. The source-order SMILES was accepted for all 35.

## Scope

These fixtures prove that CREST keeps its input order and that
`SmilesTopology` accepts a real CREST ensemble when the SMILES has the
source order. They do not prove that a given dataset SMILES has that order;
the topology check exists to catch exactly that at run time.
