# grimperium_mini

`grimperium_mini` is a distilled workflow for the biodiesel TCC comparison. It
does not replace the original Grimperium package and does not include
delta-learning, descriptor extraction for ML, CBS-QB3 integration, dashboards,
or distributed batch state.

## Scientific Scope

The original TCC used VEGA for conformational search and MOPAC for AM1, PM3,
and PM7 geometry optimization and heat-of-formation calculations. Group
contribution values were calculated externally in spreadsheets.

`grimperium_mini` keeps the same molecule/reaction comparison frame but changes
the conformational search step:

```text
TCC antigo:      SMILES/estrutura -> VEGA -> MOPAC AM1/PM3/PM7
grimperium_mini: SMILES -> RDKit 3D seed -> CREST -> MOPAC AM1/PM3/PM7
Grimperium full: CREST/PM7 + descriptors + data fusion + ML delta-learning
```

Group contribution values are loaded only as reference columns. They are not
reimplemented or recalculated.

## Commands

```bash
python -m grimperium_mini validate-data --xlsx data/grimperium_mini_pipeline_tcc.xlsx
```

```bash
python -m grimperium_mini run \
  --xlsx data/grimperium_mini_pipeline_tcc.xlsx \
  --methods AM1 PM3 PM7 \
  --work-root runs \
  --results-dir results
```

```bash
python -m grimperium_mini run \
  --xlsx data/grimperium_mini_pipeline_tcc.xlsx \
  --methods PM7 \
  --limit 2 \
  --dry-run
```

```bash
python -m grimperium_mini export-summary \
  --input results/formacao_grimperium_mini.csv \
  --reactions results/reacoes_grimperium_mini.csv \
  --output results/grimperium_mini_summary.xlsx
```

## Outputs

- `results/formacao_grimperium_mini.csv`: one selected conformer per molecule
  and semiempirical method.
- `results/conformers_detail.csv`: one row per MOPAC-optimized conformer.
- `results/reacoes_grimperium_mini.csv`: reaction enthalpies recalculated from
  the new CREST+MOPAC formation enthalpies.
- `results/grimperium_mini_summary.xlsx`: optional summary workbook generated
  from the CSV outputs.

Existing CSV outputs are backed up with a timestamped `.bak.csv` name before
being rewritten.

## Configuration

Defaults live in `grimperium_mini.config.MiniConfig`. The main external tools
can be set through CLI flags or environment variables:

- `GRIMPERIUM_MINI_CREST`
- `GRIMPERIUM_MINI_MOPAC`
- `GRIMPERIUM_MINI_XTB` — path to the xTB binary (default: `xtb`)
- `GRIMPERIUM_MINI_XTB_ENABLED` — set to `false` to skip xTB pre-optimization
- `GRIMPERIUM_MINI_THREADS`
- `GRIMPERIUM_MINI_WORK_ROOT`
- `GRIMPERIUM_MINI_RESULTS_DIR`

The MOPAC method is passed explicitly for every task and must be one of `AM1`,
`PM3`, or `PM7`.

## Known Issues

### Fedora xtb 6.7.1-4 optimizer crash

The Fedora-packaged `xtb 6.7.1-4.fc43` has a Fortran format-string bug in
`optimizer.f90:852` that causes `xtb --opt` to abort with `rc=2` and the
message:

```
Fortran runtime error: Missing comma between descriptors
(1x,""(""f7.2""%)"")
```

`grimperium-mini` detects this automatically at startup via a preflight check.
When the check fails, xTB pre-optimization is silently disabled for that run
and the pipeline falls back to `RDKit → CREST → MOPAC` (the same path used
before xTB support was added). A `WARNING` message is printed citing this issue.

**Permanent workarounds:**

```bash
# Option 1 — skip xTB via CLI flag
python -m grimperium_mini run --no-xtb --xlsx ...

# Option 2 — skip xTB via environment variable
export GRIMPERIUM_MINI_XTB_ENABLED=false

# Option 3 — install a working xTB (conda-forge build is unaffected)
conda install -c conda-forge xtb=6.7.1
export GRIMPERIUM_MINI_XTB=/path/to/conda/bin/xtb
```
