#!/usr/bin/env python3
"""Calibrate the extreme-folding filter on the ThermoCBS CHON reference geometries.

Spec: ``docs/plans/2026-10-03-conformer-selection-spec.md``, section 4.1.

Why this script exists:

* the conformer-selection protocol will discard CREST conformers whose
  molecular ends touch (dispersion or H-bond folding without a new bond);
* the thresholds ``k`` (min topological distance), ``f`` (vdW scale) and
  ``c_min`` (min contacts) must be fixed *before* any new CREST/MOPAC run,
  using only existing data: the CBS-QB3 reference geometries are taken as
  physically reasonable conformers, so a good filter must keep almost all
  of them (spec budget: at most 1% discarded);
* ``rg_ratio`` needs real ensembles and is out of scope here.

Descriptor per geometry: ``fold_contacts`` = heavy-atom pairs with
topological distance >= k and spatial distance < f * (r_vdW,i + r_vdW,j),
skipping classic intramolecular H-bond pairs (donor O-H/N-H ... acceptor
O/N with H...acceptor < ``HBOND_MAX_H_ACCEPTOR``). Intramolecular H-bonds
occur in ~23.5% of the references, so they are not folding evidence.

Topology is perceived from the geometry (RDKit ``DetermineConnectivity``)
and checked against the row SMILES; production uses the resolved molecule,
so rows whose perceived heavy-atom graph disagrees are excluded and counted.

Outputs (never overwritten): ``grid.csv`` (discard rate per grid point),
``flagged.csv`` (rows discarded by the proposed setting), ``report.md``.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdDetermineBonds, rdMolHash

# Same H...acceptor cutoff as the exploratory H-bond measurement in the spec.
HBOND_MAX_H_ACCEPTOR = 2.5  # angstrom
HBOND_ELEMENTS = frozenset({"N", "O"})

DEFAULT_K = (4, 5, 6, 7)
DEFAULT_F = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10)
DEFAULT_C_MIN = (1, 2, 3, 4)
DEFAULT_MAX_DISCARD = 0.01
EXAMPLES_IN_REPORT = 25

# Proposal bounds, from the first full run (2026-10-03) on the references:
# k=4 flags 1,5 contacts forced by substitution (ortho-tBu anisoles, 16-19%
# of nheavy 15-16), k=5 flags 1,6 contacts across planar conjugated groups
# (oximes, enamines); only k>=6 flags chain ends meeting (diacids, glycol
# ethers). f > 1 means "near", not "in contact", so it is not *extreme* folding.
PROPOSAL_MIN_K = 6
PROPOSAL_MAX_F = 1.00

INPUT_COLUMNS = [
    "mol_id",
    "smiles",
    "xyz",
    "multiplicity",
    "charge",
    "nheavy",
    "conformer_rank_by_h298",
]


@dataclass(frozen=True)
class PairTable:
    """Heavy-atom pairs of one geometry that can count as fold contacts."""

    top_distance: np.ndarray  # int, bonds
    vdw_ratio: np.ndarray  # distance / (r_vdW,i + r_vdW,j)
    atoms: np.ndarray  # (n, 2) atom indices, for the report
    distance: np.ndarray  # angstrom


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mol_from_xyz(xyz_block: str) -> Chem.Mol:
    mol = Chem.MolFromXYZBlock(xyz_block)
    if mol is None:
        raise ValueError("unparseable xyz block")
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol


def topology_matches_smiles(mol: Chem.Mol, smiles: str) -> bool:
    """Compare heavy-atom element graphs, ignoring bond orders, charges and stereo."""
    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        return False
    # The perceived geometry graph carries no stereo tags or H atoms; the
    # SMILES may (e.g. ``[H]/N=C(/N)...`` pins imine stereo).
    Chem.RemoveStereochemistry(reference)
    reference = Chem.RemoveHs(reference)
    heavy = Chem.RemoveHs(mol, sanitize=False)
    hash_fn = rdMolHash.HashFunction.ElementGraph
    return bool(
        rdMolHash.MolHash(heavy, hash_fn) == rdMolHash.MolHash(reference, hash_fn)
    )


def hbond_pairs(mol: Chem.Mol, coords: np.ndarray) -> set[tuple[int, int]]:
    """Heavy-atom (donor, acceptor) pairs joined by a classic H-bond, as sorted tuples."""
    pairs: set[tuple[int, int]] = set()
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in HBOND_ELEMENTS]
    for donor in acceptors:
        hydrogens = [
            n.GetIdx()
            for n in mol.GetAtomWithIdx(donor).GetNeighbors()
            if n.GetAtomicNum() == 1
        ]
        for h_idx, acceptor in product(hydrogens, acceptors):
            if acceptor == donor:
                continue
            if np.linalg.norm(coords[h_idx] - coords[acceptor]) < HBOND_MAX_H_ACCEPTOR:
                pairs.add((min(donor, acceptor), max(donor, acceptor)))
    return pairs


def pair_table(mol: Chem.Mol, min_k: int) -> PairTable:
    coords = mol.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    topo = Chem.GetDistanceMatrix(mol)
    table = Chem.GetPeriodicTable()
    radii = {i: table.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in heavy}
    skip = hbond_pairs(mol, coords)

    rows: list[tuple[int, int, int, float, float]] = []
    for pos, i in enumerate(heavy):
        for j in heavy[pos + 1 :]:
            bonds = int(topo[i, j])
            if bonds < min_k or (i, j) in skip:
                continue
            distance = float(np.linalg.norm(coords[i] - coords[j]))
            rows.append((i, j, bonds, distance, distance / (radii[i] + radii[j])))

    if not rows:
        empty = np.empty(0)
        return PairTable(empty.astype(int), empty, np.empty((0, 2), int), empty)
    arr = np.array(rows, dtype=float)
    return PairTable(
        top_distance=arr[:, 2].astype(int),
        vdw_ratio=arr[:, 4],
        atoms=arr[:, :2].astype(int),
        distance=arr[:, 3],
    )


def perceive(frame: pd.DataFrame, min_k: int) -> tuple[pd.DataFrame, list[PairTable]]:
    """Return per-row status and the pair tables of usable rows (same order)."""
    status: list[str] = []
    tables: list[PairTable] = []
    for xyz, smiles in zip(frame["xyz"], frame["smiles"], strict=True):
        try:
            mol = mol_from_xyz(xyz)
        except Exception:  # noqa: BLE001 - RDKit raises several types; counted below
            status.append("xyz_parse_failed")
            continue
        if len(Chem.GetMolFrags(mol)) != 1:
            status.append("disconnected_geometry")
        elif not topology_matches_smiles(mol, smiles):
            status.append("topology_mismatch")
        else:
            status.append("ok")
            tables.append(pair_table(mol, min_k))
    out = frame.drop(columns=["xyz"]).copy()
    out["perception"] = status
    return out, tables


def contact_counts(tables: list[PairTable], k: int, f: float) -> np.ndarray:
    return np.array(
        [
            int(np.count_nonzero((t.top_distance >= k) & (t.vdw_ratio < f)))
            for t in tables
        ]
    )


def evaluate_grid(
    usable: pd.DataFrame,
    tables: list[PairTable],
    ks: tuple[int, ...],
    fs: tuple[float, ...],
    c_mins: tuple[int, ...],
) -> pd.DataFrame:
    # Validation population of the spec (section 5): closed-shell neutral.
    closed_shell = ((usable["multiplicity"] == 1) & (usable["charge"] == 0)).to_numpy()
    records = []
    for k, f in product(ks, fs):
        counts = contact_counts(tables, k, f)
        for c_min in c_mins:
            flagged = counts >= c_min
            records.append(
                {
                    "k": k,
                    "f": f,
                    "c_min": c_min,
                    "flagged": int(flagged.sum()),
                    "discard_rate": float(flagged.mean()),
                    "discard_rate_closed_shell": float(flagged[closed_shell].mean()),
                }
            )
    return pd.DataFrame.from_records(records)


def propose(grid: pd.DataFrame, max_discard: float) -> pd.Series:
    """Pick the setting that detects the most contacts while keeping the budget.

    Only ``k >= PROPOSAL_MIN_K`` and ``f <= PROPOSAL_MAX_F`` qualify (see the
    constants). Within budget, prefer the largest ``f`` (so real CREST folds
    are caught), then the smallest ``k`` (shorter loops count as folds), then
    the smallest ``c_min``. Higher reference rate breaks remaining ties.
    """
    within = grid[
        (grid["discard_rate"] <= max_discard)
        & (grid["k"] >= PROPOSAL_MIN_K)
        & (grid["f"] <= PROPOSAL_MAX_F)
    ]
    if within.empty:
        raise ValueError(
            f"no grid point with k >= {PROPOSAL_MIN_K}, f <= {PROPOSAL_MAX_F} "
            f"keeps discard_rate <= {max_discard}"
        )
    ordered = within.sort_values(
        ["f", "k", "c_min", "discard_rate"], ascending=[False, True, True, False]
    )
    return ordered.iloc[0]


def frontier(grid: pd.DataFrame, max_discard: float) -> pd.DataFrame:
    """Largest f within budget for each (k, c_min)."""
    within = grid[grid["discard_rate"] <= max_discard]
    idx = within.groupby(["k", "c_min"])["f"].idxmax()
    return within.loc[idx].sort_values(["k", "c_min"]).reset_index(drop=True)


def flagged_rows(
    usable: pd.DataFrame, tables: list[PairTable], k: int, f: float, c_min: int
) -> pd.DataFrame:
    records = []
    for row, table in zip(usable.itertuples(index=False), tables, strict=True):
        mask = (table.top_distance >= k) & (table.vdw_ratio < f)
        if int(mask.sum()) < c_min:
            continue
        worst = int(np.argmin(np.where(mask, table.vdw_ratio, np.inf)))
        i, j = table.atoms[worst]
        records.append(
            {
                "mol_id": row.mol_id,
                "smiles": row.smiles,
                "nheavy": row.nheavy,
                "multiplicity": row.multiplicity,
                "fold_contacts": int(mask.sum()),
                "closest_pair": f"{i}-{j}",
                "closest_pair_bonds": int(table.top_distance[worst]),
                "closest_pair_distance_A": round(float(table.distance[worst]), 3),
                "closest_pair_vdw_ratio": round(float(table.vdw_ratio[worst]), 3),
            }
        )
    columns = [
        "mol_id",
        "smiles",
        "nheavy",
        "multiplicity",
        "fold_contacts",
        "closest_pair",
        "closest_pair_bonds",
        "closest_pair_distance_A",
        "closest_pair_vdw_ratio",
    ]
    frame = pd.DataFrame.from_records(records, columns=columns)
    return frame.sort_values(
        ["fold_contacts", "closest_pair_vdw_ratio"], ascending=[False, True]
    )


def percent(value: float) -> str:
    return f"{100 * value:.2f}%"


def render_report(
    *,
    source: Path,
    source_sha: str,
    perceived: pd.DataFrame,
    grid: pd.DataFrame,
    best: pd.Series,
    flagged: pd.DataFrame,
    usable: pd.DataFrame,
    max_discard: float,
) -> str:
    k, f, c_min = int(best["k"]), float(best["f"]), int(best["c_min"])
    lines = [
        "# Calibração do filtro de dobramento extremo",
        "",
        f"Gerado em {datetime.now(timezone.utc).isoformat()} por "
        "`scripts/calibrate_folding_filter.py`.",
        "Especificação: `docs/plans/2026-10-03-conformer-selection-spec.md`, seção 4.1.",
        "",
        "## Entrada",
        "",
        f"- Arquivo: `{source}`",
        f"- SHA-256: `{source_sha}`",
        f"- Linhas: {len(perceived)}",
        "",
        "Percepção de topologia (geometria vs SMILES da linha):",
        "",
        "| status | linhas |",
        "|---|---|",
    ]
    for status, count in perceived["perception"].value_counts().items():
        lines.append(f"| {status} | {count} |")
    lines += [
        "",
        f"Calibração usa apenas `ok` ({len(usable)} linhas).",
        "",
        "## Definição",
        "",
        "- `fold_contacts`: pares de átomos pesados com distância topológica ≥ k "
        "e distância < f · (r_vdW,i + r_vdW,j).",
        "- Raios de van der Waals: tabela periódica do RDKit "
        f"({Chem.rdBase.rdkitVersion}): C 1,70; N 1,60; O 1,55 Å.",
        "- Excluídos: pares doador/aceitador O/N com H do doador a "
        f"< {HBOND_MAX_H_ACCEPTOR} Å do aceitador (ligação H clássica).",
        "- Descartada se `fold_contacts ≥ c_min`. `rg_ratio` fora desta calibração.",
        f"- Orçamento: ≤ {percent(max_discard)} das referências descartadas.",
        "- Coluna `fechada` = multiplicidade 1 e carga 0 (população da validação).",
        "",
        "## Proposta",
        "",
        f"**k = {k}, f = {f:.2f}, c_min = {c_min}** — descarta "
        f"{int(best['flagged'])} referências ({percent(best['discard_rate'])}; "
        f"camada fechada {percent(best['discard_rate_closed_shell'])}).",
        "",
        f"Regra de escolha (no código, `propose`): só k ≥ {PROPOSAL_MIN_K} e "
        f"f ≤ {PROPOSAL_MAX_F:.2f}; dentro do orçamento, maior f, depois menor k, "
        "depois menor c_min. k = 4 sinaliza contatos 1,5 forçados por "
        "substituição (anisóis com tBu orto) e k = 5 contatos 1,6 em grupos "
        "conjugados planos (oximas, enaminas); só k ≥ 6 sinaliza pontas de cadeia "
        "se encontrando. f > 1 é proximidade, não contato extremo.",
        "",
        "Taxa de descarte por `nheavy` com a proposta:",
        "",
        "| nheavy | linhas | descartadas | taxa |",
        "|---|---|---|---|",
    ]
    flagged_ids = set(flagged["mol_id"])
    by_size = usable.assign(flag=usable["mol_id"].isin(flagged_ids)).groupby("nheavy")
    for nheavy, group in by_size:
        lines.append(
            f"| {nheavy} | {len(group)} | {int(group['flag'].sum())} | "
            f"{percent(float(group['flag'].mean()))} |"
        )
    lines += [
        "",
        "## Fronteira (maior f dentro do orçamento por k, c_min)",
        "",
        "| k | c_min | f | taxa | fechada |",
        "|---|---|---|---|---|",
    ]
    for row in frontier(grid, max_discard).itertuples(index=False):
        lines.append(
            f"| {row.k} | {row.c_min} | {row.f:.2f} | {percent(row.discard_rate)} | "
            f"{percent(row.discard_rate_closed_shell)} |"
        )
    c_values = sorted(grid["c_min"].unique())
    header = " | ".join(f"c_min={c}" for c in c_values)
    lines += [
        "",
        "## Grade completa (taxa de descarte, todas as linhas `ok`)",
        "",
        f"| k | f | {header} |",
        "|---|---|" + "---|" * len(c_values),
    ]
    pivot = grid.pivot_table(index=["k", "f"], columns="c_min", values="discard_rate")
    for (k_val, f_val), rates in pivot.iterrows():
        cells = " | ".join(percent(float(rates[c])) for c in c_values)
        lines.append(f"| {k_val} | {f_val:.2f} | {cells} |")
    lines += [
        "",
        f"## Exemplos sinalizados pela proposta (top {EXAMPLES_IN_REPORT})",
        "",
        "Ordem: mais contatos, depois menor razão vdW. Lista completa em `flagged.csv`.",
        "Índices de átomo são 0-based na ordem do bloco `xyz`.",
        "",
        "| mol_id | SMILES | nheavy | contatos | par | ligações | d (Å) | d/Σr_vdW |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in flagged.head(EXAMPLES_IN_REPORT).itertuples(index=False):
        lines.append(
            f"| {row.mol_id} | `{row.smiles}` | {row.nheavy} | {row.fold_contacts} | "
            f"{row.closest_pair} | {row.closest_pair_bonds} | "
            f"{row.closest_pair_distance_A:.3f} | {row.closest_pair_vdw_ratio:.3f} |"
        )
    lines += [
        "",
        "## Limites",
        "",
        "- Só referências CBS-QB3 (uma geometria por linha); não mede quantos "
        "conformeros CREST dobrados o filtro pegaria.",
        "- `r_max` (`rg_ratio`) exige ensembles reais; calibrado na validação.",
        "- Topologia percebida da geometria; em produção vem da molécula resolvida.",
        "",
    ]
    return "\n".join(lines)


def parse_floats(text: str) -> tuple[float, ...]:
    return tuple(float(v) for v in text.split(","))


def parse_ints(text: str) -> tuple[int, ...]:
    return tuple(int(v) for v in text.split(","))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--input", type=Path, required=True, help="thermo_cbs_chon_v2.csv"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--k", type=parse_ints, default=DEFAULT_K)
    parser.add_argument("--f", type=parse_floats, default=DEFAULT_F)
    parser.add_argument("--c-min", type=parse_ints, default=DEFAULT_C_MIN)
    parser.add_argument("--max-discard", type=float, default=DEFAULT_MAX_DISCARD)
    parser.add_argument("--limit", type=int, default=None, help="first N rows only")
    args = parser.parse_args(argv)

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        print(f"refusing to write into non-empty {args.output_dir}", file=sys.stderr)
        return 1

    rdBase.DisableLog("rdApp.*")
    frame = pd.read_csv(args.input, usecols=INPUT_COLUMNS, nrows=args.limit)
    perceived, tables = perceive(frame, min(args.k))
    usable = perceived[perceived["perception"] == "ok"].reset_index(drop=True)

    grid = evaluate_grid(usable, tables, args.k, args.f, args.c_min)
    best = propose(grid, args.max_discard)
    flagged = flagged_rows(
        usable, tables, int(best["k"]), float(best["f"]), int(best["c_min"])
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    grid.to_csv(args.output_dir / "grid.csv", index=False)
    flagged.to_csv(args.output_dir / "flagged.csv", index=False)
    report = render_report(
        source=args.input,
        source_sha=sha256_of(args.input),
        perceived=perceived,
        grid=grid,
        best=best,
        flagged=flagged,
        usable=usable,
        max_discard=args.max_discard,
    )
    (args.output_dir / "report.md").write_text(report, encoding="utf-8")
    print(
        f"usable={len(usable)}/{len(perceived)} proposal: k={int(best['k'])} "
        f"f={float(best['f']):.2f} c_min={int(best['c_min'])} "
        f"discard={percent(float(best['discard_rate']))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
