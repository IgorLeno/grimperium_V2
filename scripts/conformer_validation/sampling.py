"""Stratified, reproducible sample of the CBS reference dataset.

Population: closed-shell neutral rows (multiplicity 1, charge 0), one row
per SMILES (lowest ``conformer_rank_by_h298``, then lowest ``mol_id``),
whose reference geometry perceives to the SMILES heavy-atom graph (the
same check as ``scripts/calibrate_folding_filter.py``).

Strata: 6 ``nheavy`` bands x 5 strict rotatable-bond bands. Every
rotatable band gets the same quota; inside a band the quota is split
across ``nheavy`` bands in proportion to the population (largest
remainder). Rows are drawn per cell with one seeded generator, cells in a
fixed order, so the same dataset and seed always give the same sample.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors, rdMolHash

NHEAVY_BANDS: tuple[tuple[str, int, int], ...] = (
    ("le6", 0, 6),
    ("7-8", 7, 8),
    ("9", 9, 9),
    ("10", 10, 10),
    ("11", 11, 11),
    ("ge12", 12, 10**6),
)
ROTATABLE_BANDS: tuple[tuple[str, int, int], ...] = (
    ("0", 0, 0),
    ("1", 1, 1),
    ("2", 2, 2),
    ("3", 3, 3),
    ("ge4", 4, 10**6),
)
DEFAULT_SIZE = 200
DEFAULT_SEED = 20261004

INPUT_COLUMNS = [
    "mol_id",
    "smiles",
    "xyz",
    "multiplicity",
    "charge",
    "nheavy",
    "conformer_rank_by_h298",
]
SAMPLE_COLUMNS = [
    "mol_id",
    "smiles",
    "nheavy",
    "rotatable_bonds",
    "nheavy_band",
    "rotatable_band",
    "stratum",
    "seed",
    "dataset_sha256",
]


@dataclass(frozen=True)
class SampleReport:
    """Counts behind one sample, printed by the CLI and kept in the docs."""

    rows_total: int
    closed_shell: int
    unique_smiles: int
    reference_status: dict[str, int]
    population: int
    cell_population: dict[str, int]
    quotas: dict[str, int]


def sha256_of(path: Path) -> str:
    """Return the SHA-256 of a file, streamed."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def band_of(value: int, bands: Sequence[tuple[str, int, int]]) -> str:
    """Return the label of the band that contains ``value``."""
    for label, low, high in bands:
        if low <= value <= high:
            return label
    raise ValueError(f"{value} falls outside every band")


def strict_rotatable_bonds(smiles: str) -> int:
    """RDKit strict rotatable-bond count of ``smiles``."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES {smiles!r}")
    return int(rdMolDescriptors.CalcNumRotatableBonds(molecule, True))


def reference_status(xyz_block: str, smiles: str) -> str:
    """Whether the reference geometry perceives to the SMILES graph.

    Returns ``ok``, ``xyz_parse_failed``, ``disconnected_geometry`` or
    ``topology_mismatch`` — the statuses of the folding calibration.
    """
    try:
        geometry = Chem.MolFromXYZBlock(xyz_block)
        if geometry is None:
            return "xyz_parse_failed"
        rdDetermineBonds.DetermineConnectivity(geometry)
    except Exception:  # noqa: BLE001 - RDKit raises several types
        return "xyz_parse_failed"
    if len(Chem.GetMolFrags(geometry)) != 1:
        return "disconnected_geometry"
    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        return "topology_mismatch"
    Chem.RemoveStereochemistry(reference)
    reference = Chem.RemoveHs(reference)
    heavy = Chem.RemoveHs(geometry, sanitize=False)
    hash_fn = rdMolHash.HashFunction.ElementGraph
    if rdMolHash.MolHash(heavy, hash_fn) != rdMolHash.MolHash(reference, hash_fn):
        return "topology_mismatch"
    return "ok"


def population(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Filter ``frame`` to the sampling population and annotate its strata.

    Returns the population and the count of each reference status among
    the one-per-SMILES candidates (``ok`` rows are the population).
    """
    closed = frame[(frame["multiplicity"] == 1) & (frame["charge"] == 0)]
    unique = closed.sort_values(["conformer_rank_by_h298", "mol_id"]).drop_duplicates(
        "smiles", keep="first"
    )
    RDLogger.DisableLog("rdApp.*")
    try:
        statuses = [
            reference_status(xyz, smiles)
            for xyz, smiles in zip(unique["xyz"], unique["smiles"], strict=True)
        ]
    finally:
        RDLogger.EnableLog("rdApp.*")
    unique = unique.assign(reference_status=statuses)
    counts = {
        str(key): int(value)
        for key, value in unique["reference_status"].value_counts().sort_index().items()
    }
    usable = unique[unique["reference_status"] == "ok"].copy()
    usable["rotatable_bonds"] = [strict_rotatable_bonds(s) for s in usable["smiles"]]
    usable["nheavy_band"] = [band_of(int(n), NHEAVY_BANDS) for n in usable["nheavy"]]
    usable["rotatable_band"] = [
        band_of(int(n), ROTATABLE_BANDS) for n in usable["rotatable_bonds"]
    ]
    usable["stratum"] = usable["rotatable_band"].map(lambda r: f"rot_{r}") + (
        "|nheavy_" + usable["nheavy_band"]
    )
    return usable.sort_values("mol_id").reset_index(drop=True), counts


def largest_remainder(total: int, weights: Sequence[int]) -> list[int]:
    """Split ``total`` proportionally to ``weights``; ties go to lower index."""
    weight_sum = sum(weights)
    if weight_sum == 0:
        raise ValueError("Cannot split a quota over empty cells")
    exact = [total * w / weight_sum for w in weights]
    shares = [int(np.floor(value)) for value in exact]
    order = sorted(range(len(weights)), key=lambda i: (-(exact[i] - shares[i]), i))
    for i in order[: total - sum(shares)]:
        shares[i] += 1
    return shares


def quotas(usable: pd.DataFrame, size: int) -> dict[tuple[str, str], int]:
    """Quota of every (rotatable band, nheavy band) cell.

    Raises:
        ValueError: If ``size`` does not split evenly over the rotatable
            bands or a cell holds fewer rows than its quota.
    """
    if size % len(ROTATABLE_BANDS):
        raise ValueError(
            f"Sample size {size} does not split evenly over "
            f"{len(ROTATABLE_BANDS)} rotatable-bond bands"
        )
    per_band = size // len(ROTATABLE_BANDS)
    counts = usable.groupby(["rotatable_band", "nheavy_band"]).size()
    result: dict[tuple[str, str], int] = {}
    for rot, _, _ in ROTATABLE_BANDS:
        cells = [int(counts.get((rot, heavy), 0)) for heavy, _, _ in NHEAVY_BANDS]
        for (heavy, _, _), available, share in zip(
            NHEAVY_BANDS, cells, largest_remainder(per_band, cells), strict=True
        ):
            if share > available:
                raise ValueError(
                    f"Cell rot_{rot}|nheavy_{heavy} holds {available} rows but "
                    f"its quota is {share}"
                )
            result[(rot, heavy)] = share
    return result


def draw(
    frame: pd.DataFrame,
    *,
    size: int = DEFAULT_SIZE,
    seed: int = DEFAULT_SEED,
    dataset_sha256: str = "",
) -> tuple[pd.DataFrame, SampleReport]:
    """Draw the stratified sample from the raw dataset ``frame``."""
    usable, statuses = population(frame)
    cell_quota = quotas(usable, size)
    rng = np.random.default_rng(seed)
    picked: list[pd.DataFrame] = []
    for rot, _, _ in ROTATABLE_BANDS:
        for heavy, _, _ in NHEAVY_BANDS:
            share = cell_quota[(rot, heavy)]
            if share == 0:
                continue
            cell = usable[
                (usable["rotatable_band"] == rot) & (usable["nheavy_band"] == heavy)
            ]
            rows = rng.choice(len(cell), size=share, replace=False)
            picked.append(cell.iloc[np.sort(rows)])
    sample = pd.concat(picked).sort_values("mol_id").reset_index(drop=True)
    sample["seed"] = seed
    sample["dataset_sha256"] = dataset_sha256
    closed = frame[(frame["multiplicity"] == 1) & (frame["charge"] == 0)]
    report = SampleReport(
        rows_total=len(frame),
        closed_shell=len(closed),
        unique_smiles=int(closed["smiles"].nunique()),
        reference_status=statuses,
        population=len(usable),
        cell_population={
            f"rot_{rot}|nheavy_{heavy}": int(count)
            for (rot, heavy), count in usable.groupby(["rotatable_band", "nheavy_band"])
            .size()
            .items()
        },
        quotas={
            f"rot_{rot}|nheavy_{heavy}": share
            for (rot, heavy), share in cell_quota.items()
        },
    )
    return sample[SAMPLE_COLUMNS], report


def load_dataset(path: Path) -> pd.DataFrame:
    """Read the columns sampling needs from the CBS dataset."""
    return pd.read_csv(path, usecols=INPUT_COLUMNS)


def render_report(report: SampleReport, sample: pd.DataFrame) -> str:
    """Human-readable account of how the sample was drawn."""
    lines = [
        f"dataset rows: {report.rows_total}",
        f"closed-shell neutral rows: {report.closed_shell}",
        f"unique SMILES: {report.unique_smiles}",
        "reference check (one row per SMILES): "
        + ", ".join(f"{k}={v}" for k, v in report.reference_status.items()),
        f"population: {report.population}",
        "cells (population -> quota):",
    ]
    for cell, share in report.quotas.items():
        lines.append(f"  {cell}: {report.cell_population.get(cell, 0)} -> {share}")
    lines.append(f"sample size: {len(sample)}")
    return "\n".join(lines)


__all__ = [
    "DEFAULT_SEED",
    "DEFAULT_SIZE",
    "NHEAVY_BANDS",
    "ROTATABLE_BANDS",
    "SAMPLE_COLUMNS",
    "SampleReport",
    "band_of",
    "draw",
    "largest_remainder",
    "load_dataset",
    "population",
    "quotas",
    "reference_status",
    "render_report",
    "sha256_of",
    "strict_rotatable_bonds",
]
