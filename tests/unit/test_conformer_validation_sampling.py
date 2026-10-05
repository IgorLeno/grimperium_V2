"""Sampling of the conformer-selection validation experiment."""

from __future__ import annotations

import pandas as pd
import pytest

from tests.unit.conformer_validation_support import load

sampling = load("sampling")

# Butanoic acid, all heavy atoms bonded (perceives to its SMILES).
BUTANOIC_XYZ = """14
butanoic acid
C 2.4627 -0.4110 0.0000
C 1.2223 0.4870 0.0000
C -0.0681 -0.3365 0.0000
C -1.2996 0.5185 0.0000
O -1.3166 1.7310 0.0000
O -2.4362 -0.2190 0.0000
H 3.3640 0.2117 0.0000
H 2.4876 -1.0549 0.8860
H 2.4876 -1.0549 -0.8860
H 1.2307 1.1348 0.8820
H 1.2307 1.1348 -0.8820
H -0.1010 -0.9825 0.8800
H -0.1010 -0.9825 -0.8800
H -3.2205 0.3500 0.0000
"""
# Same atoms with the acid group pulled 4 A away: two fragments.
DISCONNECTED_XYZ = (
    BUTANOIC_XYZ.replace("-1.2996 0.5185", "-5.2996 0.5185")
    .replace("-1.3166 1.7310", "-5.3166 1.7310")
    .replace("-2.4362 -0.2190", "-6.4362 -0.2190")
    .replace("-3.2205 0.3500", "-7.2205 0.3500")
)


def test_bands_follow_the_approved_edges() -> None:
    nheavy = [
        sampling.band_of(n, sampling.NHEAVY_BANDS)
        for n in (1, 6, 7, 8, 9, 10, 11, 12, 22)
    ]
    assert nheavy == ["le6", "le6", "7-8", "7-8", "9", "10", "11", "ge12", "ge12"]
    rot = [sampling.band_of(n, sampling.ROTATABLE_BANDS) for n in (0, 1, 2, 3, 4, 9)]
    assert rot == ["0", "1", "2", "3", "ge4", "ge4"]


def test_rotatable_bonds_use_the_strict_rdkit_definition() -> None:
    assert sampling.strict_rotatable_bonds("CCCC(=O)O") == 2
    assert sampling.strict_rotatable_bonds("c1ccccc1") == 0
    # Strict mode does not count the amide C-N bond.
    assert sampling.strict_rotatable_bonds("CC(=O)NC") == 0


def test_reference_check_matches_the_calibration_statuses() -> None:
    assert sampling.reference_status(BUTANOIC_XYZ, "CCCC(=O)O") == "ok"
    assert sampling.reference_status(BUTANOIC_XYZ, "CC(C)C(=O)O") == "topology_mismatch"
    assert (
        sampling.reference_status(DISCONNECTED_XYZ, "CCCC(=O)O")
        == "disconnected_geometry"
    )
    assert sampling.reference_status("not xyz", "C") == "xyz_parse_failed"


def test_largest_remainder_is_exact_and_deterministic() -> None:
    assert sampling.largest_remainder(40, [284, 938, 980, 1367, 1348, 130]) == [
        2,
        7,
        8,
        11,
        11,
        1,
    ]
    assert sum(sampling.largest_remainder(40, [1, 1, 1])) == 40
    assert sampling.largest_remainder(2, [1, 1, 1]) == [1, 1, 0]


def test_population_keeps_closed_shell_neutral_one_row_per_smiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sampling, "reference_status", lambda xyz, smiles: xyz)
    frame = pd.DataFrame(
        [
            {
                "mol_id": "a2",
                "smiles": "CCO",
                "xyz": "ok",
                "multiplicity": 1,
                "charge": 0,
                "nheavy": 3,
                "conformer_rank_by_h298": 2,
            },
            {
                "mol_id": "a1",
                "smiles": "CCO",
                "xyz": "ok",
                "multiplicity": 1,
                "charge": 0,
                "nheavy": 3,
                "conformer_rank_by_h298": 1,
            },
            {
                "mol_id": "b",
                "smiles": "CC[O]",
                "xyz": "ok",
                "multiplicity": 2,
                "charge": 0,
                "nheavy": 3,
                "conformer_rank_by_h298": 1,
            },
            {
                "mol_id": "c",
                "smiles": "CC[NH3+]",
                "xyz": "ok",
                "multiplicity": 1,
                "charge": 1,
                "nheavy": 3,
                "conformer_rank_by_h298": 1,
            },
            {
                "mol_id": "d",
                "smiles": "CCC",
                "xyz": "topology_mismatch",
                "multiplicity": 1,
                "charge": 0,
                "nheavy": 3,
                "conformer_rank_by_h298": 1,
            },
        ]
    )

    usable, statuses = sampling.population(frame)

    assert list(usable["mol_id"]) == ["a1"]
    assert statuses == {"ok": 1, "topology_mismatch": 1}
    assert usable.loc[0, "stratum"] == "rot_0|nheavy_le6"


def _population_frame() -> pd.DataFrame:
    """Rows spread over every rotatable band and two nheavy bands."""
    families = {
        "0": ["C1CCCCC1", "C1CCCC1"],
        "1": ["CCC1CCCCC1", "CCC1CCCC1"],
        "2": ["CCCC1CCCCC1", "CCCC1CCCC1"],
        "3": ["CCCCC1CCCCC1", "CCCCC1CCCC1"],
        "ge4": ["CCCCCC1CCCCC1", "CCCCCC1CCCC1"],
    }
    rows = []
    for rot, bases in families.items():
        for base in bases:
            for n in range(12):
                rows.append(
                    {
                        "mol_id": f"{rot}-{base}-{n:02d}",
                        # Water fragments keep the SMILES distinct without rotors.
                        "smiles": base + ".O" * n,
                        "xyz": "ok",
                        "multiplicity": 1,
                        "charge": 0,
                        "nheavy": 9 if base.endswith("CCCCC1") else 10,
                        "conformer_rank_by_h298": 1,
                    }
                )
    return pd.DataFrame(rows)


def test_draw_is_deterministic_and_fills_equal_rotatable_quotas(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sampling, "reference_status", lambda xyz, smiles: xyz)
    frame = _population_frame()

    first, report = sampling.draw(frame, size=20, seed=7, dataset_sha256="abc")
    again, _ = sampling.draw(frame, size=20, seed=7, dataset_sha256="abc")
    other, _ = sampling.draw(frame, size=20, seed=8, dataset_sha256="abc")

    assert list(first.columns) == sampling.SAMPLE_COLUMNS
    assert first.equals(again)
    assert list(first["mol_id"]) != list(other["mol_id"])
    assert first["rotatable_band"].value_counts().to_dict() == {
        "0": 4,
        "1": 4,
        "2": 4,
        "3": 4,
        "ge4": 4,
    }
    assert first["mol_id"].is_unique
    assert set(first["seed"]) == {7}
    assert set(first["dataset_sha256"]) == {"abc"}
    assert sum(report.quotas.values()) == 20


def test_quota_larger_than_a_cell_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sampling, "reference_status", lambda xyz, smiles: xyz)
    with pytest.raises(ValueError, match="quota"):
        sampling.draw(_population_frame(), size=200, seed=1)


def test_size_must_split_over_rotatable_bands() -> None:
    with pytest.raises(ValueError, match="split evenly"):
        sampling.quotas(pd.DataFrame(columns=["rotatable_band", "nheavy_band"]), 21)
