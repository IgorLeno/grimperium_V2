from __future__ import annotations

from pathlib import Path

from semi_imperium.training.geometry_store import GeometryRecord, GeometryStore


def _record(mol_id: str = "mol_00001", hof: float = -12.5) -> GeometryRecord:
    return GeometryRecord(
        mol_id=mol_id,
        smiles="CO",
        hamiltonian="PM7",
        protocol="2",
        signature="abc123",
        calculation_id="calc-1",
        conformer_index=3,
        elements=("C", "O", "H"),
        coordinates=((0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (-0.5, 0.9, 0.0)),
        hof_kcal_mol=hof,
        lowest_frequency_cm1=212.4,
        minimum_confirmed=True,
    )


def test_reading_a_missing_database_does_not_create_it(tmp_path: Path) -> None:
    path = tmp_path / "training.sqlite3"
    store = GeometryStore(path)
    assert store.count() == 0
    assert store.get("mol_00001") is None
    assert store.mol_ids() == set()
    assert not path.exists()


def test_save_then_get_round_trips(tmp_path: Path) -> None:
    store = GeometryStore(tmp_path / "training.sqlite3")
    store.save(_record())
    loaded = store.get("mol_00001")
    assert loaded is not None
    assert loaded.elements == ("C", "O", "H")
    assert loaded.coordinates[1] == (1.4, 0.0, 0.0)
    assert loaded.hof_kcal_mol == -12.5
    assert loaded.minimum_confirmed is True
    assert loaded.created_at  # filled on save


def test_saving_the_same_molecule_replaces_the_row(tmp_path: Path) -> None:
    store = GeometryStore(tmp_path / "training.sqlite3")
    store.save(_record(hof=-1.0))
    store.save(_record(hof=-2.0))
    assert store.count() == 1
    loaded = store.get("mol_00001")
    assert loaded is not None and loaded.hof_kcal_mol == -2.0


def test_xyz_export_keeps_atom_order_and_count(tmp_path: Path) -> None:
    lines = _record().to_xyz().splitlines()
    assert lines[0] == "3"
    assert "mol_00001" in lines[1]
    assert [line.split()[0] for line in lines[2:]] == ["C", "O", "H"]


def test_mol_ids_lists_every_saved_molecule(tmp_path: Path) -> None:
    store = GeometryStore(tmp_path / "training.sqlite3")
    store.save(_record("mol_00001"))
    store.save(_record("mol_00002"))
    assert store.mol_ids() == {"mol_00001", "mol_00002"}
