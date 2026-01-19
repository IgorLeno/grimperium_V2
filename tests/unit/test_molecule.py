"""
Tests for Molecule dataclass hierarchy.
"""

from datetime import datetime

import pytest

from grimperium.core.molecule import (
    Molecule,
    MoleculeIdentity,
    MoleculeMeta,
    MoleculeProperties,
    MoleculeResults,
    MoleculeStatus,
)


class TestMoleculeCreation:
    """Tests for Molecule creation."""

    def test_minimal_creation(self):
        """Create molecule with minimal fields."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="CCO"),
            properties=MoleculeProperties(nheavy=3),
        )

        assert mol.mol_id == "mol_001"
        assert mol.smiles == "CCO"
        assert mol.properties.nheavy == 3

    def test_default_status_is_pending(self):
        """Default status is PENDING."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
        )

        assert mol.status == MoleculeStatus.PENDING

    def test_default_reruns_is_zero(self):
        """Default reruns is 0."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
        )

        assert mol.reruns == 0


class TestMoleculeStatusChecks:
    """Tests for status check properties."""

    def test_is_pending(self):
        """is_pending returns True for PENDING status."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.PENDING),
        )

        assert mol.is_pending is True
        assert mol.is_complete is False
        assert mol.is_failed is False

    def test_is_complete(self):
        """is_complete returns True for OK status."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.OK),
        )

        assert mol.is_complete is True
        assert mol.is_pending is False

    def test_is_failed(self):
        """is_failed returns True for RERUN status."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.RERUN),
        )

        assert mol.is_failed is True

    def test_is_skipped(self):
        """is_skipped returns True for SKIP status."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.SKIP),
        )

        assert mol.is_skipped is True


class TestCanRerun:
    """Tests for can_rerun method."""

    def test_can_rerun_when_eligible(self):
        """can_rerun returns True when reruns < max_reruns."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.RERUN, reruns=1),
        )

        assert mol.can_rerun(max_reruns=3) is True

    def test_cannot_rerun_when_exhausted(self):
        """can_rerun returns False when reruns >= max_reruns."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.RERUN, reruns=3),
        )

        assert mol.can_rerun(max_reruns=3) is False

    def test_cannot_rerun_if_not_failed(self):
        """can_rerun returns False for non-RERUN status."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
            meta=MoleculeMeta(status=MoleculeStatus.PENDING, reruns=0),
        )

        assert mol.can_rerun(max_reruns=3) is False


class TestFromCsvDict:
    """Tests for from_csv_dict factory method."""

    def test_creates_molecule_from_dict(self):
        """from_csv_dict creates Molecule from dict."""
        row = {
            "mol_id": "mol_001",
            "smiles": "CCO",
            "nheavy": 3,
            "status": "Pending",
            "charge": 0,
            "multiplicity": 1,
        }

        mol = Molecule.from_csv_dict(row)

        assert mol.mol_id == "mol_001"
        assert mol.smiles == "CCO"
        assert mol.properties.nheavy == 3

    def test_missing_mol_id_raises(self):
        """from_csv_dict raises for missing mol_id."""
        row = {
            "smiles": "CCO",
            "nheavy": 3,
            "status": "Pending",
        }

        with pytest.raises(ValueError) as exc:
            Molecule.from_csv_dict(row)

        assert "mol_id" in str(exc.value)

    def test_handles_optional_fields(self):
        """from_csv_dict handles optional fields."""
        row = {
            "mol_id": "mol_001",
            "smiles": "CCO",
            "nheavy": 3,
            "status": "Pending",
            "H298_cbs": "",  # Empty optional field
            "delta_1": None,  # None optional field
        }

        mol = Molecule.from_csv_dict(row)

        assert mol.properties.H298_cbs is None
        assert mol.results.delta_1 is None


class TestToCsvDict:
    """Tests for to_csv_dict method."""

    def test_converts_to_dict(self):
        """to_csv_dict returns dict with all fields."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="CCO"),
            properties=MoleculeProperties(nheavy=3, H298_cbs=-100.5),
            meta=MoleculeMeta(status=MoleculeStatus.OK),
        )

        d = mol.to_csv_dict()

        assert d["mol_id"] == "mol_001"
        assert d["smiles"] == "CCO"
        assert d["nheavy"] == 3
        assert d["H298_cbs"] == -100.5
        assert d["status"] == "OK"

    def test_empty_values_as_empty_string(self):
        """to_csv_dict uses empty string for None values."""
        mol = Molecule(
            identity=MoleculeIdentity(mol_id="mol_001", smiles="C"),
            properties=MoleculeProperties(nheavy=1),
        )

        d = mol.to_csv_dict()

        assert d["H298_cbs"] == ""  # None -> ""
        assert d["delta_1"] == ""
