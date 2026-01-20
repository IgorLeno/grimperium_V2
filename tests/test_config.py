# tests/test_config.py
import pytest
from src.grimperium.crest_pm7.config import MoleculeStatus

def test_molecule_status_enum_exists():
    """Test that MoleculeStatus enum has required values"""
    assert hasattr(MoleculeStatus, 'PENDING')
    assert hasattr(MoleculeStatus, 'CURRENT_BATCH')
    assert hasattr(MoleculeStatus, 'RUNNING')
    assert hasattr(MoleculeStatus, 'OK')
    assert hasattr(MoleculeStatus, 'RERUN')
    assert hasattr(MoleculeStatus, 'SKIP')

def test_molecule_status_values():
    """Test MoleculeStatus enum values"""
    assert MoleculeStatus.PENDING.value == "Pending"
    assert MoleculeStatus.CURRENT_BATCH.value == "Current Batch"  # Verify actual value
    assert MoleculeStatus.RUNNING.value == "Running"  # Verify actual value
    assert MoleculeStatus.OK.value == "OK"
    assert MoleculeStatus.RERUN.value == "Rerun"
    assert MoleculeStatus.SKIP.value == "Skip"
