"""Tests for CLI mock data module."""

from grimperium.cli import mock_data


def test_cbs_original_not_in_databases() -> None:
    """
    Verifies CBS Original database removed from mock databases.

    CBS Original file was deleted in Phase C BATCH 12 (Bug #1).
    This test confirms it no longer appears in the database menu.
    """
    database_names = [db.name for db in mock_data.DATABASES]

    assert (
        "CBS Reference (Original)" not in database_names
    ), "CBS Reference (Original) should be removed - the CSV file no longer exists"


def test_cbs_chon_still_exists() -> None:
    """Verify CBS CHON-only database still exists after CBS Original removal."""
    database_names = [db.name for db in mock_data.DATABASES]
    assert "CBS Reference (CHON-only)" in database_names


def test_crest_pm7_still_exists() -> None:
    """Verify CREST PM7 database still exists after CBS Original removal."""
    database_names = [db.name for db in mock_data.DATABASES]
    assert "CREST PM7" in database_names


def test_databases_count_after_cbs_original_removal() -> None:
    """After removing CBS Original, should have exactly 3 databases."""
    # CHON-only, CREST PM7, NIST Experimental = 3 total
    assert (
        len(mock_data.DATABASES) == 3
    ), f"Expected 3 databases after CBS Original removal, got {len(mock_data.DATABASES)}"
