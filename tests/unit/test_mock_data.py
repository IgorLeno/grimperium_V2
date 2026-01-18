"""Tests for CLI mock data module."""

from grimperium.cli import mock_data


def test_cbs_original_not_in_databases() -> None:
    """
    Bug #1: CBS Original database file was deleted but reference remains.

    Verify that "CBS Reference (Original)" is NOT in the DATABASES list
    since the file data/phase_a/CBS_Original.csv was removed.
    """
    database_names = [db.name for db in mock_data.DATABASES]

    # This test should FAIL initially (CBS Original is still there)
    # After fix, it should PASS (CBS Original removed)
    assert "CBS Reference (Original)" not in database_names, (
        "CBS Reference (Original) should be removed - the CSV file no longer exists"
    )


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
    assert len(mock_data.DATABASES) == 3, (
        f"Expected 3 databases after CBS Original removal, got {len(mock_data.DATABASES)}"
    )
