"""Enums for batch processing of CREST PM7 pipeline.

This module defines status enums, sorting strategies, and failure policies
used throughout the batch processing system.
"""

from enum import Enum


class MoleculeStatus(str, Enum):
    """Molecule processing status in batch pipeline.

    State transitions:
        PENDING -> SELECTED (via select_batch)
        SELECTED -> RUNNING (via mark_running)
        RUNNING -> OK (via mark_success)
        RUNNING -> RERUN (via mark_rerun, if retry_count < max_retries)
        RUNNING -> SKIP (via mark_skip, if retry_count >= max_retries)
        RERUN -> SELECTED (via select_batch, next batch)
        OK -> TERMINAL (no changes allowed)
        SKIP -> TERMINAL (no changes allowed)
        * -> PENDING (via reset_batch, only for ALL_OR_NOTHING policy)

    Values:
        - PENDING: "Pending" (Title Case)
        - SELECTED: "Selected" (Title Case)
        - RUNNING: "Running" (Title Case)
        - OK: "OK" (UPPERCASE - backward compatible with old CSV format)
        - RERUN: "Rerun" (Title Case)
        - SKIP: "Skip" (Title Case)

    BREAKING CHANGE (Round 2):
    The OK enum value was changed from "Ok" to "OK" for consistency
    and backward compatibility with older CSV formats.

    If you have existing CSVs with "Ok" status (title case), they will
    NOT match the new enum value "OK" (uppercase). You must normalize:

        df['status'] = df['status'].replace({'Ok': 'OK'})
        df.to_csv(csv_path, index=False)
    """

    PENDING = "Pending"
    """Initial state - molecule not yet processed."""

    SELECTED = "Selected"
    """Molecule selected for current batch, awaiting execution."""

    RUNNING = "Running"
    """Molecule currently being processed."""

    OK = "OK"
    """Processing completed successfully. Value: "OK" (uppercase)."""

    RERUN = "Rerun"
    """Processing failed, awaiting retry in next batch."""

    SKIP = "Skip"
    """Processing aborted - max retries exhausted."""


class BatchSortingStrategy(str, Enum):
    """Strategy for selecting and ordering molecules in a batch.

    Determines which molecules are picked first and in what order
    they appear in the batch.
    """

    RERUN_FIRST_THEN_EASY = "rerun_first_then_easy"
    """Prioritize RERUN molecules, then PENDING sorted by nheavy ascending."""

    RANDOM = "random"
    """Random selection from available molecules."""

    BY_NHEAVY = "by_nheavy"
    """Sort by number of heavy atoms (ascending - smallest first)."""

    BY_NROTBONDS = "by_nrotbonds"
    """Sort by number of rotatable bonds (ascending - least flexible first)."""


class BatchFailurePolicy(str, Enum):
    """Policy for handling batch failures.

    Determines what happens when some molecules in a batch fail.
    """

    PARTIAL_OK = "partial_ok"
    """Accept partial success - failed molecules marked RERUN/SKIP individually."""

    ALL_OR_NOTHING = "all_or_nothing"
    """If any molecule fails, reset entire batch to PENDING for retry."""
