"""Test fixtures for Grimperium."""

from tests.fixtures.real_data import (
    HAS_PM7_DATASET,
    HAS_REAL_DATASET,
    HAS_REQUIRED_DATASETS,
    get_dataset_stats,
    load_real_subset,
    load_real_train_test_split,
)

__all__ = [
    "HAS_PM7_DATASET",
    "HAS_REAL_DATASET",
    "HAS_REQUIRED_DATASETS",
    "load_real_subset",
    "load_real_train_test_split",
    "get_dataset_stats",
]
