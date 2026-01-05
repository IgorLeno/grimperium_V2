"""Test fixtures for Grimperium."""

from tests.fixtures.real_data import (
    get_dataset_stats,
    load_real_subset,
    load_real_train_test_split,
)

__all__ = [
    "load_real_subset",
    "load_real_train_test_split",
    "get_dataset_stats",
]