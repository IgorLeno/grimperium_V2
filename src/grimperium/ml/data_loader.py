"""ML-specific data loading: filter CSV for training-ready rows."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from grimperium import MatrixFloat
from grimperium.core.csv_data_loader import CSVDataLoader

logger = logging.getLogger(__name__)


def load_ml_data(
    csv_path: Path,
    *,
    strict: bool = False,
) -> tuple[pd.DataFrame, MatrixFloat, MatrixFloat]:
    """Load CSV and return only training-ready rows.

    Filters:
        1. status == 'OK' only
        2. Drops rows where H298_cbs or H298_pm7 is NaN

    Args:
        csv_path: Path to thermo_pm7.csv
        strict: If True, raise on validation errors; False = skip invalid rows

    Returns:
        (df, y_cbs, y_pm7) where df has feature columns and targets extracted

    Raises:
        ValueError: If no valid rows remain after filtering
    """
    loader = CSVDataLoader(csv_path, strict=strict)
    df = loader.load_dataframe()

    # Filter: status == OK
    df = df[df["status"] == "OK"].copy()

    # Filter: both targets present
    df = df.dropna(subset=["H298_cbs", "H298_pm7"])

    if len(df) == 0:
        msg = "No valid rows after filtering (status=OK, non-null targets)"
        raise ValueError(msg)

    logger.info("Loaded %d training-ready rows from %s", len(df), csv_path)

    y_cbs: MatrixFloat = df["H298_cbs"].to_numpy(dtype=np.float64)
    y_pm7: MatrixFloat = df["H298_pm7"].to_numpy(dtype=np.float64)

    return df, y_cbs, y_pm7
