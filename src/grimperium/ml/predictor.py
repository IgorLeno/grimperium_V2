"""Batch prediction for Phase F.

Applies a trained DeltaLearner model to all eligible molecules in a CSV,
writing ``H298_predicted`` and ``delta_correction`` columns back to disk.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal, overload

import numpy as np
import pandas as pd

from grimperium import DictStrAny

logger = logging.getLogger(__name__)


@overload
def predict_batch(
    csv_path: Path,
    model_path: Path,
    *,
    return_stats: Literal[False] = ...,
) -> pd.DataFrame:
    ...


@overload
def predict_batch(
    csv_path: Path,
    model_path: Path,
    *,
    return_stats: Literal[True],
) -> tuple[pd.DataFrame, DictStrAny]:
    ...


def predict_batch(
    csv_path: Path,
    model_path: Path,
    *,
    return_stats: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, DictStrAny]:
    """Apply a trained model to all eligible molecules in a CSV.

    Eligible rows have ``status == 'OK'`` AND ``cbs_quality_flag == 'OK'``.
    Non-eligible rows get NaN for both new columns.

    Args:
        csv_path: Path to thermo_pm7.csv (will be overwritten with new columns).
        model_path: Path to the trained ``.joblib`` model bundle.
        return_stats: If True, return ``(df, stats_dict)`` instead of just ``df``.

    Returns:
        DataFrame with ``H298_predicted`` and ``delta_correction`` columns added.
        If ``return_stats=True``, returns ``(df, stats)`` where stats has keys:
        ``n_predicted``, ``mean_delta``, ``std_delta``, ``min_predicted``,
        ``max_predicted``.

    Raises:
        FileNotFoundError: If model or CSV file does not exist.
    """
    from grimperium.ml.persistence import load_model

    # 1. Validate paths
    csv_path = Path(csv_path)
    model_path = Path(model_path)

    if not csv_path.exists():
        msg = f"CSV file not found: {csv_path}"
        raise FileNotFoundError(msg)

    # load_model raises FileNotFoundError("Model file not found: ...") itself

    # 2. Load model and data
    learner, pipeline = load_model(model_path)
    df = pd.read_csv(csv_path)

    # 3. Init new columns as NaN
    df["H298_predicted"] = np.nan
    df["delta_correction"] = np.nan

    # 4. Build eligible mask
    eligible_mask = (df["status"] == "OK") & (df["cbs_quality_flag"] == "OK")
    eligible_df = df.loc[eligible_mask].copy()

    if len(eligible_df) > 0:
        # 5. Transform features and extract y_pm7
        x_eligible = pipeline.transform(eligible_df)
        y_pm7 = eligible_df["H298_pm7"].to_numpy(dtype=np.float64)

        # 6. Predict
        h298_predicted = learner.predict(x_eligible, y_pm7)
        delta_correction = h298_predicted - y_pm7

        # 7. Assign back
        df.loc[eligible_mask, "H298_predicted"] = h298_predicted
        df.loc[eligible_mask, "delta_correction"] = delta_correction

    # 8. Overwrite CSV
    df.to_csv(csv_path, index=False)

    n_predicted = int(eligible_mask.sum())
    logger.info("Batch prediction complete: %d molecules predicted", n_predicted)

    if return_stats:
        predicted_vals = df.loc[eligible_mask, "H298_predicted"]
        delta_vals = df.loc[eligible_mask, "delta_correction"]
        stats: DictStrAny = {
            "n_predicted": n_predicted,
            "mean_delta": float(delta_vals.mean()) if n_predicted > 0 else 0.0,
            "std_delta": float(delta_vals.std(ddof=0)) if n_predicted > 0 else 0.0,
            "min_predicted": float(predicted_vals.min()) if n_predicted > 0 else 0.0,
            "max_predicted": float(predicted_vals.max()) if n_predicted > 0 else 0.0,
        }
        return df, stats

    return df
