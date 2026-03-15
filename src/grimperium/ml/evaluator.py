"""Model evaluation on full dataset."""

from __future__ import annotations

from pathlib import Path

from grimperium import DictStrAny, MatrixFloat
from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.data_loader import load_ml_data
from grimperium.ml.features import FeaturePipeline
from grimperium.ml.trainer import DEFAULT_FEATURE_COLS


def evaluate(
    learner: DeltaLearner,
    csv_path: Path,
    *,
    feature_cols: list[str] | None = None,
) -> DictStrAny:
    """Evaluate a trained DeltaLearner on the full dataset.

    Args:
        learner (DeltaLearner): A fitted DeltaLearner instance.
        csv_path (str): Path to thermo_pm7.csv.
        feature_cols (list[str], optional): Feature column names. If None, uses
            a defensive copy of DEFAULT_FEATURE_COLS.

    Returns:
        DictStrAny: Metric dictionary with keys `rmse`, `mae`, `r2`, `mape`,
        and `max_error` computed on the full dataset.

    Example:
        >>> from pathlib import Path
        >>> from grimperium.ml.trainer import train
        >>> learner, _, _ = train(Path("thermo_pm7.csv"))
        >>> metrics = evaluate(learner, Path("thermo_pm7.csv"))
        >>> sorted(metrics.keys())
        ['mae', 'mape', 'max_error', 'r2', 'rmse']
    """
    if feature_cols is None:
        feature_cols = list(DEFAULT_FEATURE_COLS)

    df, y_cbs, y_pm7 = load_ml_data(csv_path)

    pipeline = FeaturePipeline(feature_cols)
    X: MatrixFloat = pipeline.fit_transform(df)

    return learner.evaluate(X, y_cbs, y_pm7)
