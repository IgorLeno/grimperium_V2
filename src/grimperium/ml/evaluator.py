"""Model evaluation on full dataset."""

from __future__ import annotations

from pathlib import Path

from grimperium import DictStrAny, MatrixFloat
from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.data_loader import load_ml_data
from grimperium.ml.features import FeaturePipeline


def evaluate(
    learner: DeltaLearner,
    csv_path: Path,
    *,
    pipeline: FeaturePipeline | None = None,
) -> DictStrAny:
    """Evaluate a trained DeltaLearner on the full dataset.

    Args:
        learner (DeltaLearner): A fitted DeltaLearner instance.
        csv_path (Path): Path to thermo_pm7.csv.
        pipeline (FeaturePipeline | None): Pre-fitted feature pipeline obtained
            during training. Required to avoid fitting on evaluation data.

    Returns:
        DictStrAny: Metric dictionary with keys `rmse`, `mae`, `r2`, `mape`,
        and `max_error` computed on the full dataset.

    Example:
        >>> from pathlib import Path
        >>> from grimperium.ml.trainer import train
        >>> learner, _, _, pipeline = train(
        ...     Path("thermo_pm7.csv"),
        ...     return_pipeline=True,
        ... )
        >>> metrics = evaluate(
        ...     learner,
        ...     Path("thermo_pm7.csv"),
        ...     pipeline=pipeline,
        ... )
        >>> sorted(metrics.keys())
        ['mae', 'mape', 'max_error', 'r2', 'rmse']
    """
    if pipeline is None:
        msg = (
            "evaluate requires a pre-fitted FeaturePipeline to avoid "
            "evaluation leakage. Call train(..., return_pipeline=True) "
            "and pass pipeline=..."
        )
        raise ValueError(msg)

    df, y_cbs, y_pm7 = load_ml_data(csv_path)
    X: MatrixFloat = pipeline.transform(df)

    return learner.evaluate(X, y_cbs, y_pm7)
