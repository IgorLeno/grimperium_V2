"""Training orchestration with data-leakage-free split."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, overload

from sklearn.model_selection import train_test_split

from grimperium import DictStrAny, MatrixFloat
from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.data_loader import load_ml_data
from grimperium.ml.features import FeaturePipeline

# Default feature columns matching thermo_pm7.csv
DEFAULT_FEATURE_COLS = (
    "nheavy",
    "rdkit_nrotbonds",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
)


@overload
def train(
    csv_path: Path,
    *,
    feature_cols: list[str] | None = None,
    test_size: float = 0.2,
    random_state: int = 42,
    return_pipeline: Literal[False] = False,
) -> tuple[DeltaLearner, DictStrAny, DictStrAny]: ...


@overload
def train(
    csv_path: Path,
    *,
    feature_cols: list[str] | None = None,
    test_size: float = 0.2,
    random_state: int = 42,
    return_pipeline: Literal[True],
) -> tuple[DeltaLearner, DictStrAny, DictStrAny, FeaturePipeline]: ...


def train(
    csv_path: Path,
    *,
    feature_cols: list[str] | None = None,
    test_size: float = 0.2,
    random_state: int = 42,
    return_pipeline: bool = False,
) -> (
    tuple[DeltaLearner, DictStrAny, DictStrAny]
    | tuple[DeltaLearner, DictStrAny, DictStrAny, FeaturePipeline]
):
    """Train a DeltaLearner with proper train/test split.

    CRITICAL: Split happens on the DataFrame BEFORE fit_transform
    to prevent data leakage through imputation statistics.

    Args:
        csv_path: Path to `thermo_pm7.csv`.
        feature_cols: Feature column names. If None, uses
            `DEFAULT_FEATURE_COLS`.
        test_size: Fraction of data for the test split. Must satisfy
            0 < test_size < 1.
        random_state: Random seed for reproducible splitting.
        return_pipeline: If True, also returns the fitted FeaturePipeline for
            leakage-safe downstream evaluation.

    Returns:
        Tuple with `(learner, train_metrics, test_metrics)` and, when
        `return_pipeline=True`, an additional fitted `FeaturePipeline`.

    Example:
        >>> from pathlib import Path
        >>> learner, train_metrics, test_metrics = train(
        ...     Path("thermo_pm7.csv"),
        ...     test_size=0.2,
        ...     random_state=42,
        ... )
        >>> train_metrics["rmse"] >= 0
        True
    """
    if not 0 < test_size < 1:
        msg = f"test_size must be between 0 and 1 (exclusive), got {test_size}"
        raise ValueError(msg)

    if feature_cols is None:
        feature_cols = list(DEFAULT_FEATURE_COLS)

    # Step 1: Load data
    df, y_cbs, y_pm7 = load_ml_data(csv_path)

    # Step 2: Split DataFrame BEFORE feature engineering (no leakage!)
    (
        df_train,
        df_test,
        y_cbs_train,
        y_cbs_test,
        y_pm7_train,
        y_pm7_test,
    ) = train_test_split(
        df,
        y_cbs,
        y_pm7,
        test_size=test_size,
        random_state=random_state,
    )

    # Step 3: Feature engineering — imputer fitted on train only
    pipeline = FeaturePipeline(feature_cols)
    X_train: MatrixFloat = pipeline.train(df_train)
    X_test: MatrixFloat = pipeline.transform(df_test)

    # Step 4: Train DeltaLearner
    learner = DeltaLearner()
    learner.fit(X_train, y_cbs_train, y_pm7_train)

    # Step 5: Evaluate on both splits
    train_metrics: DictStrAny = learner.evaluate(
        X_train, y_cbs_train, y_pm7_train
    )
    test_metrics: DictStrAny = learner.evaluate(X_test, y_cbs_test, y_pm7_test)

    if return_pipeline:
        return learner, train_metrics, test_metrics, pipeline
    return learner, train_metrics, test_metrics
