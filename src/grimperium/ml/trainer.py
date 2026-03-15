"""Training orchestration with data-leakage-free split."""

from __future__ import annotations

from pathlib import Path

from sklearn.model_selection import train_test_split

from grimperium import DictStrAny, MatrixFloat
from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.data_loader import load_ml_data
from grimperium.ml.features import FeaturePipeline

# Default feature columns matching thermo_pm7.csv
DEFAULT_FEATURE_COLS = [
    "nheavy",
    "rdkit_nrotbonds",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
]


def train(
    csv_path: Path,
    *,
    feature_cols: list[str] | None = None,
    test_size: float = 0.2,
    random_state: int = 42,
) -> tuple[DeltaLearner, DictStrAny, DictStrAny]:
    """Train a DeltaLearner with proper train/test split.

    CRITICAL: Split happens on the DataFrame BEFORE fit_transform
    to prevent data leakage through imputation statistics.

    Parameters:
        csv_path: Path to thermo_pm7.csv
        feature_cols: Feature column names (defaults to DEFAULT_FEATURE_COLS)
        test_size: Fraction of data for test set
        random_state: Random seed for reproducibility

    Returns:
        (learner, train_metrics, test_metrics)
    """
    if feature_cols is None:
        feature_cols = DEFAULT_FEATURE_COLS

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
    X_train: MatrixFloat = pipeline.fit_transform(df_train)
    X_test: MatrixFloat = pipeline.transform(df_test)

    # Step 4: Train DeltaLearner
    learner = DeltaLearner()
    learner.fit(X_train, y_cbs_train, y_pm7_train)

    # Step 5: Evaluate on both splits
    train_metrics: DictStrAny = learner.evaluate(X_train, y_cbs_train, y_pm7_train)
    test_metrics: DictStrAny = learner.evaluate(X_test, y_cbs_test, y_pm7_test)

    return learner, train_metrics, test_metrics
