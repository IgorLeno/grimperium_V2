"""
Data fusion module for combining Chemperium and semiempirical data.

This module provides the DataFusion class for:
    - Merging Chemperium CBS data with PM7 calculations
    - Computing delta values (H298_CBS - H298_PM7)
    - Creating task-specific data views
    - Aligning and validating merged datasets

The delta-learning approach:
    delta = H298_CBS - H298_PM7
    prediction = H298_PM7 + model.predict(features)

Example:
    >>> from grimperium.data import DataFusion
    >>> fusion = DataFusion()
    >>> merged = fusion.merge(chemperium_df, pm7_df)
    >>> deltas = fusion.compute_deltas(merged)

"""

import warnings
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd

# Supported task types for task views
TaskName = Literal["enthalpy", "entropy", "heat_capacity"]


class DataFusion:
    """
    Fuse Chemperium and semiempirical datasets.

    Handles merging of high-accuracy CBS data with fast semiempirical
    calculations (PM7), computes delta values for ML training, and
    creates task-specific data views.

    Attributes:
        merged_data: Merged DataFrame (None until merge() is called)
        delta_column: Name of computed delta column

    Example:
        >>> fusion = DataFusion()
        >>> merged = fusion.merge(chemperium_df, pm7_df, on="smiles")
        >>> print(f"Delta mean: {merged['delta_pm7'].mean():.2f} kcal/mol")
        >>>
        >>> # Or use task views directly
        >>> X, y = fusion.select_task_view(df, task="enthalpy")

    """

    # Default feature columns for ML
    FEATURE_COLUMNS = ["nheavy", "charge", "multiplicity"]

    def __init__(
        self,
        target_column: str = "H298_cbs",
        semiempirical_column: str = "H298_pm7",
        delta_column: str = "delta_pm7",
    ) -> None:
        """
        Initialize DataFusion.

        Args:
            target_column: CBS target column name
            semiempirical_column: Semiempirical column name
            delta_column: Output delta column name

        """
        self.target_column = target_column
        self.semiempirical_column = semiempirical_column
        self.delta_column = delta_column
        self.merged_data: Optional[pd.DataFrame] = None

    def merge(
        self,
        chemperium_df: pd.DataFrame,
        semiempirical_df: pd.DataFrame,
        on: str = "smiles",
        how: str = "inner",
        validate_merge: bool = True,
    ) -> pd.DataFrame:
        """
        Merge Chemperium and semiempirical DataFrames.

        Args:
            chemperium_df: DataFrame with CBS data
            semiempirical_df: DataFrame with PM7 data
            on: Column to merge on
            how: Merge type ('inner', 'left', 'outer')
            validate_merge: Whether to validate merge result

        Returns:
            Merged DataFrame

        Raises:
            ValueError: If merge validation fails

        """
        # Perform merge
        merged = chemperium_df.merge(semiempirical_df, on=on, how=how)

        # Validate if requested
        if validate_merge:
            expected_rows = len(chemperium_df) if how == "inner" else None
            self._validate_merge(merged, expected_rows)

        # Store reference
        self.merged_data = merged

        return merged

    def compute_deltas(
        self,
        df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Compute delta values: delta = CBS - semiempirical.

        Args:
            df: DataFrame with both CBS and semiempirical columns

        Returns:
            DataFrame with added delta column

        Raises:
            ValueError: If required columns are missing

        """
        if df is None:
            df = self.merged_data
        if df is None:
            raise ValueError("No data available. Call merge() first or provide df.")

        # Check required columns exist
        if self.target_column not in df.columns:
            raise ValueError(f"Target column '{self.target_column}' not found")
        if self.semiempirical_column not in df.columns:
            raise ValueError(
                f"Semiempirical column '{self.semiempirical_column}' not found"
            )

        # Compute delta
        df = df.copy()
        df[self.delta_column] = df[self.target_column] - df[self.semiempirical_column]

        # Update stored data
        self.merged_data = df

        return df

    def select_task_view(
        self,
        df: pd.DataFrame,
        task: TaskName = "enthalpy",
    ) -> tuple[pd.DataFrame, Union[pd.Series, pd.DataFrame]]:
        """
        Create a task-specific view of the data.

        Args:
            df: Source DataFrame
            task: Task type - "enthalpy", "entropy", or "heat_capacity"

        Returns:
            Tuple of (X, y) where:
                - X: Feature DataFrame
                - y: Target Series (enthalpy/entropy) or DataFrame (heat_capacity)

        Raises:
            ValueError: If task is invalid or required columns missing

        """
        if task == "enthalpy":
            target_col = "H298_cbs"
            if target_col not in df.columns:
                raise ValueError(f"Column '{target_col}' required for enthalpy task")
            X = df[self._default_feature_columns(df, exclude=[target_col])].copy()
            y = df[target_col].copy()
            return X, y

        if task == "entropy":
            target_col = "S298"
            if target_col not in df.columns:
                raise ValueError(f"Column '{target_col}' required for entropy task")
            X = df[self._default_feature_columns(df, exclude=[target_col])].copy()
            y = df[target_col].copy()
            return X, y

        if task == "heat_capacity":
            cp_columns = [c for c in df.columns if c.startswith("cp_")]
            if not cp_columns:
                raise ValueError("No cp_* columns found for heat_capacity task")
            X = df[self._default_feature_columns(df, exclude=cp_columns)].copy()
            Y = df[cp_columns].copy()
            return X, Y

        raise ValueError(f"Invalid task: {task}. Must be one of: enthalpy, entropy, heat_capacity")

    def get_training_data(
        self,
        df: Optional[pd.DataFrame] = None,
    ) -> tuple[pd.DataFrame, np.ndarray]:
        """
        Get feature DataFrame and delta targets for training.

        Args:
            df: Source DataFrame (uses merged_data if None)

        Returns:
            Tuple of (features_df, delta_array)

        Raises:
            ValueError: If no data or delta column missing

        """
        if df is None:
            df = self.merged_data
        if df is None:
            raise ValueError("No data available. Call merge() first or provide df.")

        if self.delta_column not in df.columns:
            raise ValueError(
                f"Delta column '{self.delta_column}' not found. "
                "Call compute_deltas() first."
            )

        # Get feature columns
        feature_cols = self._default_feature_columns(df, exclude=[self.delta_column])
        features = df[feature_cols].copy()
        deltas = df[self.delta_column].values

        return features, deltas

    def analyze_deltas(
        self,
        df: Optional[pd.DataFrame] = None,
    ) -> dict[str, float]:
        """
        Compute statistics on delta distribution.

        Args:
            df: DataFrame with delta column

        Returns:
            Dict with mean, std, min, max, median

        Raises:
            ValueError: If no data or delta column missing

        """
        if df is None:
            df = self.merged_data
        if df is None:
            raise ValueError("No data available. Call merge() first or provide df.")

        if self.delta_column not in df.columns:
            raise ValueError(
                f"Delta column '{self.delta_column}' not found. "
                "Call compute_deltas() first."
            )

        deltas = df[self.delta_column]

        return {
            "mean": float(deltas.mean()),
            "std": float(deltas.std()),
            "min": float(deltas.min()),
            "max": float(deltas.max()),
            "median": float(deltas.median()),
        }

    def _validate_merge(
        self,
        merged: pd.DataFrame,
        expected_rows: Optional[int],
    ) -> None:
        """
        Validate merge result.

        Args:
            merged: Merged DataFrame
            expected_rows: Expected number of rows (None to skip check)

        Raises:
            ValueError: If merge produced empty result

        """
        if len(merged) == 0:
            raise ValueError("Merge produced empty result. Check merge columns.")

        # Log warning if significant row loss (for inner joins)
        if expected_rows is not None and len(merged) < expected_rows * 0.5:
            warnings.warn(
                f"Merge retained only {len(merged)}/{expected_rows} rows. "
                "Check for mismatched merge keys.",
                stacklevel=2,
            )

    def _default_feature_columns(
        self,
        df: pd.DataFrame,
        exclude: Optional[list[str]] = None,
    ) -> list[str]:
        """
        Get default feature columns from DataFrame.

        Args:
            df: Source DataFrame
            exclude: Columns to exclude from features

        Returns:
            List of feature column names

        """
        exclude_set = set(exclude or [])
        return [c for c in self.FEATURE_COLUMNS if c in df.columns and c not in exclude_set]

    def __repr__(self) -> str:
        """String representation."""
        n_rows = len(self.merged_data) if self.merged_data is not None else 0
        return f"DataFusion(n_merged={n_rows}, delta='{self.delta_column}')"
