"""
Chemperium dataset loader.

This module provides the ChemperiumLoader class for loading and
preprocessing the Chemperium dataset containing ~52k molecules
with CBS-level thermodynamic properties.

Dataset columns:
    - smiles: SMILES string representation
    - xyz: Cartesian coordinates (optimized geometry)
    - charge: Total molecular charge
    - multiplicity: Spin multiplicity
    - nheavy: Number of heavy atoms
    - H298_cbs: Enthalpy at 298K (CBS reference)
    - H298_b3: Enthalpy at 298K (B3LYP)
    - S298: Entropy at 298K
    - cp_1 to cp_45: Heat capacity at different temperatures

Example:
    >>> from grimperium.data import ChemperiumLoader
    >>> loader = ChemperiumLoader()
    >>> df = loader.load("chemperium.csv")
    >>> train, test = loader.split(df, test_size=0.2)

"""

import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from sklearn.model_selection import train_test_split

# Module-level constants for dataset paths
THERMO_CBS_OPT_PATH = "data/thermo_cbs_opt.csv"  # Original: 52,837 molecules
THERMO_CBS_CLEAN_PATH = "data/thermo_cbs_clean.csv"  # Filtered: 30,026 molecules (Phase A)


class ChemperiumLoader:
    """
    Loader for Chemperium thermodynamic dataset.

    Handles loading, validation, and preprocessing of the Chemperium
    dataset containing CBS-level thermodynamic properties for ~52k molecules.

    Attributes:
        data: Loaded DataFrame (None until load() is called)
        validate: Whether to validate data on load

    Example:
        >>> loader = ChemperiumLoader()
        >>> df = loader.load("path/to/chemperium.csv")
        >>> print(f"Loaded {len(df)} molecules")
        >>> train, test = loader.split(df, test_size=0.2)

    """

    # Expected columns in Chemperium dataset
    REQUIRED_COLUMNS = [
        "smiles",
        "charge",
        "multiplicity",
        "nheavy",
        "H298_cbs",
    ]

    OPTIONAL_COLUMNS = [
        "xyz",
        "H298_b3",
        "S298",
        "A",
        "B",
    ]

    # Heat capacity columns (cp_1 to cp_45)
    CP_COLUMNS = [f"cp_{i}" for i in range(1, 46)]

    # Feature columns for ML (excluding targets)
    FEATURE_COLUMNS = ["nheavy", "charge", "multiplicity"]

    def __init__(self, validate: bool = True) -> None:
        """
        Initialize ChemperiumLoader.

        Args:
            validate: Whether to validate data on load

        """
        self.validate = validate
        self.data: Optional[pd.DataFrame] = None

    def load(
        self,
        path: Union[str, Path],
        columns: Optional[list[str]] = None,
        max_nheavy: Optional[int] = None,
        allowed_elements: Optional[Iterable[str]] = None,
    ) -> pd.DataFrame:
        """
        Load Chemperium dataset from file.

        Args:
            path: Path to CSV or Parquet file
            columns: Specific columns to load (None for all)
            max_nheavy: Filter molecules with nheavy <= this value
            allowed_elements: Filter by allowed elements (stub, requires RDKit)

        Returns:
            Loaded DataFrame

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If required columns are missing

        """
        path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")

        # Read based on file extension
        if path.suffix == ".csv":
            df = pd.read_csv(path, usecols=columns)
        elif path.suffix in {".parquet", ".pq"}:
            df = pd.read_parquet(path, columns=columns)
        else:
            raise ValueError(f"Unsupported file type: {path.suffix}")

        # Validate if enabled
        if self.validate:
            self._validate_dataframe(df)

        # Apply filters
        df = self._apply_filters(df, max_nheavy, allowed_elements)

        # Store reference
        self.data = df

        return df

    def split(
        self,
        df: Optional[pd.DataFrame] = None,
        test_size: float = 0.2,
        random_state: int = 42,
        stratify_by: Optional[str] = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Split data into train and test sets.

        Args:
            df: DataFrame to split (uses self.data if None)
            test_size: Fraction for test set
            random_state: Random seed for reproducibility
            stratify_by: Column for stratified split

        Returns:
            Tuple of (train_df, test_df)

        Raises:
            ValueError: If no data available

        """
        if df is None:
            df = self.data
        if df is None:
            raise ValueError("No data available. Call load() first or provide df.")

        stratify = df[stratify_by] if stratify_by else None

        train_df, test_df = train_test_split(
            df,
            test_size=test_size,
            random_state=random_state,
            stratify=stratify,
            shuffle=True,
        )

        return train_df, test_df

    def train_val_test_split(
        self,
        df: Optional[pd.DataFrame] = None,
        test_size: float = 0.2,
        val_size: float = 0.1,
        random_state: int = 42,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Split data into train, validation, and test sets.

        Args:
            df: DataFrame to split (uses self.data if None)
            test_size: Fraction for test set
            val_size: Fraction for validation set
            random_state: Random seed for reproducibility

        Returns:
            Tuple of (train_df, val_df, test_df)

        Raises:
            ValueError: If no data available or invalid split sizes

        """
        if df is None:
            df = self.data
        if df is None:
            raise ValueError("No data available. Call load() first or provide df.")

        if test_size + val_size >= 1.0:
            raise ValueError("test_size + val_size must be less than 1.0")

        # First split: separate test set
        train_val_df, test_df = train_test_split(
            df,
            test_size=test_size,
            random_state=random_state,
            shuffle=True,
        )

        # Second split: separate validation from train
        # Adjust val_size relative to remaining data
        relative_val_size = val_size / (1.0 - test_size)

        train_df, val_df = train_test_split(
            train_val_df,
            test_size=relative_val_size,
            random_state=random_state,
            shuffle=True,
        )

        return train_df, val_df, test_df

    def get_features(
        self,
        df: Optional[pd.DataFrame] = None,
        include_cp: bool = False,
    ) -> pd.DataFrame:
        """
        Extract feature columns from dataset.

        Args:
            df: Source DataFrame (uses self.data if None)
            include_cp: Whether to include heat capacity columns

        Returns:
            DataFrame with feature columns only

        Raises:
            ValueError: If no data available

        """
        if df is None:
            df = self.data
        if df is None:
            raise ValueError("No data available. Call load() first or provide df.")

        # Collect feature columns that exist in df
        feature_cols = [c for c in self.FEATURE_COLUMNS if c in df.columns]

        if include_cp:
            cp_cols = [c for c in self.CP_COLUMNS if c in df.columns]
            feature_cols.extend(cp_cols)

        return df[feature_cols].copy()

    def get_targets(
        self,
        df: Optional[pd.DataFrame] = None,
        target: str = "H298_cbs",
    ) -> pd.Series:
        """
        Extract target column from dataset.

        Args:
            df: Source DataFrame (uses self.data if None)
            target: Target column name

        Returns:
            Series with target values

        Raises:
            ValueError: If no data available or target column missing

        """
        if df is None:
            df = self.data
        if df is None:
            raise ValueError("No data available. Call load() first or provide df.")

        if target not in df.columns:
            raise ValueError(f"Target column '{target}' not found in DataFrame")

        return df[target].copy()

    def _validate_dataframe(self, df: pd.DataFrame) -> None:
        """
        Validate DataFrame has required columns.

        Args:
            df: DataFrame to validate

        Raises:
            ValueError: If required columns are missing

        """
        missing = [c for c in self.REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

    def _apply_filters(
        self,
        df: pd.DataFrame,
        max_nheavy: Optional[int],
        allowed_elements: Optional[Iterable[str]],
    ) -> pd.DataFrame:
        """
        Apply optional filters to DataFrame.

        Args:
            df: DataFrame to filter
            max_nheavy: Maximum number of heavy atoms
            allowed_elements: Allowed chemical elements (stub)

        Returns:
            Filtered DataFrame

        """
        if max_nheavy is not None and "nheavy" in df.columns:
            df = df[df["nheavy"] <= max_nheavy]

        if allowed_elements is not None:
            # Stub: Element filtering requires RDKit for SMILES parsing
            # Will be implemented in a future batch with RDKit dependency
            pass

        return df

    def __repr__(self) -> str:
        """String representation."""
        n_rows = len(self.data) if self.data is not None else 0
        return f"ChemperiumLoader(n_molecules={n_rows})"

    @classmethod
    def load_thermo_cbs_opt(
        cls,
        path: Union[str, Path] = THERMO_CBS_CLEAN_PATH,
        max_nheavy: Optional[int] = None,
        validate: bool = True,
    ) -> pd.DataFrame:
        """Load thermodynamic CBS dataset (DEPRECATED - use load_thermo_cbs_clean).

        This method loads the original unfiltered CBS dataset containing 52,837 molecules.

        DEPRECATION NOTE:
        - Default path changed 2026-01-10: now points to cleaned dataset by default
        - Use load_thermo_cbs_clean() for new code (recommended)
        - For original unfiltered data, pass file_path="data/thermo_cbs_opt.csv"
        - Deprecation timeline: 2026-01-10 (warning) → 2026-06-10 (official) → 2026-12-10 (removal TBD)

        Args:
            path: Path to CBS dataset CSV file.
                  Default: THERMO_CBS_CLEAN_PATH (cleaned, 30,026 molecules)
                  For original: "data/thermo_cbs_opt.csv" (unfiltered, 52,837 molecules)
            max_nheavy: Filter molecules with nheavy <= this value
            validate: Whether to validate data on load

        Returns:
            pd.DataFrame: Thermodynamic dataset with columns:
                          smiles, charge, multiplicity, nheavy, H298_cbs
                          Shape: (N, M) where N = molecules, M = 5 + optional columns
                          See REQUIRED_COLUMNS class constant for exact column definitions.

        Raises:
            FileNotFoundError: If path does not exist
            ValueError: If CSV format is invalid

        Warnings:
            DeprecationWarning: When using default path (now points to clean dataset)

        See Also:
            load_thermo_cbs_clean: Recommended method for Phase A onwards
            docs/DATASET_MIGRATION.md: Complete migration guide and dataset comparison
        """
        # Issue deprecation warning if using new default
        if path == THERMO_CBS_CLEAN_PATH:
            warnings.warn(
                "load_thermo_cbs_opt() default behavior changed (2026-01-10). "
                "Old default was 'data/thermo_cbs_opt.csv' (original dataset). "
                "New default is 'data/thermo_cbs_clean.csv' (filtered dataset). "
                "For new code, use load_thermo_cbs_clean() explicitly. "
                "To load original dataset, pass path='data/thermo_cbs_opt.csv'.",
                DeprecationWarning,
                stacklevel=2,
            )

        loader = cls(validate=validate)
        df = loader.load(path, max_nheavy=max_nheavy)
        return df

    @classmethod
    def load_thermo_cbs_clean(
        cls,
        max_nheavy: Optional[int] = None,
        validate: bool = True,
    ) -> pd.DataFrame:
        """Load cleaned thermodynamic CBS dataset (PRIMARY - Phase A onwards).

        Loads filtered dataset with 30,026 molecules (56.8% of original).
        Removes halogenated and sulfur-containing molecules for Phase A scope.

        Dataset Characteristics:
        - Total molecules: 30,026 (filtered from original 52,837)
        - Elements: C, H, O, N only (no halogens, no sulfur)
        - Filtering: 22,811 molecules removed (43.2% of original)
        - Use case: Phase A validation and forward development

        Filtering Applied:
        1. Halogenated molecules removed (F, Cl, Br, I)
        2. Sulfur-containing molecules removed (S)
        3. Non-essential columns removed

        Chemical Implications:
        - Homogeneous chemical space (no highly polarizable atoms)
        - Well-conditioned delta-learning problem
        - Smaller energy delta variance vs. full dataset

        Args:
            max_nheavy: Filter molecules with nheavy <= this value
            validate: Whether to validate data on load

        Returns:
            pd.DataFrame: Cleaned thermodynamic dataset with columns:
                          smiles, charge, multiplicity, nheavy, H298_cbs
                          Shape: (30026, M) where M = 5 + optional columns

        Raises:
            FileNotFoundError: If data/thermo_cbs_clean.csv does not exist
            ValueError: If CSV format is invalid

        See Also:
            load_thermo_cbs_opt: Legacy method (deprecated)
            docs/DATASET_MIGRATION.md: Dataset filtering details and statistics

        Example:
            >>> loader = ChemperiumLoader()
            >>> df = loader.load_thermo_cbs_clean()
            >>> assert len(df) == 30026
            >>> assert set(df.columns) >= {'smiles', 'charge', 'multiplicity', 'nheavy', 'H298_cbs'}
        """
        loader = cls(validate=validate)
        df = loader.load(THERMO_CBS_CLEAN_PATH, max_nheavy=max_nheavy)
        return df
