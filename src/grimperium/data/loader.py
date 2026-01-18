"""Chemperium dataset loader.

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

import pandas as pd
from sklearn.model_selection import train_test_split

# Module-level constants for dataset paths
# PRIMARY dataset (CHON molecules only - delta-learning optimized)
THERMO_CBS_CHON_PATH = "data/thermo_cbs_chon.csv"  # Primary: 29,568 molecules (CHON only)

# SECONDARY dataset (PM7 optimization results)
THERMO_PM7_PATH = "data/thermo_pm7.csv"  # Secondary: PM7-optimized results from CREST pipeline

# DEPRECATED paths (kept for reference only - DO NOT USE)
# THERMO_CBS_OPT_PATH = "data/thermo_cbs_opt.csv"  # ❌ REMOVED (superseded by THERMO_CBS_CHON_PATH)
# THERMO_CBS_CLEAN_PATH = "data/thermo_cbs_clean.csv"  # ❌ REMOVED (superseded by THERMO_CBS_CHON_PATH)


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
        self.data: pd.DataFrame | None = None

    def load(
        self,
        path: str | Path,
        columns: list[str] | None = None,
        max_nheavy: int | None = None,
        allowed_elements: Iterable[str] | None = None,
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
        df: pd.DataFrame | None = None,
        test_size: float = 0.2,
        random_state: int = 42,
        stratify_by: str | None = None,
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
        df: pd.DataFrame | None = None,
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
        df: pd.DataFrame | None = None,
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
        df: pd.DataFrame | None = None,
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
        max_nheavy: int | None,
        allowed_elements: Iterable[str] | None,
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
    def load_thermo_pm7(
        cls,
        max_nheavy: int | None = None,
        validate: bool = True,
    ) -> pd.DataFrame:
        """Load PM7-optimized results dataset (secondary).
        
        File: thermo_pm7.csv
        Size: Results from CREST + PM7 optimization pipeline
        
        Properties:
        - PM7 enthalpy (semiempirical optimization)
        - Conformer details
        - Quality grades
        - Batch processing metadata
        
        Delta-learning context:
        Used as experimental validation target or alternative cheap baseline
        when validating delta-learning models trained on CBS/B3LYP delta.
        
        Dataset Characteristics:
        - Source: CREST conformer generation + PM7 optimization
        - Method: Semiempirical (PM7)
        - Use case: Experimental validation, semiempirical baseline
        
        Args:
            max_nheavy: Filter molecules with nheavy <= this value
            validate: Whether to validate data on load
        
        Returns:
            pd.DataFrame: PM7 properties and conformer data
        
        Raises:
            FileNotFoundError: If data/thermo_pm7.csv not found
            ValueError: If CSV format is invalid
        
        See Also:
            load_thermo_cbs_chon: Primary CHON dataset (recommended)
            docs/DATASETS.md: Comprehensive dataset documentation
            docs/CREST_INTEGRATION.md: CREST + PM7 pipeline details
        
        Example:
            >>> df = ChemperiumLoader.load_thermo_pm7()
            >>> assert 'pm7_enthalpy' in df.columns or 'H298_pm7' in df.columns
        """
        if not Path(THERMO_PM7_PATH).exists():
            raise FileNotFoundError(
                f"Secondary dataset not found: {THERMO_PM7_PATH}\n"
                f"This file contains PM7-optimized results from CREST pipeline\n"
                f"Note: thermo_batch_final.csv was renamed to thermo_pm7.csv"
            )
        
        loader = cls(validate=validate)
        df = loader.load(THERMO_PM7_PATH, max_nheavy=max_nheavy)
        return df

    @classmethod
    def load_thermo_cbs_chon(
        cls,
        max_nheavy: int | None = None,
        validate: bool = True,
    ) -> pd.DataFrame:
        """Load CHON thermochemistry dataset (primary).
        
        File: thermo_cbs_chon.csv
        Size: 29,568 molecules
        Composition: C, H, O, N only (no halogens, sulfur, or rare heteroatoms)
        
        Properties:
        - CBS-level enthalpy (high-accuracy reference)
        - B3LYP-level enthalpy (cheap alternative for delta-learning)
        - Molecular properties: mass, n_heavy_atoms, etc.
        
        Delta-learning context:
        This dataset is specifically curated for delta-learning models where
        the target is the correction (delta) between B3LYP (cheap) and CBS (accurate).
        The CHON-only constraint ensures:
        - Homogeneous electronic physics
        - No exotic valence states or relativistic effects
        - Learnable correction patterns
        
        Dataset Characteristics:
        - Total molecules: 29,568 (filtered from original 52,837)
        - Elements: C, H, O, N only
        - Removed: Halogens (F, Cl, Br, I), sulfur (S), rare heteroatoms (B, P, As, Ge)
        - Use case: Delta-learning model training and validation
        
        Chemical Implications:
        - Homogeneous chemical space (no highly polarizable atoms)
        - Well-conditioned delta-learning problem
        - Smaller energy delta variance vs. full dataset
        - No hypervalence extremes or relativistic effects
        
        Args:
            max_nheavy: Filter molecules with nheavy <= this value
            validate: Whether to validate data on load
        
        Returns:
            pd.DataFrame: Thermochemistry dataset with columns:
                          smiles, charge, multiplicity, nheavy, H298_cbs, H298_b3
                          Shape: (29568, M) where M = 5+ columns
        
        Raises:
            FileNotFoundError: If data/thermo_cbs_chon.csv not found
            ValueError: If CSV format is invalid
        
        See Also:
            load_thermo_pm7: PM7-optimized results (secondary dataset)
            docs/DATASETS.md: Comprehensive dataset documentation
        
        Example:
            >>> df = ChemperiumLoader.load_thermo_cbs_chon()
            >>> assert set(df.columns) >= {'smiles', 'charge', 'multiplicity', 'nheavy', 'H298_cbs'}
            >>> assert len(df) == 29568  # Full dataset (before max_nheavy filter)
        """
        if not Path(THERMO_CBS_CHON_PATH).exists():
            raise FileNotFoundError(
                f"Primary dataset not found: {THERMO_CBS_CHON_PATH}\n"
                f"Expected: 29,568 CHON molecules with CBS and B3LYP enthalpies\n"
                f"Note: thermo_cbs_clean.csv and thermo_cbs_opt.csv are no longer used"
            )
        
        loader = cls(validate=validate)
        df = loader.load(THERMO_CBS_CHON_PATH, max_nheavy=max_nheavy)
        return df
