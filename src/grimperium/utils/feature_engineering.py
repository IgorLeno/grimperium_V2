"""
Feature engineering for molecular data.

This module provides feature extraction from SMILES:
    - Morgan Fingerprints (circular fingerprints)
    - RDKit molecular descriptors
    - Tabular features from dataset

Feature vector composition:
    [tabular (3)] + [Morgan FP (256)] + [RDKit descriptors (10+)]
    Total: ~270 features

Example:
    >>> from grimperium.utils import FeatureEngineer
    >>> fe = FeatureEngineer(morgan_bits=256)
    >>> features = fe.transform(["CCO", "CC(=O)O"])
    >>> print(features.shape)  # (2, ~270)

"""

import pandas as pd

from grimperium import MatrixFloat


class FeatureEngineer:
    """
    Feature engineering pipeline for molecular data.

    Combines multiple feature sources:
        - Tabular: nheavy, charge, multiplicity
        - Morgan FP: Circular fingerprints from SMILES
        - RDKit: Molecular descriptors (MW, TPSA, LogP, etc.)

    Attributes:
        morgan_bits: Number of Morgan fingerprint bits
        morgan_radius: Morgan fingerprint radius
        rdkit_descriptors: List of RDKit descriptor names
        tabular_features: List of tabular feature names

    Example:
        >>> fe = FeatureEngineer(morgan_bits=512)
        >>> X = fe.fit_transform(smiles_list, df)
        >>> print(f"Features shape: {X.shape}")

    """

    # Default RDKit descriptors
    DEFAULT_RDKIT_DESCRIPTORS = [
        "MolWt",
        "TPSA",
        "MolLogP",
        "NumRotatableBonds",
        "NumHDonors",
        "NumHAcceptors",
        "NumHeteroatoms",
        "NumAromaticRings",
        "FractionCSP3",
        "HeavyAtomMolWt",
    ]

    def __init__(
        self,
        morgan_bits: int = 256,
        morgan_radius: int = 2,
        rdkit_descriptors: list[str] | None = None,
        tabular_features: list[str] | None = None,
        use_morgan: bool = True,
        use_rdkit: bool = True,
        use_tabular: bool = True,
    ) -> None:
        """
        Initialize FeatureEngineer.

        Args:
            morgan_bits: Number of bits for Morgan fingerprints
            morgan_radius: Radius for Morgan fingerprints
            rdkit_descriptors: RDKit descriptors to compute
            tabular_features: Tabular features to include
            use_morgan: Whether to include Morgan FP
            use_rdkit: Whether to include RDKit descriptors
            use_tabular: Whether to include tabular features

        """
        self.morgan_bits = morgan_bits
        self.morgan_radius = morgan_radius
        self.rdkit_descriptors = rdkit_descriptors or self.DEFAULT_RDKIT_DESCRIPTORS
        self.tabular_features = tabular_features or ["nheavy", "charge", "multiplicity"]
        self.use_morgan = use_morgan
        self.use_rdkit = use_rdkit
        self.use_tabular = use_tabular

        self._n_features: int | None = None
        self._feature_names: list[str] | None = None

    def fit(
        self,
        smiles: list[str],
        df: pd.DataFrame | None = None,
    ) -> "FeatureEngineer":
        """
        Fit the feature engineer (compute feature names).

        Args:
            smiles: List of SMILES strings
            df: Optional DataFrame with tabular features

        Returns:
            self for method chaining

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def transform(
        self,
        smiles: list[str],
        df: pd.DataFrame | None = None,
    ) -> MatrixFloat:
        """
        Transform SMILES to feature matrix.

        Args:
            smiles: List of SMILES strings
            df: Optional DataFrame with tabular features

        Returns:
            Feature matrix of shape (n_samples, n_features)

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def fit_transform(
        self,
        smiles: list[str],
        df: pd.DataFrame | None = None,
    ) -> MatrixFloat:
        """
        Fit and transform in one step.

        Args:
            smiles: List of SMILES strings
            df: Optional DataFrame with tabular features

        Returns:
            Feature matrix of shape (n_samples, n_features)

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def get_feature_names(self) -> list[str]:
        """Get list of feature names."""
        raise NotImplementedError("Will be implemented in Batch 3")

    @property
    def n_features(self) -> int:
        """Number of features."""
        if self._n_features is None:
            raise ValueError("Must call fit() first")
        return self._n_features

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"FeatureEngineer(morgan={self.morgan_bits}bits, "
            f"rdkit={len(self.rdkit_descriptors)}, "
            f"tabular={len(self.tabular_features)})"
        )


def compute_morgan_fingerprints(
    smiles: str | list[str],
    n_bits: int = 256,
    radius: int = 2,
) -> MatrixFloat:
    """
    Compute Morgan fingerprints from SMILES.

    Args:
        smiles: Single SMILES or list of SMILES
        n_bits: Number of fingerprint bits
        radius: Fingerprint radius

    Returns:
        Fingerprint array of shape (n_samples, n_bits)

    Example:
        >>> fps = compute_morgan_fingerprints(["CCO", "CC(=O)O"])
        >>> print(fps.shape)  # (2, 256)

    """
    raise NotImplementedError("Will be implemented in Batch 3")


def compute_rdkit_descriptors(
    smiles: str | list[str],
    descriptors: list[str] | None = None,
) -> MatrixFloat:
    """
    Compute RDKit molecular descriptors from SMILES.

    Args:
        smiles: Single SMILES or list of SMILES
        descriptors: List of descriptor names (None for defaults)

    Returns:
        Descriptor array of shape (n_samples, n_descriptors)

    Available descriptors:
        - MolWt: Molecular weight
        - TPSA: Topological polar surface area
        - MolLogP: Wildman-Crippen LogP
        - NumRotatableBonds: Number of rotatable bonds
        - NumHDonors: Number of H-bond donors
        - NumHAcceptors: Number of H-bond acceptors
        - And many more...

    Example:
        >>> descs = compute_rdkit_descriptors(["CCO"], ["MolWt", "TPSA"])
        >>> print(descs)  # [[46.07, 20.23]]

    """
    raise NotImplementedError("Will be implemented in Batch 3")


def extract_tabular_features(
    df: pd.DataFrame,
    features: list[str] | None = None,
) -> MatrixFloat:
    """
    Extract tabular features from DataFrame.

    Args:
        df: Source DataFrame
        features: List of column names to extract

    Returns:
        Feature array of shape (n_samples, n_features)

    """
    raise NotImplementedError("Will be implemented in Batch 3")
