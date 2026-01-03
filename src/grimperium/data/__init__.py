"""
Data loading and fusion module for Grimperium.

This module provides classes for:
    - Loading Chemperium dataset
    - Handling semiempirical (PM7) calculations
    - Fusing multiple data sources and computing deltas

Classes:
    ChemperiumLoader: Load and preprocess Chemperium dataset
    SemiempiricalHandler: Interface for PM7 calculations via MOPAC
    DataFusion: Merge datasets and compute delta values

"""

from grimperium.data.fusion import DataFusion
from grimperium.data.loader import ChemperiumLoader
from grimperium.data.semiempirical import SemiempiricalHandler

__all__ = [
    "ChemperiumLoader",
    "DataFusion",
    "SemiempiricalHandler",
]
