"""
Grimperium - ML Ensemble Framework for Molecular Thermodynamic Property Prediction.

A production-ready framework implementing delta-learning for correction of
semiempirical methods (PM7) using ensemble ML models (KRR + XGBoost).

Features:
    - Delta-learning: Corrects PM7 predictions to CBS-level accuracy
    - Hybrid features: Tabular + Morgan Fingerprints + RDKit descriptors
    - Ensemble models: Kernel Ridge Regression + XGBoost
    - Production-ready: Full CI/CD, type hints, comprehensive testing

Example:
    >>> from grimperium import GrimperiumAPI
    >>> api = GrimperiumAPI()
    >>> predictions = api.predict(smiles_list)

"""

__version__ = "0.2.0"
__author__ = "Grimperium Team"
__license__ = "MIT"

from grimperium.api import GrimperiumAPI
from grimperium.config import GrimperiumConfig

__all__ = [
    "GrimperiumAPI",
    "GrimperiumConfig",
    "__version__",
]
