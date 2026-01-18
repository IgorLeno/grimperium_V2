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

# Standard library imports
from collections import deque
from typing import Any

# Third-party imports
import numpy as np
from numpy.typing import NDArray

# Type aliases for mypy --strict compliance
MatrixFloat = NDArray[np.floating[Any]]
MatrixInt = NDArray[np.signedinteger[Any]]
Vector = NDArray[np.floating[Any]]
DictStrAny = dict[str, Any]
ListAny = list[Any]
DequeAny = deque[Any]

# Local module imports (must be after type aliases to avoid circular imports)
from grimperium.api import GrimperiumAPI  # noqa: E402
from grimperium.config import GrimperiumConfig  # noqa: E402

__all__ = [
    "GrimperiumAPI",
    "GrimperiumConfig",
    "__version__",
    # Type aliases
    "MatrixFloat",
    "MatrixInt",
    "Vector",
    "DictStrAny",
    "ListAny",
    "DequeAny",
]
