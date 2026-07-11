"""ViewModels for GRIMPERIUM CLI (UI-facing data structures)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class PredictionResult:
    """Represents a prediction result from the real pipeline."""

    smiles: str
    h298_pm7: float  # Raw MOPAC HOF, kcal/mol
    delta_correction: float  # Model-predicted delta, kcal/mol
    h298_corrected: float  # h298_pm7 + delta_correction
    model_name: str
    model_version: str  # From metadata
    execution_time: float  # Total seconds
    n_conformers: int  # Conformers from CREST
