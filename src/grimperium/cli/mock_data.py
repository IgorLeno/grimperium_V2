"""
Mock data for GRIMPERIUM CLI MVP.

This module provides mock data for the CLI interface before
real integration with the ML models and databases.
"""

from dataclasses import dataclass
from datetime import date, datetime
from random import Random


@dataclass
class Database:
    """Represents a molecular database."""

    name: str
    description: str
    molecules: int
    last_updated: datetime | date | None
    status: str  # "ready" | "in_development"
    properties: list[str]


@dataclass
class Model:
    """Represents a trained ML model."""

    name: str
    algorithm: str
    mae: float | None  # Mean Absolute Error
    r2: float | None  # R² score
    training_date: date | None
    status: str  # "ready" | "in_development"
    hyperparameters: dict[str, str | int | float]
    file_size: str | None  # e.g., "12.5 MB"


@dataclass
class PredictionResult:
    """Represents a prediction result."""

    smiles: str
    molecule_name: str
    property_name: str
    predicted_value: float
    unit: str
    confidence: float
    model_used: str


@dataclass
class DivergenceStats:
    """Statistics for CBS vs PM7 divergence."""

    severity: str  # LOW, MEDIUM, HIGH, CRITICAL
    range_min: float
    range_max: float
    count: int
    percentage: float


# Mock databases
# NOTE: CREST PM7 is now loaded dynamically from phase_a_results.json
# by DatabasesView.load_real_phase_a_results(). This entry serves as
# fallback when the file doesn't exist (dev/demo mode).
DATABASES: list[Database] = [
    Database(
        name="CBS Reference",
        description=(
            "Complete Basis Set (CBS) reference energies from "
            "Chemperium"
        ),
        molecules=30026,
        last_updated=date(2026, 1, 11),
        status="ready",
        properties=["H298_cbs", "H298_b3", "smiles", "charge", "multiplicity"],
    ),
    Database(
        name="CREST PM7",
        description="CREST conformer search with PM7 optimization",
        molecules=0,  # Real count comes from phase_a_results.json
        last_updated=date(2026, 1, 1),
        status="in_development",  # Reflects that no calculations done yet
        properties=["H298_pm7", "conformers", "smiles", "quality_grade"],
    ),
    Database(
        name="NIST Experimental",
        description="NIST experimental thermochemistry data",
        molecules=0,
        last_updated=date(2026, 1, 1),
        status="in_development",
        properties=["H298_exp", "uncertainty", "reference"],
    ),
]

# Mock models
MODELS: list[Model] = [
    Model(
        name="DeltaRF_v1.0",
        algorithm="Random Forest",
        mae=2.34,
        r2=0.945,
        training_date=date(2026, 1, 12),
        status="ready",
        hyperparameters={
            "n_estimators": 500,
            "max_depth": 20,
            "min_samples_split": 5,
            "random_state": 42,
        },
        file_size="45.2 MB",
    ),
    Model(
        name="DeltaXGB_v1.0",
        algorithm="XGBoost",
        mae=2.12,
        r2=0.952,
        training_date=date(2026, 1, 12),
        status="ready",
        hyperparameters={
            "n_estimators": 300,
            "max_depth": 8,
            "learning_rate": 0.05,
            "subsample": 0.8,
        },
        file_size="28.7 MB",
    ),
    Model(
        name="DeltaKRR_v1.0",
        algorithm="Kernel Ridge Regression",
        mae=2.89,
        r2=0.931,
        training_date=date(2026, 1, 10),
        status="ready",
        hyperparameters={
            "alpha": 1.0,
            "kernel": "rbf",
            "gamma": 0.01,
        },
        file_size="156.3 MB",
    ),
    Model(
        name="DeltaNN_v1.0",
        algorithm="Neural Network",
        mae=None,
        r2=None,
        training_date=None,
        status="in_development",
        hyperparameters={
            "layers": "[256, 128, 64]",
            "activation": "ReLU",
            "optimizer": "Adam",
        },
        file_size=None,
    ),
    Model(
        name="DeltaEnsemble_v1.0",
        algorithm="Ensemble (RF + XGB + KRR)",
        mae=1.98,
        r2=0.961,
        training_date=date(2026, 1, 12),
        status="ready",
        hyperparameters={
            "weights": "[0.35, 0.45, 0.20]",
            "method": "weighted_average",
        },
        file_size="230.2 MB",
    ),
]

# Mock divergence statistics
DIVERGENCE_STATS: list[DivergenceStats] = [
    DivergenceStats(
        severity="LOW",
        range_min=0,
        range_max=10,
        count=22519,
        percentage=75.0,
    ),
    DivergenceStats(
        severity="MEDIUM",
        range_min=10,
        range_max=25,
        count=5405,
        percentage=18.0,
    ),
    DivergenceStats(
        severity="HIGH",
        range_min=25,
        range_max=50,
        count=1501,
        percentage=5.0,
    ),
    DivergenceStats(
        severity="CRITICAL",
        range_min=50,
        range_max=100,
        count=601,
        percentage=2.0,
    ),
]

# Common molecule names for mock predictions
MOLECULE_NAMES: dict[str, str] = {
    "CCO": "Ethanol",
    "CC(=O)O": "Acetic Acid",
    "c1ccccc1": "Benzene",
    "CC": "Ethane",
    "C": "Methane",
    "CO": "Methanol",
    "CC(C)C": "Isobutane",
    "C1CCCCC1": "Cyclohexane",
    "CC(=O)C": "Acetone",
    "CCN": "Ethylamine",
    "c1ccc(O)cc1": "Phenol",
    "CC(C)O": "Isopropanol",
}

# Default model for predictions
DEFAULT_MODEL = "DeltaXGB_v1.0"


def get_molecule_name(smiles: str) -> str:
    """Get molecule name from SMILES or return 'Unknown'."""
    return MOLECULE_NAMES.get(smiles, "Unknown Molecule")


def mock_predict(
    smiles: str,
    model_name: str = DEFAULT_MODEL,
) -> PredictionResult:
    """
    Generate a mock prediction for a given SMILES string.

    In MVP, this returns plausible mock values.
    Will be replaced with real model inference later.
    """
    # Generate deterministic but varied values based on SMILES
    # Use a local RNG to avoid mutating global random state
    seed = sum(ord(c) for c in smiles)
    rnd = Random(seed)

    predicted_value = rnd.uniform(-100, 50)
    confidence = rnd.uniform(0.85, 0.99)

    return PredictionResult(
        smiles=smiles,
        molecule_name=get_molecule_name(smiles),
        property_name="Heat of Formation (HOF)",
        predicted_value=round(predicted_value, 2),
        unit="kcal/mol",
        confidence=round(confidence, 3),
        model_used=model_name,
    )


def get_database_by_name(name: str) -> Database | None:
    """Get a database by name."""
    for db in DATABASES:
        if db.name == name:
            return db
    return None


def get_model_by_name(name: str) -> Model | None:
    """Get a model by name."""
    for model in MODELS:
        if model.name == name:
            return model
    return None


def get_ready_models() -> list[Model]:
    """Get all models that are ready for use."""
    return [m for m in MODELS if m.status == "ready"]


def get_ready_databases() -> list[Database]:
    """Get all databases that are ready for use."""
    return [db for db in DATABASES if db.status == "ready"]


# System status info
SYSTEM_INFO = {
    "version": "1.0.0-beta",
    "build_date": date(2026, 1, 12),
    "python_version": "3.10+",
    "databases_ready": len(get_ready_databases()),
    "databases_total": len(DATABASES),
    "models_ready": len(get_ready_models()),
    "models_total": len(MODELS),
    "default_model": DEFAULT_MODEL,
}
