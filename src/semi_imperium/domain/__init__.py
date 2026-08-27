"""Traceable domain model for Semi-Imperium calculations.

Public surface of the domain layer: molecular identity, effective
configuration and its signature, explicit lifecycle states, and the
records that carry timestamps and scientific provenance.
"""

from semi_imperium.domain.configuration import (
    SIGNATURE_VERSION,
    CalculationSignature,
    ConformerSearchSettings,
    ConformerSelectionSettings,
    EffectiveConfiguration,
    SemiempiricalSettings,
    VerificationSettings,
)
from semi_imperium.domain.enums import (
    ALLOWED_CALCULATION_TRANSITIONS,
    ALLOWED_RUN_TRANSITIONS,
    COHERENT_VERIFICATION_OUTCOMES,
    DEFAULT_REUSABLE_STATES,
    RESULT_CALCULATION_STATES,
    TERMINAL_CALCULATION_STATES,
    TERMINAL_RUN_STATES,
    VERIFIED_ONLY_REUSABLE_STATES,
    CalculationState,
    RunState,
    VerificationOutcome,
    VerificationPolicy,
)
from semi_imperium.domain.hashing import DIGEST_ALGORITHM, stable_digest
from semi_imperium.domain.identity import MolecularIdentity, MoleculeInputType
from semi_imperium.domain.records import (
    RECORD_SCHEMA_VERSION,
    CalculationRecord,
    CalculationResultData,
    RunRecord,
    ScientificProvenance,
    Timestamps,
    build_calculation_id,
    build_reuse_key,
    utc_now,
)

__all__ = [
    "ALLOWED_CALCULATION_TRANSITIONS",
    "ALLOWED_RUN_TRANSITIONS",
    "COHERENT_VERIFICATION_OUTCOMES",
    "DEFAULT_REUSABLE_STATES",
    "DIGEST_ALGORITHM",
    "RECORD_SCHEMA_VERSION",
    "RESULT_CALCULATION_STATES",
    "SIGNATURE_VERSION",
    "TERMINAL_CALCULATION_STATES",
    "TERMINAL_RUN_STATES",
    "VERIFIED_ONLY_REUSABLE_STATES",
    "CalculationRecord",
    "CalculationResultData",
    "CalculationSignature",
    "CalculationState",
    "ConformerSearchSettings",
    "ConformerSelectionSettings",
    "EffectiveConfiguration",
    "MolecularIdentity",
    "MoleculeInputType",
    "RunRecord",
    "RunState",
    "ScientificProvenance",
    "SemiempiricalSettings",
    "Timestamps",
    "VerificationOutcome",
    "VerificationPolicy",
    "VerificationSettings",
    "build_calculation_id",
    "build_reuse_key",
    "stable_digest",
    "utc_now",
]
