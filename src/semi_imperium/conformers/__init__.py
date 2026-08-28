"""Optional conformer search and bounded conformer selection.

Public surface of the conformer stage: the ensemble types, the adapter
protocols that hold CREST and CONFPASS at arm's length, the two
selection strategies, and the workflow that chooses between the CREST
route and the initial-3D route.

:class:`semi_imperium.conformers.initial_structure.RDKitInitialStructure`
is imported from its own module on purpose: it is the only piece that
pulls RDKit in, and the rest of the stage stays free of it.
"""

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConformerRequest,
    ConformerSearchBackend,
    ConfPassBackend,
    ConfPassCandidate,
    ConfPassRanking,
    InitialStructureBackend,
)
from semi_imperium.conformers.confpass import (
    AdaptedStructure,
    ConfPassSelector,
    MoleculeTopology,
    UnavailableConfPass,
    build_confpass_candidates,
    read_sd_record,
    to_sd_record,
)
from semi_imperium.conformers.crest import (
    HARTREE_TO_KCAL_MOL,
    CrestConformerSearch,
    CrestRun,
    CrestRunner,
    parse_crest_ensemble,
)
from semi_imperium.conformers.ensemble import (
    UNKNOWN_PROGRAM_VERSION,
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerSearchProvenance,
)
from semi_imperium.conformers.selection import (
    FORBIDDEN_EVIDENCE_TERMS,
    PAS_COMPLETENESS_LABEL_KEY,
    ConformerSelector,
    EnergyTopNSelector,
    SelectionResult,
)
from semi_imperium.conformers.workflow import ConformerPreparation, ConformerWorkflow

__all__ = [
    "FORBIDDEN_EVIDENCE_TERMS",
    "HARTREE_TO_KCAL_MOL",
    "PAS_COMPLETENESS_LABEL_KEY",
    "UNKNOWN_PROGRAM_VERSION",
    "AdaptedStructure",
    "ConfPassBackend",
    "ConfPassCandidate",
    "ConfPassRanking",
    "ConfPassSelector",
    "Conformer",
    "ConformerBackendError",
    "ConformerEnsemble",
    "ConformerGeometry",
    "ConformerPreparation",
    "ConformerRequest",
    "ConformerSearchBackend",
    "ConformerSearchProvenance",
    "ConformerSelector",
    "ConformerWorkflow",
    "CrestConformerSearch",
    "CrestRun",
    "CrestRunner",
    "EnergyTopNSelector",
    "InitialStructureBackend",
    "MoleculeTopology",
    "SelectionResult",
    "UnavailableConfPass",
    "build_confpass_candidates",
    "parse_crest_ensemble",
    "read_sd_record",
    "to_sd_record",
]
