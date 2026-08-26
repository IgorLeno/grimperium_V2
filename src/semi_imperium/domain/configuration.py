"""Effective configuration and its reproducible calculation signature.

Reuse is the whole point of this module. A stored result may only be
reused when the *effective* configuration that produced it is identical
to the one now being requested, which means the signature has to cover
every knob that can change the number:

* the conformer search (CREST),
* the conformer selection strategy,
* the semiempirical Hamiltonian and its MOPAC settings,
* the minimum-verification policy.

Anything outside that list (thread counts, timeouts, output directories)
is execution detail: it changes how long the calculation takes, not what
it computes, so it is deliberately kept out of the digest.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from semi_imperium.domain.enums import VerificationPolicy
from semi_imperium.domain.hashing import DIGEST_ALGORITHM, stable_digest

if TYPE_CHECKING:  # pragma: no cover - typing only
    from grimperium.crest_pm7.config import PM7Config

#: Bumping this invalidates every previously computed signature on purpose.
#: Do it whenever the meaning of a signature field changes.
SIGNATURE_VERSION = 1

SIGNATURE_ID_LENGTH = 16


@dataclass(frozen=True)
class ConformerSearchSettings:
    """CREST conformer search settings that affect the sampled ensemble."""

    program: str = "crest"
    method: str = "gfn2"
    quick_mode: str = "off"
    opt_level: int = 2
    rmsd_threshold: float = 0.125
    energy_window_kcal_mol: float = 6.0
    max_conformers: int = 10
    nci: bool = False
    use_v3: bool = True
    preoptimizer: str = "xtb"
    """Pre-optimizer name, or ``"none"`` when pre-optimization is disabled."""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "program": self.program,
            "method": self.method,
            "quick_mode": self.quick_mode,
            "opt_level": self.opt_level,
            "rmsd_threshold": self.rmsd_threshold,
            "energy_window_kcal_mol": self.energy_window_kcal_mol,
            "max_conformers": self.max_conformers,
            "nci": self.nci,
            "use_v3": self.use_v3,
            "preoptimizer": self.preoptimizer,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ConformerSearchSettings:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            program=str(payload["program"]),
            method=str(payload["method"]),
            quick_mode=str(payload["quick_mode"]),
            opt_level=int(payload["opt_level"]),
            rmsd_threshold=float(payload["rmsd_threshold"]),
            energy_window_kcal_mol=float(payload["energy_window_kcal_mol"]),
            max_conformers=int(payload["max_conformers"]),
            nci=bool(payload["nci"]),
            use_v3=bool(payload["use_v3"]),
            preoptimizer=str(payload["preoptimizer"]),
        )


@dataclass(frozen=True)
class ConformerSelectionSettings:
    """Which conformers from the search are carried into the final result."""

    strategy: str = "lowest_pm7_hof_within_crest_subset"
    subset_size: int = 3
    energy_window_kcal_mol: float | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "strategy": self.strategy,
            "subset_size": self.subset_size,
            "energy_window_kcal_mol": self.energy_window_kcal_mol,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ConformerSelectionSettings:
        """Deserialize from JSON-compatible primitives."""
        window = payload.get("energy_window_kcal_mol")
        return cls(
            strategy=str(payload["strategy"]),
            subset_size=int(payload["subset_size"]),
            energy_window_kcal_mol=None if window is None else float(window),
        )


@dataclass(frozen=True)
class SemiempiricalSettings:
    """MOPAC Hamiltonian and the settings that change the reported energy."""

    program: str = "mopac"
    hamiltonian: str = "PM7"
    keywords: tuple[str, ...] = ("PRECISE",)
    scf_convergence: float = 1.0e-4
    solvent: str = "gas"
    """Solvent model name; ``"gas"`` means gas phase, never an empty value."""

    def __post_init__(self) -> None:
        if not self.hamiltonian.strip():
            raise ValueError("SemiempiricalSettings.hamiltonian must not be empty")
        # Keyword order is an input-file detail, not a scientific difference.
        object.__setattr__(
            self, "keywords", tuple(sorted(word.upper() for word in self.keywords))
        )
        object.__setattr__(self, "hamiltonian", self.hamiltonian.strip().upper())

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "program": self.program,
            "hamiltonian": self.hamiltonian,
            "keywords": list(self.keywords),
            "scf_convergence": self.scf_convergence,
            "solvent": self.solvent,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> SemiempiricalSettings:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            program=str(payload["program"]),
            hamiltonian=str(payload["hamiltonian"]),
            keywords=tuple(str(word) for word in payload["keywords"]),
            scf_convergence=float(payload["scf_convergence"]),
            solvent=str(payload["solvent"]),
        )


@dataclass(frozen=True)
class VerificationSettings:
    """The minimum-verification policy applied to every calculation."""

    policy: VerificationPolicy = VerificationPolicy.NONE
    imaginary_frequency_threshold_cm1: float = -10.0
    """Frequencies below this value count as imaginary (saddle evidence)."""

    def __post_init__(self) -> None:
        if self.imaginary_frequency_threshold_cm1 > 0:
            raise ValueError(
                "VerificationSettings.imaginary_frequency_threshold_cm1 must be "
                f"<= 0, got {self.imaginary_frequency_threshold_cm1}"
            )

    @property
    def requires_minimum(self) -> bool:
        """Whether a result is only acceptable when proven to be a minimum."""
        return self.policy is VerificationPolicy.REQUIRE_MINIMUM

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "policy": self.policy.value,
            "imaginary_frequency_threshold_cm1": (
                self.imaginary_frequency_threshold_cm1
            ),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> VerificationSettings:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            policy=VerificationPolicy(str(payload["policy"])),
            imaginary_frequency_threshold_cm1=float(
                payload["imaginary_frequency_threshold_cm1"]
            ),
        )


@dataclass(frozen=True)
class CalculationSignature:
    """A reproducible digest of one effective configuration.

    Two signatures are equal exactly when the two configurations would
    compute the same number, which is the precondition for reuse.
    """

    digest: str
    algorithm: str = DIGEST_ALGORITHM
    version: int = SIGNATURE_VERSION

    def __post_init__(self) -> None:
        if not self.digest.strip():
            raise ValueError("CalculationSignature.digest must not be empty")

    @property
    def short(self) -> str:
        """Short, human-quotable form used in filesystem paths and logs."""
        return self.digest[:SIGNATURE_ID_LENGTH]

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "digest": self.digest,
            "algorithm": self.algorithm,
            "version": self.version,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> CalculationSignature:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            digest=str(payload["digest"]),
            algorithm=str(payload["algorithm"]),
            version=int(payload["version"]),
        )


@dataclass(frozen=True)
class EffectiveConfiguration:
    """The complete, resolved configuration behind one calculation.

    "Effective" means every default has already been applied: nothing in
    here is inherited at read time, so a persisted configuration can be
    replayed years later without consulting today's defaults.
    """

    method_id: str
    method_version: str
    property_id: str
    conformer_search: ConformerSearchSettings = field(
        default_factory=ConformerSearchSettings
    )
    conformer_selection: ConformerSelectionSettings = field(
        default_factory=ConformerSelectionSettings
    )
    semiempirical: SemiempiricalSettings = field(default_factory=SemiempiricalSettings)
    verification: VerificationSettings = field(default_factory=VerificationSettings)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "method_id": self.method_id,
            "method_version": self.method_version,
            "property_id": self.property_id,
            "conformer_search": self.conformer_search.to_dict(),
            "conformer_selection": self.conformer_selection.to_dict(),
            "semiempirical": self.semiempirical.to_dict(),
            "verification": self.verification.to_dict(),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> EffectiveConfiguration:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            method_id=str(payload["method_id"]),
            method_version=str(payload["method_version"]),
            property_id=str(payload["property_id"]),
            conformer_search=ConformerSearchSettings.from_dict(
                payload["conformer_search"]
            ),
            conformer_selection=ConformerSelectionSettings.from_dict(
                payload["conformer_selection"]
            ),
            semiempirical=SemiempiricalSettings.from_dict(payload["semiempirical"]),
            verification=VerificationSettings.from_dict(payload["verification"]),
        )

    def signature_payload(self) -> dict[str, Any]:
        """Return exactly the fields that participate in the signature."""
        payload = self.to_dict()
        payload["signature_version"] = SIGNATURE_VERSION
        return payload

    def signature(self) -> CalculationSignature:
        """Return the reproducible signature of this configuration."""
        return CalculationSignature(
            digest=stable_digest(self.signature_payload()),
            algorithm=DIGEST_ALGORITHM,
            version=SIGNATURE_VERSION,
        )

    @classmethod
    def from_pm7_config(
        cls,
        config: PM7Config,
        *,
        method_id: str = "crest_pm7",
        method_version: str = "1.0",
        property_id: str = "standard_enthalpy_of_formation",
        verification: VerificationSettings | None = None,
        conformer_selection: ConformerSelectionSettings | None = None,
    ) -> EffectiveConfiguration:
        """Project Grimperium's :class:`PM7Config` onto this model.

        This reads Grimperium's configuration; it never writes to it and
        never reinterprets its stored scientific data. Execution-only
        fields (threads, timeouts, paths, monitoring thresholds) are
        intentionally dropped because they cannot change the result.
        """
        keywords = ("PRECISE",) if config.mopac_precise_scf else ()
        return cls(
            method_id=method_id,
            method_version=method_version,
            property_id=property_id,
            conformer_search=ConformerSearchSettings(
                method=config.crest_method,
                quick_mode=config.crest_quick_mode,
                opt_level=config.crest_opt_level,
                rmsd_threshold=config.crest_rmsd_threshold,
                energy_window_kcal_mol=config.energy_window,
                max_conformers=config.max_conformers,
                nci=config.crest_nci,
                use_v3=config.crest_v3,
                preoptimizer="xtb" if config.xtb_preopt else "none",
            ),
            conformer_selection=(
                conformer_selection
                if conformer_selection is not None
                else ConformerSelectionSettings()
            ),
            semiempirical=SemiempiricalSettings(
                keywords=keywords,
                scf_convergence=config.mopac_scf_threshold,
            ),
            verification=(
                verification if verification is not None else VerificationSettings()
            ),
        )


__all__ = [
    "SIGNATURE_VERSION",
    "CalculationSignature",
    "ConformerSearchSettings",
    "ConformerSelectionSettings",
    "EffectiveConfiguration",
    "SemiempiricalSettings",
    "VerificationSettings",
]
