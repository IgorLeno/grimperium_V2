"""Contract for Semi-Imperium's traceable domain and persistence model.

The tests are grouped by the three properties the model has to hold:

1. identity, timestamps and provenance survive a round trip and never
   touch Grimperium's own scientific data;
2. reuse is decided by molecular identity plus a reproducible signature
   that separates CREST, conformer selection, MOPAC and verification;
3. every lifecycle state is explicit — no blank cell ever stands in for
   "pending", "unverified" or "saddle".
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from grimperium.crest_pm7.config import PM7Config
from semi_imperium.domain import (
    DEFAULT_REUSABLE_STATES,
    VERIFIED_ONLY_REUSABLE_STATES,
    CalculationRecord,
    CalculationResultData,
    CalculationState,
    ConformerSearchSettings,
    ConformerSelectionSettings,
    EffectiveConfiguration,
    MolecularIdentity,
    RunRecord,
    RunState,
    ScientificProvenance,
    SemiempiricalSettings,
    Timestamps,
    VerificationOutcome,
    VerificationPolicy,
    VerificationSettings,
    build_calculation_id,
)
from semi_imperium.persistence import SemiImperiumStore, StoreIntegrityError

UTC = timezone.utc
T0 = datetime(2026, 8, 26, 12, 0, 0, tzinfo=UTC)


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def make_molecule(
    canonical_smiles: str = "CCO", **overrides: object
) -> MolecularIdentity:
    """Build a molecular identity without needing RDKit."""
    fields: dict[str, object] = {
        "canonical_smiles": canonical_smiles,
        "charge": 0,
        "multiplicity": 1,
    }
    fields.update(overrides)
    return MolecularIdentity(**fields)  # type: ignore[arg-type]


def make_configuration(**overrides: object) -> EffectiveConfiguration:
    """Build a fully resolved effective configuration."""
    fields: dict[str, object] = {
        "method_id": "crest_pm7",
        "method_version": "1.0",
        "property_id": "standard_enthalpy_of_formation",
        "conformer_search": ConformerSearchSettings(),
        "conformer_selection": ConformerSelectionSettings(),
        "semiempirical": SemiempiricalSettings(),
        "verification": VerificationSettings(),
    }
    fields.update(overrides)
    return EffectiveConfiguration(**fields)  # type: ignore[arg-type]


def make_provenance(**overrides: object) -> ScientificProvenance:
    """Build provenance describing which software produced a number."""
    fields: dict[str, object] = {
        "method_id": "crest_pm7",
        "method_version": "1.0",
        "property_id": "standard_enthalpy_of_formation",
        "semi_imperium_version": "0.2.0",
        "grimperium_version": "0.2.0",
        "program_versions": {"crest": "3.0.1", "mopac": "23.1.2"},
        "hostname": "workstation-01",
    }
    fields.update(overrides)
    return ScientificProvenance(**fields)  # type: ignore[arg-type]


def make_record(
    *,
    run_id: str = "run-0001",
    molecule: MolecularIdentity | None = None,
    configuration: EffectiveConfiguration | None = None,
    created_at: datetime = T0,
) -> CalculationRecord:
    """Build a freshly registered (``PENDING``) calculation record."""
    configuration = configuration or make_configuration()
    return CalculationRecord(
        run_id=run_id,
        molecule=molecule or make_molecule(),
        signature=configuration.signature(),
        provenance=make_provenance(),
        timestamps=Timestamps(created_at=created_at),
    )


def make_result() -> CalculationResultData:
    """Build a plausible scientific payload for a finished calculation."""
    return CalculationResultData(
        energy_hof_kcal_mol=-56.12,
        conformer_index=2,
        conformers_evaluated=3,
        lowest_frequency_cm1=42.7,
        artifact_paths=("artifacts/mol/mopac.out",),
    )


def finish_verified(record: CalculationRecord) -> CalculationRecord:
    """Drive a record through RUNNING to a confirmed minimum."""
    running = record.transition_to(CalculationState.RUNNING, at=T0)
    return running.transition_to(
        CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        result=make_result(),
        at=T0 + timedelta(minutes=5),
    )


# ---------------------------------------------------------------------------
# 1. Identity, timestamps and provenance
# ---------------------------------------------------------------------------


def test_molecule_id_ignores_labels_but_tracks_charge_and_multiplicity() -> None:
    base = make_molecule()
    relabelled = make_molecule(display_name="ethanol", inchikey="ABCDEFGHIJKLMN-X")
    cation = make_molecule(charge=1)
    triplet = make_molecule(multiplicity=3)

    assert base.molecule_id == relabelled.molecule_id
    assert base.molecule_id != cation.molecule_id
    assert base.molecule_id != triplet.molecule_id
    assert base.molecule_id.startswith("mol-")


def test_molecule_id_is_canonical_across_smiles_spellings() -> None:
    pytest.importorskip("rdkit")

    written_one_way = MolecularIdentity.from_smiles("OCC")
    written_another_way = MolecularIdentity.from_smiles("CCO", display_name="ethanol")

    assert written_one_way.canonical_smiles == written_another_way.canonical_smiles
    assert written_one_way.molecule_id == written_another_way.molecule_id


def test_calculation_identity_is_deterministic_per_run_molecule_and_signature() -> None:
    configuration = make_configuration()
    record = make_record(configuration=configuration)

    expected = build_calculation_id(
        run_id="run-0001",
        molecule_id=record.molecule.molecule_id,
        signature_digest=configuration.signature().digest,
    )

    assert record.calculation_id == expected
    assert record.reuse_key == (
        f"{record.molecule.molecule_id}/{configuration.signature().digest}"
    )
    # A second run over the same molecule is a different calculation.
    assert make_record(run_id="run-0002").calculation_id != record.calculation_id


def test_calculation_id_mismatch_is_rejected() -> None:
    record = make_record()

    with pytest.raises(ValueError, match="calculation_id does not match"):
        CalculationRecord(
            run_id=record.run_id,
            molecule=record.molecule,
            signature=record.signature,
            provenance=record.provenance,
            calculation_id="calc-deadbeefdeadbeef",
        )


def test_timestamps_require_timezone_and_ordering() -> None:
    with pytest.raises(ValueError, match="timezone-aware"):
        Timestamps(created_at=datetime(2026, 8, 26, 12, 0, 0))

    with pytest.raises(ValueError, match="must not precede started_at"):
        Timestamps(
            created_at=T0,
            started_at=T0 + timedelta(minutes=5),
            completed_at=T0 + timedelta(minutes=1),
        )

    timestamps = Timestamps(
        created_at=T0, started_at=T0, completed_at=T0 + timedelta(seconds=90)
    )
    assert timestamps.duration_seconds == 90.0


def test_transitions_record_start_and_completion_timestamps() -> None:
    record = make_record()
    assert record.timestamps.started_at is None
    assert record.timestamps.completed_at is None

    finished = finish_verified(record)

    assert finished.timestamps.started_at == T0
    assert finished.timestamps.completed_at == T0 + timedelta(minutes=5)
    assert finished.timestamps.duration_seconds == 300.0


def test_calculation_round_trip_preserves_identity_and_provenance(
    tmp_path: Path,
) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    record = finish_verified(make_record())

    path = store.save_calculation(record)
    loaded = store.list_calculations(run_id=record.run_id)

    assert path.is_file()
    assert loaded == [record]
    assert loaded[0].provenance.program_versions == {
        "crest": "3.0.1",
        "mopac": "23.1.2",
    }
    assert loaded[0].result is not None
    assert loaded[0].result.energy_hof_kcal_mol == pytest.approx(-56.12)


def test_run_round_trip_preserves_effective_configuration(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    configuration = make_configuration(
        verification=VerificationSettings(policy=VerificationPolicy.REQUIRE_MINIMUM)
    )
    run = RunRecord(
        run_id="run-0001",
        configuration=configuration,
        provenance=make_provenance(),
        timestamps=Timestamps(created_at=T0),
        label="nightly",
        molecule_ids=(make_molecule().molecule_id,),
    )

    store.save_run(run)
    loaded = store.load_run("run-0001")

    assert loaded == run
    assert loaded.configuration.verification.requires_minimum is True
    assert loaded.signature.digest == configuration.signature().digest
    assert store.list_runs() == [loaded]


def test_store_writes_only_inside_its_own_root(tmp_path: Path) -> None:
    """Semi-Imperium never writes into Grimperium's data directories."""
    store_root = tmp_path / "store"
    untouched = tmp_path / "grimperium_data"
    untouched.mkdir()
    dataset = untouched / "thermo_pm7.csv"
    dataset.write_text("smiles,H298_pm7\nCCO,-56.0\n")
    before = dataset.read_text()

    store = SemiImperiumStore(store_root)
    run = RunRecord(
        run_id="run-0001",
        configuration=make_configuration(),
        provenance=make_provenance(),
        timestamps=Timestamps(created_at=T0),
    )
    store.save_run(run)
    store.save_calculation(finish_verified(make_record()))

    written = [
        path for path in tmp_path.rglob("*") if path.is_file() and path != dataset
    ]
    assert written
    assert all(path.is_relative_to(store_root) for path in written)
    assert dataset.read_text() == before
    assert sorted(p.name for p in store_root.iterdir()) == ["calculations", "runs"]


def test_tampered_molecule_identity_is_rejected_on_read(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    record = finish_verified(make_record())
    path = store.save_calculation(record)

    payload = json.loads(path.read_text())
    payload["molecule"]["canonical_smiles"] = "CCCO"
    path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="molecule_id does not match"):
        store.list_calculations()


def test_tampered_run_configuration_is_rejected_on_read(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    run = RunRecord(
        run_id="run-0001",
        configuration=make_configuration(),
        provenance=make_provenance(),
        timestamps=Timestamps(created_at=T0),
    )
    path = store.save_run(run)

    payload = json.loads(path.read_text())
    payload["configuration"]["semiempirical"]["hamiltonian"] = "AM1"
    path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="signature does not match"):
        store.load_run("run-0001")


# ---------------------------------------------------------------------------
# 2. Signature and reuse
# ---------------------------------------------------------------------------


def test_signature_is_reproducible_for_equal_configurations() -> None:
    first = make_configuration().signature()
    second = make_configuration().signature()

    assert first == second
    assert first.digest == second.digest
    assert first.algorithm == "sha256"
    assert len(first.digest) == 64


def test_signature_separates_every_scientifically_relevant_group() -> None:
    base = make_configuration()
    variants = {
        "crest": make_configuration(
            conformer_search=ConformerSearchSettings(method="gfnff")
        ),
        "selection": make_configuration(
            conformer_selection=ConformerSelectionSettings(subset_size=5)
        ),
        "hamiltonian": make_configuration(
            semiempirical=SemiempiricalSettings(hamiltonian="AM1")
        ),
        "mopac_settings": make_configuration(
            semiempirical=SemiempiricalSettings(scf_convergence=1.0e-6)
        ),
        "verification": make_configuration(
            verification=VerificationSettings(
                policy=VerificationPolicy.REQUIRE_MINIMUM
            )
        ),
    }

    digests = {name: config.signature().digest for name, config in variants.items()}
    digests["base"] = base.signature().digest

    assert len(set(digests.values())) == len(digests), digests


def test_signature_ignores_mopac_keyword_ordering() -> None:
    one = make_configuration(
        semiempirical=SemiempiricalSettings(keywords=("PRECISE", "GNORM=0.01"))
    )
    other = make_configuration(
        semiempirical=SemiempiricalSettings(keywords=("gnorm=0.01", "precise"))
    )

    assert one.signature().digest == other.signature().digest


def test_signature_from_pm7_config_tracks_science_not_execution_details() -> None:
    baseline = EffectiveConfiguration.from_pm7_config(PM7Config())
    faster_machine = EffectiveConfiguration.from_pm7_config(
        PM7Config(crest_threads=16, crest_timeout=900.0, mopac_timeout_base=600.0)
    )
    stricter_science = EffectiveConfiguration.from_pm7_config(
        PM7Config(),
        verification=VerificationSettings(
            policy=VerificationPolicy.REQUIRE_MINIMUM
        ),
    )
    different_hamiltonian_settings = EffectiveConfiguration.from_pm7_config(
        PM7Config(mopac_precise_scf=False)
    )

    assert baseline.signature() == faster_machine.signature()
    assert baseline.signature() != stricter_science.signature()
    assert baseline.signature() != different_hamiltonian_settings.signature()
    assert baseline.semiempirical.hamiltonian == "PM7"


def test_reuse_requires_matching_molecule_and_signature(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    configuration = make_configuration()
    molecule = make_molecule()
    store.save_calculation(
        finish_verified(make_record(molecule=molecule, configuration=configuration))
    )

    same_request = store.find_reusable(molecule, configuration.signature())
    other_molecule = store.find_reusable(
        make_molecule("CCCO"), configuration.signature()
    )
    other_settings = store.find_reusable(
        molecule,
        make_configuration(
            semiempirical=SemiempiricalSettings(hamiltonian="AM1")
        ).signature(),
    )

    assert same_request.can_reuse is True
    assert same_request.record is not None
    assert same_request.record.state is CalculationState.VERIFIED
    assert other_molecule.can_reuse is False
    assert other_settings.can_reuse is False
    assert "no stored calculation" in other_settings.reason


def test_reuse_respects_the_minimum_verification_requirement(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    configuration = make_configuration()
    molecule = make_molecule()
    running = make_record(molecule=molecule, configuration=configuration).transition_to(
        CalculationState.RUNNING, at=T0
    )
    store.save_calculation(
        running.transition_to(
            CalculationState.UNVERIFIED,
            verification=VerificationOutcome.NOT_REQUESTED,
            result=make_result(),
            at=T0 + timedelta(minutes=1),
        )
    )

    lenient = store.find_reusable(
        molecule, configuration.signature(), accepted_states=DEFAULT_REUSABLE_STATES
    )
    strict = store.find_reusable(
        molecule,
        configuration.signature(),
        accepted_states=VERIFIED_ONLY_REUSABLE_STATES,
    )

    assert lenient.can_reuse is True
    assert strict.can_reuse is False


def test_failed_calculations_are_never_reused(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    configuration = make_configuration()
    molecule = make_molecule()
    running = make_record(molecule=molecule, configuration=configuration).transition_to(
        CalculationState.RUNNING, at=T0
    )
    store.save_calculation(
        running.transition_to(
            CalculationState.FAILED,
            error_message="MOPAC SCF did not converge",
            at=T0 + timedelta(minutes=2),
        )
    )

    decision = store.find_reusable(molecule, configuration.signature())

    assert decision.can_reuse is False
    assert store.list_calculations()[0].state is CalculationState.FAILED


def test_reuse_prefers_the_most_recently_completed_attempt(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    configuration = make_configuration()
    molecule = make_molecule()

    older = finish_verified(
        make_record(
            run_id="run-0001", molecule=molecule, configuration=configuration
        )
    )
    newer_base = make_record(
        run_id="run-0002",
        molecule=molecule,
        configuration=configuration,
        created_at=T0 + timedelta(days=1),
    )
    newer = newer_base.transition_to(
        CalculationState.RUNNING, at=T0 + timedelta(days=1)
    ).transition_to(
        CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        result=make_result(),
        at=T0 + timedelta(days=1, minutes=5),
    )
    store.save_calculation(older)
    store.save_calculation(newer)

    decision = store.find_reusable(molecule, configuration.signature())

    assert decision.record is not None
    assert decision.record.run_id == "run-0002"
    assert decision.reuse_key == older.reuse_key


# ---------------------------------------------------------------------------
# 3. Explicit states
# ---------------------------------------------------------------------------


def terminal_records() -> dict[CalculationState, CalculationRecord]:
    """One persisted-ready record per lifecycle state, all distinct runs."""
    pending = make_record(run_id="run-pending")
    running = make_record(run_id="run-running").transition_to(
        CalculationState.RUNNING, at=T0
    )
    return {
        CalculationState.PENDING: pending,
        CalculationState.RUNNING: running,
        CalculationState.VERIFIED: finish_verified(make_record(run_id="run-verified")),
        CalculationState.UNVERIFIED: make_record(run_id="run-unverified")
        .transition_to(CalculationState.RUNNING, at=T0)
        .transition_to(
            CalculationState.UNVERIFIED,
            verification=VerificationOutcome.INCONCLUSIVE,
            result=make_result(),
            at=T0 + timedelta(minutes=1),
        ),
        CalculationState.SADDLE: make_record(run_id="run-saddle")
        .transition_to(CalculationState.RUNNING, at=T0)
        .transition_to(
            CalculationState.SADDLE,
            verification=VerificationOutcome.SADDLE_POINT,
            result=CalculationResultData(
                energy_hof_kcal_mol=-51.0,
                conformer_index=0,
                conformers_evaluated=3,
                lowest_frequency_cm1=-145.2,
            ),
            at=T0 + timedelta(minutes=1),
        ),
        CalculationState.FAILED: make_record(run_id="run-failed")
        .transition_to(CalculationState.RUNNING, at=T0)
        .transition_to(
            CalculationState.FAILED,
            error_message="CREST timed out",
            at=T0 + timedelta(minutes=1),
        ),
    }


def test_every_lifecycle_state_persists_explicitly(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    records = terminal_records()
    assert set(records) == set(CalculationState)

    for state, record in records.items():
        path = store.save_calculation(record)
        payload = json.loads(path.read_text())

        assert payload["state"] == state.value
        assert payload["state"].strip() != ""
        assert payload["verification"].strip() != ""
        assert CalculationRecord.from_dict(payload).state is state


def test_saddle_is_a_result_state_not_a_failure(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    saddle = terminal_records()[CalculationState.SADDLE]

    store.save_calculation(saddle)
    stored = store.list_calculations(run_id="run-saddle")[0]

    assert stored.state.has_result is True
    assert stored.result is not None
    assert stored.result.lowest_frequency_cm1 == pytest.approx(-145.2)
    assert stored.verification is VerificationOutcome.SADDLE_POINT
    assert stored.error_message is None


def test_state_and_verification_must_agree() -> None:
    running = make_record().transition_to(CalculationState.RUNNING, at=T0)

    with pytest.raises(ValueError, match="incoherent"):
        running.transition_to(
            CalculationState.VERIFIED,
            verification=VerificationOutcome.SADDLE_POINT,
            result=make_result(),
            at=T0,
        )


def test_result_states_require_a_result_and_failures_require_a_reason() -> None:
    running = make_record().transition_to(CalculationState.RUNNING, at=T0)

    with pytest.raises(ValueError, match="must carry a result"):
        running.transition_to(
            CalculationState.VERIFIED,
            verification=VerificationOutcome.MINIMUM_CONFIRMED,
            at=T0,
        )

    with pytest.raises(ValueError, match="must carry an error_message"):
        running.transition_to(CalculationState.FAILED, at=T0)


def test_terminal_states_reject_further_transitions() -> None:
    finished = finish_verified(make_record())

    with pytest.raises(ValueError, match="Illegal calculation transition"):
        finished.transition_to(CalculationState.RUNNING, at=T0)


def test_store_refuses_to_rewrite_a_terminal_calculation(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    record = make_record()
    running = record.transition_to(CalculationState.RUNNING, at=T0)
    verified = running.transition_to(
        CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        result=make_result(),
        at=T0 + timedelta(minutes=5),
    )

    store.save_calculation(record)
    store.save_calculation(running)
    store.save_calculation(verified)
    store.save_calculation(verified)  # idempotent re-save is allowed

    contradiction = CalculationRecord(
        run_id=record.run_id,
        molecule=record.molecule,
        signature=record.signature,
        provenance=record.provenance,
        state=CalculationState.FAILED,
        verification=VerificationOutcome.FAILED,
        timestamps=verified.timestamps,
        error_message="rewritten after the fact",
    )
    with pytest.raises(StoreIntegrityError, match="already terminal"):
        store.save_calculation(contradiction)

    assert store.list_calculations()[0].state is CalculationState.VERIFIED


def test_store_refuses_to_rewrite_a_terminal_run(tmp_path: Path) -> None:
    store = SemiImperiumStore(tmp_path / "store")
    run = RunRecord(
        run_id="run-0001",
        configuration=make_configuration(),
        provenance=make_provenance(),
        timestamps=Timestamps(created_at=T0),
    )
    completed = run.transition_to(RunState.RUNNING, at=T0).transition_to(
        RunState.COMPLETED, at=T0 + timedelta(hours=1)
    )
    store.save_run(completed)

    reopened = RunRecord(
        run_id="run-0001",
        configuration=make_configuration(),
        provenance=make_provenance(),
        state=RunState.RUNNING,
        timestamps=Timestamps(created_at=T0),
    )
    with pytest.raises(StoreIntegrityError, match="already terminal"):
        store.save_run(reopened)

    assert store.load_run("run-0001").state is RunState.COMPLETED
