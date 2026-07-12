"""Tests for typed CLI session context."""

from __future__ import annotations

from pathlib import Path

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.controller import CliController
from grimperium.cli.session import (
    AnalysisSourceRef,
    DatasetRef,
    MethodExecutionOverrides,
    ModelState,
    SessionContext,
)


def test_session_context_derives_analysis_path_from_dataset() -> None:
    dataset = DatasetRef(
        database_id="official.crest_pm7",
        alias="PM7",
        name="CREST PM7",
        path=Path("data/thermo_pm7.csv"),
        role="crest_pm7",
        capabilities=frozenset({"readable", "analysis_input"}),
    )
    session = SessionContext(dataset=dataset)

    assert session.analysis_path == Path("data/thermo_pm7.csv")

    session.analysis_source = AnalysisSourceRef(
        source_type="csv",
        name="adhoc",
        path=Path("tmp/adhoc.csv"),
    )
    assert session.analysis_path == Path("tmp/adhoc.csv")


def test_method_execution_overrides_are_schema_only() -> None:
    overrides = MethodExecutionOverrides(
        n_conformers=3,
        crest_timeout_minutes=10.0,
    )

    assert overrides.n_conformers == 3
    assert overrides.mopac_timeout_minutes is None


def test_controller_current_csv_path_is_derived_from_session_dataset() -> None:
    ctrl = CliController()
    dataset = DatasetRef(
        database_id="user.123",
        alias="USER",
        name="User DB",
        path=Path("data/user.csv"),
        role="analysis",
        capabilities=frozenset({"readable", "analysis_input"}),
    )

    ctrl.set_dataset(dataset)
    assert ctrl.current_csv_path == Path("data/user.csv")
    assert ctrl.session_summary()["dataset"] == "USER"

    ctrl.set_csv(Path("data/adhoc.csv"))
    assert ctrl.current_csv_path == Path("data/adhoc.csv")
    assert ctrl.session_summary()["dataset"] == "Ad-hoc: adhoc"

    ctrl.set_csv(None)
    assert ctrl.current_csv_path == Path("data/user.csv")


def test_model_state_reconciles_when_method_changes(tmp_path: Path) -> None:
    ctrl = CliController()
    missing_model = tmp_path / "missing.joblib"
    ctrl.set_model("stale", missing_model)

    method_b = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    ctrl.set_method(method_b)
    assert ctrl.session.model.state is ModelState.MISSING
    assert ctrl.status == "Model required"
    assert ctrl.session_summary()["model"] == "Missing: stale"

    method_a = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )
    ctrl.set_method(method_a)
    assert ctrl.session.model.state is ModelState.NOT_REQUIRED
    assert ctrl.session_summary()["model"] == "Not required"
