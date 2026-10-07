"""Read-only dashboard API on small fixture files written per test."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from grimperium.webapi.app import create_app
from grimperium.webapi.config import DashboardPaths
from semi_imperium.domain import (
    CalculationRecord,
    ConformerSearchSettings,
    ConformerSelectionSettings,
    EffectiveConfiguration,
    MolecularIdentity,
    RunRecord,
    ScientificProvenance,
    SemiempiricalSettings,
    Timestamps,
    VerificationSettings,
)
from semi_imperium.persistence import SemiImperiumStore

EXP_ROWS = {
    "source": ["atct", "nist:gas", "bains:Pedley", "bains:Yaws"],
    "uncertainty_kind": ["reported", "reported", "none", "none"],
    "conflict_flag": ["False", "True", "False", "True"],
    "validation_only": ["False", "True", "False", "False"],
    "low_priority_only": ["False", "False", "False", "True"],
    "stewart_reference": ["False", "False", "True", "False"],
}

COMPARISON_ROWS = {
    "exp_id": ["exp_00000", "exp_00001"],
    "mol_id": ["cbs_1", "cbs_2"],
    "inchikey": ["VNWKTOKETHGBQD-UHFFFAOYSA-N", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"],
    "smiles": ["C", "CCO"],
    "nheavy": [1, 3],
    "source": ["atct", "nist:gas"],
    "H298_exp": [-17.8, -56.2],
    "uncertainty": [0.01, None],
    "H298_cbs": [-16.8, -53.2],
    "cbs_minus_exp": [1.0, 3.0],
    "conflict_flag": [False, True],
    "low_priority_only": [False, False],
    "validation_only": [False, True],
}


@pytest.fixture
def paths(tmp_path: Path) -> DashboardPaths:
    pd.DataFrame(
        {
            "mol_id": ["a", "b", "c", "d"],
            "status": ["OK", "OK", "Pending", "Skip"],
            "cbs_quality_flag": ["OK", "SUSPECT", "", ""],
        }
    ).to_csv(tmp_path / "pm7.csv", index=False)
    pd.DataFrame(EXP_ROWS).to_csv(tmp_path / "exp.csv", index=False)
    pd.DataFrame(COMPARISON_ROWS).to_csv(tmp_path / "vs_cbs.csv", index=False)
    return DashboardPaths(
        pm7_csv=tmp_path / "pm7.csv",
        semi_imperium_store=tmp_path / "store",
        exp_csv=tmp_path / "exp.csv",
        exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
    )


def _client(paths: DashboardPaths) -> TestClient:
    return TestClient(create_app(paths))


def _seed_store(root: Path) -> None:
    """One run with one pending PM7 calculation (the summary needs the run)."""
    configuration = EffectiveConfiguration(
        method_id="crest_pm7",
        method_version="1.0",
        property_id="standard_enthalpy_of_formation",
        conformer_search=ConformerSearchSettings(),
        conformer_selection=ConformerSelectionSettings(),
        semiempirical=SemiempiricalSettings(),
        verification=VerificationSettings(),
    )
    provenance = ScientificProvenance(
        method_id="crest_pm7",
        method_version="1.0",
        property_id="standard_enthalpy_of_formation",
        semi_imperium_version="0.2.0",
        grimperium_version="0.2.0",
    )
    created = Timestamps(created_at=datetime(2026, 10, 7, tzinfo=timezone.utc))
    store = SemiImperiumStore(root)
    store.save_run(
        RunRecord(
            run_id="run-0001",
            configuration=configuration,
            provenance=provenance,
            timestamps=created,
        )
    )
    store.save_calculation(
        CalculationRecord(
            run_id="run-0001",
            molecule=MolecularIdentity(
                canonical_smiles="CCO", charge=0, multiplicity=1
            ),
            signature=configuration.signature(),
            provenance=provenance,
            timestamps=created,
        )
    )


def test_health_reports_read_only(paths: DashboardPaths) -> None:
    response = _client(paths).get("/api/health")
    assert response.status_code == 200
    assert response.json() == {"status": "ok", "read_only": True}


def test_progress_counts_legacy_rows_and_labels_protocol(
    paths: DashboardPaths,
) -> None:
    body = _client(paths).get("/api/dataset/progress").json()

    legacy = body["legacy"]
    assert legacy["protocol"] == "crest_pm7"
    assert legacy["total"] == 4
    assert legacy["status"] == {"OK": 2, "Pending": 1, "Skip": 1}
    assert legacy["cbs_quality_flag"] == {"unset": 2, "OK": 1, "SUSPECT": 1}


def test_progress_treats_missing_store_as_empty(paths: DashboardPaths) -> None:
    store = _client(paths).get("/api/dataset/progress").json()["semi_imperium"]

    assert store["exists"] is False
    assert (store["molecules"], store["runs"], store["calculations"]) == (0, 0, 0)
    assert set(store["states"].values()) == {0}
    assert not paths.semi_imperium_store.exists()


def test_progress_summarises_store_records(paths: DashboardPaths) -> None:
    _seed_store(paths.semi_imperium_store)
    before = sorted(p.stat().st_mtime_ns for p in paths.semi_imperium_store.rglob("*"))

    store = _client(paths).get("/api/dataset/progress").json()["semi_imperium"]

    assert store["exists"] is True
    assert (store["molecules"], store["runs"], store["calculations"]) == (1, 1, 1)
    assert store["states"]["pending"] == 1
    after = sorted(p.stat().st_mtime_ns for p in paths.semi_imperium_store.rglob("*"))
    assert after == before


def test_experimental_summary_counts_sources_and_flags(paths: DashboardPaths) -> None:
    body = _client(paths).get("/api/experimental/summary").json()

    assert body["rows"] == 4
    assert body["by_source"]["nist:gas"] == 1
    assert body["by_source_family"] == {"bains": 2, "atct": 1, "nist": 1}
    assert body["flags"] == {
        "conflict_flag": 2,
        "validation_only": 1,
        "low_priority_only": 1,
        "stewart_reference": 1,
    }
    assert body["uncertainty_kind"] == {"reported": 2, "none": 2}


def test_experimental_comparison_returns_points_and_stats(
    paths: DashboardPaths,
) -> None:
    body = _client(paths).get("/api/experimental/comparison").json()

    assert [p["exp_id"] for p in body["points"]] == ["exp_00000", "exp_00001"]
    methane, ethanol = body["points"]
    assert methane["uncertainty"] == pytest.approx(0.01)
    assert ethanol["uncertainty"] is None
    assert ethanol["validation_only"] is True and ethanol["conflict_flag"] is True
    assert ethanol["nheavy"] == 3
    stats = body["stats"]
    assert stats["n"] == 2
    assert stats["bias"] == pytest.approx(2.0)
    assert stats["mae"] == pytest.approx(2.0)
    assert stats["rmse"] == pytest.approx((5.0) ** 0.5)


@pytest.mark.parametrize(
    ("endpoint", "attribute"),
    [
        ("/api/dataset/progress", "pm7_csv"),
        ("/api/experimental/summary", "exp_csv"),
        ("/api/experimental/comparison", "exp_vs_cbs_csv"),
    ],
)
def test_missing_file_is_an_explicit_503(
    paths: DashboardPaths, endpoint: str, attribute: str
) -> None:
    missing = getattr(paths, attribute)
    missing.unlink()

    response = _client(paths).get(endpoint)

    assert response.status_code == 503
    assert response.json() == {"detail": "file not found", "path": str(missing)}


def test_missing_column_is_reported(paths: DashboardPaths) -> None:
    pd.DataFrame({"status": ["OK"]}).to_csv(paths.pm7_csv, index=False)

    response = _client(paths).get("/api/dataset/progress")

    assert response.status_code == 503
    assert "cbs_quality_flag" in response.json()["detail"]


def test_non_boolean_flag_is_reported(paths: DashboardPaths) -> None:
    rows = dict(EXP_ROWS, conflict_flag=["False", "maybe", "False", "True"])
    pd.DataFrame(rows).to_csv(paths.exp_csv, index=False)

    response = _client(paths).get("/api/experimental/summary")

    assert response.status_code == 503
    assert "conflict_flag" in response.json()["detail"]


def test_non_numeric_value_is_reported(paths: DashboardPaths) -> None:
    rows = dict(COMPARISON_ROWS, H298_cbs=["-16.8", "n/a"])
    pd.DataFrame(rows).to_csv(paths.exp_vs_cbs_csv, index=False)

    response = _client(paths).get("/api/experimental/comparison")

    assert response.status_code == 503
    assert "H298_cbs" in response.json()["detail"]


def test_only_get_routes_are_exposed(paths: DashboardPaths) -> None:
    app = create_app(paths)
    methods = {
        method
        for route in app.routes
        if getattr(route, "path", "").startswith("/api/")
        for method in getattr(route, "methods", set())
    }
    assert methods <= {"GET", "HEAD"}
