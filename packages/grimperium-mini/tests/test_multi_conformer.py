"""Tests for grimperium_mini.multi_conformer."""

from __future__ import annotations

from pathlib import Path

import pytest

from grimperium_mini.config import MiniConfig
from grimperium_mini.models import MopacResult
from grimperium_mini.multi_conformer import (
    _build_columns,
    _safe_name,
    process_molecule_multi,
    run_multi_conformer_pipeline,
    select_existing_conformers,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_config(tmp_path: Path) -> MiniConfig:
    return MiniConfig(
        work_root=tmp_path / "runs",
        results_dir=tmp_path / "results",
        mopac_executable="mopac",
    )


def _mopac_ok(hof_kJmol: float) -> MopacResult:
    return MopacResult(
        "success",
        hof_kcalmol=hof_kJmol / 4.184,
        hof_kJmol=hof_kJmol,
    )


def _mopac_fail() -> MopacResult:
    return MopacResult("failed", error_message="test failure")


def _make_conformer_files(run_dir: Path, n: int) -> list[Path]:
    """Create minimal valid XYZ conformer files in *run_dir*."""
    run_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for i in range(1, n + 1):
        p = run_dir / f"conformer_{i:04d}.xyz"
        p.write_text(
            "5\ncomment\nC  0.0  0.0  0.0\nH  1.0  0.0  0.0\n"
            "H -1.0  0.0  0.0\nH  0.0  1.0  0.0\nH  0.0 -1.0  0.0\n",
            encoding="utf-8",
        )
        paths.append(p)
    return paths


def _default_mol(hf_exp: float | None = -100.0) -> dict[str, object]:
    return {
        "mol_id": "test_001",
        "molecule": "Methane",
        "smiles": "C",
        "charge": 0,
        "multiplicity": 1,
        "hf_exp_kJmol": hf_exp,
    }


# ---------------------------------------------------------------------------
# select_existing_conformers
# ---------------------------------------------------------------------------


class TestSelectExistingConformers:
    def test_empty_when_dir_missing(self, tmp_path: Path) -> None:
        assert select_existing_conformers(tmp_path / "nonexistent") == []

    def test_empty_when_no_conformer_files(self, tmp_path: Path) -> None:
        (tmp_path / "other.xyz").write_text("x")
        assert select_existing_conformers(tmp_path) == []

    def test_returns_sorted_ascending(self, tmp_path: Path) -> None:
        for name in ["conformer_0003.xyz", "conformer_0001.xyz", "conformer_0002.xyz"]:
            (tmp_path / name).write_text("x")
        result = select_existing_conformers(tmp_path)
        assert [p.name for p in result] == [
            "conformer_0001.xyz",
            "conformer_0002.xyz",
            "conformer_0003.xyz",
        ]

    def test_limits_to_max_n(self, tmp_path: Path) -> None:
        for i in range(1, 15):
            (tmp_path / f"conformer_{i:04d}.xyz").write_text("x")
        assert len(select_existing_conformers(tmp_path, max_n=10)) == 10

    def test_returns_all_when_fewer_than_max(self, tmp_path: Path) -> None:
        for i in range(1, 4):
            (tmp_path / f"conformer_{i:04d}.xyz").write_text("x")
        assert len(select_existing_conformers(tmp_path, max_n=10)) == 3

    def test_custom_max_n(self, tmp_path: Path) -> None:
        for i in range(1, 6):
            (tmp_path / f"conformer_{i:04d}.xyz").write_text("x")
        assert len(select_existing_conformers(tmp_path, max_n=3)) == 3


# ---------------------------------------------------------------------------
# _safe_name
# ---------------------------------------------------------------------------


def test_safe_name_replaces_spaces() -> None:
    assert _safe_name("hello world") == "hello_world"


def test_safe_name_keeps_alnum_dash_underscore() -> None:
    assert _safe_name("cbt-030_ABC") == "cbt-030_ABC"


# ---------------------------------------------------------------------------
# _build_columns
# ---------------------------------------------------------------------------


class TestBuildColumns:
    def test_meta_columns_present(self) -> None:
        cols = _build_columns(10)
        for c in [
            "mol_id",
            "molecule",
            "smiles",
            "hf_exp_kJmol",
            "n_conformers_used",
            "status",
            "mopac_status",
            "total_time_s",
            "total_time_min",
        ]:
            assert c in cols

    def test_method_columns_all_three_methods(self) -> None:
        cols = _build_columns(10)
        for m in ["am1", "pm3", "pm7"]:
            assert f"hof_{m}_conf01_kJmol" in cols
            assert f"hof_{m}_conf10_kJmol" in cols
            assert f"absdev_{m}_conf01_kJmol" in cols
            assert f"absdev_{m}_conf10_kJmol" in cols
            assert f"best_conf_index_{m}" in cols
            assert f"best_hof_{m}_kJmol" in cols

    def test_column_count_for_10_conformers(self) -> None:
        cols = _build_columns(10)
        # 9 meta + (10 hof + 10 absdev + 2) × 3 methods = 9 + 66 = 75
        assert len(cols) == 75

    def test_column_count_for_5_conformers(self) -> None:
        cols = _build_columns(5)
        # 9 meta + (5 + 5 + 2) × 3 = 9 + 36 = 45
        assert len(cols) == 45

    def test_no_duplicates(self) -> None:
        cols = _build_columns(10)
        assert len(cols) == len(set(cols))


# ---------------------------------------------------------------------------
# process_molecule_multi
# ---------------------------------------------------------------------------


class TestProcessMoleculeMulti:
    def test_no_conformers_status(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        row = process_molecule_multi(_default_mol(), config)
        assert row["status"] == "no_conformers"
        assert row["n_conformers_used"] == 0
        assert row["mopac_status"] == "not_attempted"
        # All HoF / absdev columns should be blank
        assert row["hof_pm7_conf01_kJmol"] == ""
        assert row["best_conf_index_pm7"] == ""

    def test_hof_columns_populated(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        config = _make_config(tmp_path)
        run_dir = config.work_root / "test_001"
        _make_conformer_files(run_dir, 3)

        call_seq = [-91.0, -92.0, -93.0] * 3  # 9 calls total (3 methods × 3 confs)
        call_iter = iter(call_seq)

        monkeypatch.setattr(
            "grimperium_mini.multi_conformer.run_mopac",
            lambda *a, **kw: _mopac_ok(next(call_iter)),
        )

        row = process_molecule_multi(_default_mol(), config, max_conformers=10)

        assert row["status"] == "success"
        assert row["n_conformers_used"] == 3
        assert row["hof_am1_conf01_kJmol"] != ""
        assert row["hof_am1_conf03_kJmol"] != ""
        # Columns beyond n_confs must be blank
        assert row["hof_am1_conf04_kJmol"] == ""
        assert row["hof_pm7_conf10_kJmol"] == ""

    def test_best_conformer_selected_by_min_absdev(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """best_conf_index should point to the conformer with HoF closest to exp."""
        config = _make_config(tmp_path)
        run_dir = config.work_root / "test_001"
        _make_conformer_files(run_dir, 3)

        # exp = -100.0; PM7: [-105, -100, -80] → conf 2 is closest
        per_method: dict[str, list[float]] = {
            "AM1": [-110.0, -105.0, -90.0],
            "PM3": [-120.0, -115.0, -110.0],
            "PM7": [-105.0, -100.0, -80.0],
        }
        counters: dict[str, int] = {"AM1": 0, "PM3": 0, "PM7": 0}

        def mock_run_mopac(
            xyz_file: Path, work_dir: Path, method: str, cfg: object, **kw: object
        ) -> MopacResult:
            idx = counters[method]
            counters[method] += 1
            return _mopac_ok(per_method[method][idx])

        monkeypatch.setattr("grimperium_mini.multi_conformer.run_mopac", mock_run_mopac)

        row = process_molecule_multi(_default_mol(hf_exp=-100.0), config)

        assert row["best_conf_index_pm7"] == 2
        assert row["best_hof_pm7_kJmol"] == -100.0
        # PM3: all far; closest is conf 3 (-110, dev=10) or conf 1 (-120, dev=20)?
        # conf 3 = -110, dev=10 → conf 3 wins
        assert row["best_conf_index_pm3"] == 3

    def test_no_exp_data_yields_blank_best_and_absdev(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        config = _make_config(tmp_path)
        run_dir = config.work_root / "test_001"
        _make_conformer_files(run_dir, 2)

        monkeypatch.setattr(
            "grimperium_mini.multi_conformer.run_mopac",
            lambda *a, **kw: _mopac_ok(-95.0),
        )

        row = process_molecule_multi(_default_mol(hf_exp=None), config)

        # HoFs still recorded
        assert row["hof_pm7_conf01_kJmol"] != ""
        # Deviations and best columns blank (no experimental reference)
        assert row["absdev_pm7_conf01_kJmol"] == ""
        assert row["best_conf_index_pm7"] == ""
        assert row["best_hof_pm7_kJmol"] == ""

    def test_failed_mopac_excluded_from_best(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A conformer whose MOPAC call fails should not affect best selection."""
        config = _make_config(tmp_path)
        run_dir = config.work_root / "test_001"
        _make_conformer_files(run_dir, 3)

        # Every method: conf1 OK (-105), conf2 FAIL, conf3 OK (-80)
        counters: dict[str, int] = {"AM1": 0, "PM3": 0, "PM7": 0}
        hofs = [-105.0, None, -80.0]

        def mock_run_mopac(
            xyz_file: Path, work_dir: Path, method: str, cfg: object, **kw: object
        ) -> MopacResult:
            idx = counters[method]
            counters[method] += 1
            h = hofs[idx]
            return _mopac_ok(h) if h is not None else _mopac_fail()

        monkeypatch.setattr("grimperium_mini.multi_conformer.run_mopac", mock_run_mopac)

        row = process_molecule_multi(_default_mol(hf_exp=-100.0), config)

        # mopac_status must be partial (not all succeeded)
        assert row["mopac_status"] == "partial"
        # conf2 hof should be blank
        assert row["hof_pm7_conf02_kJmol"] == ""
        # best must be among {conf1, conf3} — conf1 (-105) is closer to -100
        assert row["best_conf_index_pm7"] == 1

    def test_all_mopac_fail_gives_failed_status(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        config = _make_config(tmp_path)
        run_dir = config.work_root / "test_001"
        _make_conformer_files(run_dir, 2)

        monkeypatch.setattr(
            "grimperium_mini.multi_conformer.run_mopac",
            lambda *a, **kw: _mopac_fail(),
        )

        row = process_molecule_multi(_default_mol(), config)
        assert row["status"] == "failed"
        assert row["mopac_status"] == "failed"
        assert row["best_conf_index_pm7"] == ""


# ---------------------------------------------------------------------------
# run_multi_conformer_pipeline (integration-level, MOPAC mocked)
# ---------------------------------------------------------------------------


@pytest.fixture()
def data_xlsx() -> Path:
    return Path("data/grimperium_mini_pipeline_tcc.xlsx")


@pytest.mark.skipif(
    not Path("data/grimperium_mini_pipeline_tcc.xlsx").exists(),
    reason="xlsx fixture not available",
)
def test_pipeline_writes_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, data_xlsx: Path
) -> None:
    config = _make_config(tmp_path)

    monkeypatch.setattr(
        "grimperium_mini.multi_conformer.run_mopac",
        lambda *a, **kw: _mopac_ok(-90.0),
    )

    rows = run_multi_conformer_pipeline(data_xlsx, config, limit=2, max_conformers=10)

    csv_path = config.results_dir / "grimperium_mini_multiconf_summary.csv"
    assert csv_path.exists(), "CSV output must be created"
    assert len(rows) == 2

    # All rows have 'no_conformers' status because tmp_path/runs is empty
    for row in rows:
        assert row["status"] == "no_conformers"


@pytest.mark.skipif(
    not Path("data/grimperium_mini_pipeline_tcc.xlsx").exists(),
    reason="xlsx fixture not available",
)
def test_pipeline_progress_callbacks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, data_xlsx: Path
) -> None:
    config = _make_config(tmp_path)
    monkeypatch.setattr(
        "grimperium_mini.multi_conformer.run_mopac",
        lambda *a, **kw: _mopac_ok(-90.0),
    )

    events: list[str] = []
    run_multi_conformer_pipeline(
        data_xlsx,
        config,
        limit=2,
        on_progress=lambda event, _data: events.append(event),
    )

    assert events.count("start") == 2
    assert events.count("done") == 2


# ---------------------------------------------------------------------------
# run_multi_conformer_pipeline — resume / skip behaviour
# ---------------------------------------------------------------------------


def test_pipeline_resume_skips_successful(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Molecules with status=success in an existing CSV must not be reprocessed."""
    import csv as _csv

    import grimperium_mini.multi_conformer as mc

    config = _make_config(tmp_path)
    output_path = config.results_dir / "grimperium_mini_multiconf_summary.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Pre-populate CSV: mol_a already succeeded.
    cols = mc._build_columns(10)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = _csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        existing: dict[str, object] = dict.fromkeys(cols, "")
        existing["mol_id"] = "mol_a"
        existing["status"] = "success"
        existing["mopac_status"] = "success"
        existing["n_conformers_used"] = "3"
        writer.writerow(existing)

    mol_a = {
        "mol_id": "mol_a",
        "molecule": "A",
        "smiles": "C",
        "charge": 0,
        "multiplicity": 1,
        "hf_exp_kJmol": None,
    }
    mol_b = {
        "mol_id": "mol_b",
        "molecule": "B",
        "smiles": "CC",
        "charge": 0,
        "multiplicity": 1,
        "hf_exp_kJmol": None,
    }
    monkeypatch.setattr(mc, "load_molecules", lambda *a, **kw: [mol_a, mol_b])

    processed: list[str] = []

    def tracking_process(
        mol: dict[str, object], cfg: object, **kw: object
    ) -> dict[str, object]:
        processed.append(str(mol["mol_id"]))
        blank: dict[str, object] = dict.fromkeys(cols, "")
        blank.update(
            {
                "mol_id": mol["mol_id"],
                "status": "no_conformers",
                "n_conformers_used": 0,
                "mopac_status": "not_attempted",
                "total_time_s": 0.0,
                "total_time_min": 0.0,
            }
        )
        return blank

    monkeypatch.setattr(mc, "process_molecule_multi", tracking_process)

    rows = mc.run_multi_conformer_pipeline(
        Path("dummy.xlsx"), config, max_conformers=10
    )

    # mol_a must not have been processed again
    assert "mol_a" not in processed
    assert "mol_b" in processed

    by_id = {str(r["mol_id"]): r for r in rows}
    assert by_id["mol_a"]["status"] == "success", "existing row must be preserved"
    assert by_id["mol_b"]["status"] == "no_conformers"


def test_pipeline_resume_progress_events_for_skipped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Skipped molecules emit a 'done' event (no 'start') with status=skipped."""
    import csv as _csv

    import grimperium_mini.multi_conformer as mc

    config = _make_config(tmp_path)
    output_path = config.results_dir / "grimperium_mini_multiconf_summary.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cols = mc._build_columns(10)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = _csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        row: dict[str, object] = dict.fromkeys(cols, "")
        row["mol_id"] = "mol_a"
        row["status"] = "success"
        writer.writerow(row)

    monkeypatch.setattr(
        mc,
        "load_molecules",
        lambda *a, **kw: [
            {
                "mol_id": "mol_a",
                "molecule": "",
                "smiles": "",
                "charge": 0,
                "multiplicity": 1,
                "hf_exp_kJmol": None,
            }
        ],
    )
    monkeypatch.setattr(mc, "process_molecule_multi", lambda *a, **kw: {})

    events: list[tuple[str, str]] = []
    mc.run_multi_conformer_pipeline(
        Path("dummy.xlsx"),
        config,
        on_progress=lambda ev, data: events.append(
            (ev, str(data.get("mol_id", "")) if isinstance(data, dict) else "")
        ),
    )

    starts = [e for e in events if e[0] == "start"]
    dones = [e for e in events if e[0] == "done"]
    assert starts == [], "skipped molecule must not fire a 'start' event"
    assert len(dones) == 1
    assert dones[0][1] == "mol_a"
