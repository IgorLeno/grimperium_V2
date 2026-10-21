"""Cell package: the zip the dashboard hands out, its token and its contents."""

from __future__ import annotations

import importlib.metadata
import io
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest
from fastapi.testclient import TestClient
from packaging.requirements import Requirement

from grimperium.webapi import cell_package
from grimperium.webapi.app import create_app
from grimperium.webapi.cell_package import (
    CELL_DISTRIBUTIONS,
    PackageError,
    WheelBuild,
    suggest_central_url,
)
from grimperium.webapi.config import DashboardPaths
from semi_imperium.training.cells import CellRegistry
from semi_imperium.training.db import database_path
from tests.training.conftest import row

LOCAL = {"Origin": "http://127.0.0.1:5173"}
URL = "http://100.64.0.1:8100"
PROGRAMS = {
    "crest": "3.0.2",
    "crest_sha256": "abc",
    "crest_path": "/home/central/micromamba/envs/grimperium-quimica/bin/crest",
    "environment_sha256": "def",
}


def fake_wheel(out_dir: Path) -> WheelBuild:
    wheel = out_dir / "grimperium-0.2.0-py3-none-any.whl"
    wheel.write_bytes(b"wheel bytes")
    return WheelBuild(path=wheel, commit="f" * 40, dirty=False)


class Env:
    def __init__(self, tmp_path: Path, builder=fake_wheel) -> None:  # type: ignore[no-untyped-def]
        pd.DataFrame([row("a", nheavy=1)]).to_csv(tmp_path / "cbs.csv", index=False)
        self.paths = DashboardPaths(
            semi_imperium_store=tmp_path / "store",
            exp_csv=tmp_path / "exp.csv",
            exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
            cbs_csv=tmp_path / "cbs.csv",
        )
        self.cells = CellRegistry(database_path(self.paths.semi_imperium_store))
        self.client = TestClient(
            create_app(
                self.paths,
                launcher=lambda _: None,
                wheel_builder=builder,
                central_programs=lambda: dict(PROGRAMS),
            )
        )

    def create(self, name: str = "PC do Prof. Silva") -> tuple[str, str]:
        body = self.client.post("/api/cells", json={"name": name}, headers=LOCAL)
        assert body.status_code == 200
        return body.json()["cell"]["cell_id"], body.json()["token"]

    def package(self, cell_id: str, url: str = URL) -> object:
        return self.client.post(
            f"/api/cells/{cell_id}/package", json={"central_url": url}, headers=LOCAL
        )


@pytest.fixture
def env(tmp_path: Path) -> Env:
    return Env(tmp_path)


def _zip(response: object) -> zipfile.ZipFile:
    return zipfile.ZipFile(io.BytesIO(response.content))  # type: ignore[attr-defined]


def test_package_carries_a_new_token_and_the_old_one_stops(env: Env) -> None:
    cell_id, old = env.create()
    response = env.package(cell_id, URL + "/")
    assert response.status_code == 200  # type: ignore[attr-defined]
    headers = response.headers  # type: ignore[attr-defined]
    assert headers["content-type"] == "application/zip"
    assert headers["cache-control"] == "no-store"
    folder = f"celula-pc-do-prof-silva-{cell_id.removeprefix('cell-')}"
    assert f'filename="{folder}.zip"' in headers["content-disposition"]

    archive = _zip(response)
    config = json.loads(archive.read(f"{folder}/config.json"))
    assert config["central_url"] == URL
    token = config["token"]
    assert token != old
    assert env.cells.authenticate(old) is None
    found = env.cells.authenticate(token)
    assert found is not None and found.cell_id == cell_id


def test_package_contents_and_file_modes(env: Env) -> None:
    cell_id, _ = env.create()
    archive = _zip(env.package(cell_id))
    names = {
        name.split("/", 1)[1]: info
        for name, info in ((i.filename, i) for i in archive.infolist())
    }
    assert set(names) == {
        "LEIA-ME.txt",
        "instalar.sh",
        "windows/configurar-windows.ps1",
        "config.json",
        "pacote.json",
        "pacote.env",
        "conda-lock.yml",
        "environment.yml",
        "requirements-celula.txt",
        "constraints-central.txt",
        "wheels/grimperium-0.2.0-py3-none-any.whl",
    }
    assert names["config.json"].external_attr >> 16 & 0o777 == 0o600
    assert names["instalar.sh"].external_attr >> 16 & 0o777 == 0o755
    folder = archive.namelist()[0].split("/")[0]

    manifest = json.loads(archive.read(f"{folder}/pacote.json"))
    assert manifest["cell_id"] == cell_id and manifest["commit"] == "f" * 40
    assert manifest["python"] == f"{sys.version_info.major}.{sys.version_info.minor}"
    # Install paths differ between machines; only versions and hashes compare.
    assert manifest["central_programs"] == {
        k: v for k, v in PROGRAMS.items() if not k.endswith("_path")
    }
    repo = Path(cell_package.__file__).resolve().parents[3]
    assert (
        archive.read(f"{folder}/conda-lock.yml")
        == (repo / "conda-lock.yml").read_bytes()
    )

    requirements = archive.read(f"{folder}/requirements-celula.txt").decode().split()
    assert [r.split("==")[0] for r in requirements] == list(CELL_DISTRIBUTIONS)
    for line in requirements:
        name, version = line.split("==")
        assert importlib.metadata.version(name) == version
    constraints = archive.read(f"{folder}/constraints-central.txt").decode()
    assert "grimperium==" not in constraints.lower()

    # The token lives only in config.json.
    token = json.loads(archive.read(f"{folder}/config.json"))["token"]
    for name in archive.namelist():
        if not name.endswith("config.json"):
            assert token.encode() not in archive.read(name)


def test_failed_build_keeps_the_installed_cells_token(tmp_path: Path) -> None:
    def broken(_out_dir: Path) -> WheelBuild:
        raise PackageError("poetry não encontrado")

    env = Env(tmp_path, builder=broken)
    cell_id, token = env.create()
    response = env.package(cell_id)
    assert response.status_code == 503  # type: ignore[attr-defined]
    assert "poetry" in response.json()["detail"]  # type: ignore[attr-defined]
    assert env.cells.authenticate(token) is not None


def test_revoked_unknown_and_invalid_requests(env: Env) -> None:
    cell_id, _ = env.create()
    assert env.package("nope").status_code == 404  # type: ignore[attr-defined]
    for bad in ("100.64.0.1:8100", "ftp://x", "http://a b"):
        assert env.package(cell_id, bad).status_code == 422  # type: ignore[attr-defined]
    foreign = env.client.post(
        f"/api/cells/{cell_id}/package",
        json={"central_url": URL},
        headers={"Origin": "http://evil.example"},
    )
    assert foreign.status_code == 403
    env.client.post(f"/api/cells/{cell_id}/revoke", headers=LOCAL)
    assert env.package(cell_id).status_code == 409  # type: ignore[attr-defined]


def test_admission_by_hand(env: Env) -> None:
    cell_id, _ = env.create()
    path = f"/api/cells/{cell_id}/admission"
    response = env.client.post(path, json={"admission": "passed"}, headers=LOCAL)
    assert response.status_code == 200 and response.json()["admission"] == "passed"
    assert env.cells.get(cell_id).admission == "passed"
    bad = env.client.post(path, json={"admission": "yes"}, headers=LOCAL)
    assert bad.status_code == 422
    foreign = env.client.post(
        path, json={"admission": "passed"}, headers={"Origin": "http://evil.example"}
    )
    assert foreign.status_code == 403
    unknown = env.client.post(
        "/api/cells/nope/admission", json={"admission": "passed"}, headers=LOCAL
    )
    assert unknown.status_code == 404


def test_suggested_central_url(env: Env, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GRIMPERIUM_CELL_URL", "http://100.64.0.9:8100/")
    response = env.client.get("/api/cells/package-defaults")
    assert response.json() == {"central_url": "http://100.64.0.9:8100"}
    monkeypatch.setattr(cell_package.shutil, "which", lambda _name: None)
    assert suggest_central_url({}) == ""


# --- what the cell imports must be what the installer installs ---------------


def _requirement_closure(names: tuple[str, ...]) -> set[str]:
    """Distributions pip installs for ``names`` (no extras), normalized."""
    seen: set[str] = set()
    pending = list(names)
    while pending:
        name = re.sub(r"[-_.]+", "-", pending.pop()).lower()
        if name in seen:
            continue
        seen.add(name)
        for spec in importlib.metadata.requires(name) or []:
            requirement = Requirement(spec)
            if requirement.marker is None or requirement.marker.evaluate({"extra": ""}):
                pending.append(requirement.name)
    return seen


#: Imported only inside ``try/except ImportError`` by libraries the cell uses.
OPTIONAL_IMPORTS = {"click", "charset-normalizer"}


def test_cell_imports_only_what_the_package_installs() -> None:
    script = (
        "import sys, importlib.metadata as md\n"
        "import semi_imperium.cell.__main__, semi_imperium.cell.loop\n"
        "import semi_imperium.training.worker, semi_imperium.training.__main__\n"
        "owners = md.packages_distributions()\n"
        "tops = {m.split('.')[0] for m in sys.modules} - set(sys.stdlib_module_names)\n"
        "for top in sorted(tops):\n"
        "    for dist in owners.get(top, []):\n"
        "        print(dist)\n"
    )
    output = subprocess.run(  # noqa: S603
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    ).stdout
    imported = {re.sub(r"[-_.]+", "-", d).lower() for d in output.split()}
    # Namespace packages (sphinxcontrib) get pulled in by a .pth, not an import.
    imported = {d for d in imported if not d.startswith("sphinxcontrib-")}
    missing = imported - _requirement_closure(CELL_DISTRIBUTIONS) - OPTIONAL_IMPORTS
    assert not missing, f"add to CELL_DISTRIBUTIONS: {sorted(missing)}"


def test_installer_scripts_parse() -> None:
    templates = cell_package.TEMPLATES
    subprocess.run(  # noqa: S603
        ["bash", "-n", str(templates / "instalar.sh")], check=True  # noqa: S607
    )
    windows = (templates / "configurar-windows.ps1").read_text(encoding="utf-8")
    for placeholder in ("@DISTRO@", "@LINUX_USER@", "@INICIAR@"):
        assert placeholder in windows
