"""``celula start | status | update``."""

from __future__ import annotations

import fcntl
import json
import logging
import secrets
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from semi_imperium.cell import __main__ as cli
from semi_imperium.cell.loop import STATE_FILENAME
from semi_imperium.cell.outbox import Outbox, Reply


@pytest.fixture(autouse=True)
def _detach_log_files() -> Iterator[None]:
    """``start`` attaches a file handler; never leak it into other tests."""
    yield
    logger = logging.getLogger("semi_imperium.cell")
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()


@pytest.fixture
def home(tmp_path: Path) -> Path:
    return tmp_path / "grimperium-celula"


@pytest.fixture
def token() -> str:
    return secrets.token_urlsafe(32)


@pytest.fixture
def config(home: Path, token: str) -> Path:
    home.mkdir()
    path = home / "config.json"
    path.write_text(
        json.dumps({"central_url": "http://100.64.0.1:8100", "token": token})
    )
    path.chmod(0o600)
    return path


@pytest.fixture
def running(home: Path, config: Path) -> Iterator[None]:
    """Another ``celula start`` holding the lock of this folder."""
    with (home / cli.LOCK_FILENAME).open("w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


class FakeLoop:
    instances: list[FakeLoop] = []

    def __init__(self, config: Any, client: Any, **_kwargs: Any) -> None:
        self.config, self.client = config, client
        FakeLoop.instances.append(self)

    def run(self) -> int:
        logging.getLogger("semi_imperium.cell.loop").info("rodando")
        return 0


def test_start_runs_the_loop_with_the_config_and_logs_to_a_file(
    config: Path, home: Path, token: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    FakeLoop.instances = []
    monkeypatch.setattr(cli, "CellLoop", FakeLoop)

    assert cli.main(["--config", str(config), "start"]) == 0

    (loop,) = FakeLoop.instances
    assert loop.config.root == home and loop.config.token == token
    log = (home / "celula.log").read_text()
    assert "rodando" in log and token not in log


def test_start_refuses_a_config_others_can_read(
    config: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    config.chmod(0o644)

    assert cli.main(["--config", str(config), "start"]) == 2
    assert "chmod 600" in capsys.readouterr().err


@pytest.mark.usefixtures("running")
def test_a_second_start_in_the_same_folder_is_refused(
    config: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert cli.main(["--config", str(config), "start"]) == 3
    assert "já está rodando" in capsys.readouterr().err


def test_status_before_the_first_start(
    config: Path, capsys: pytest.CaptureFixture[str], token: str
) -> None:
    assert cli.main(["--config", str(config), "status"]) == 0

    out = capsys.readouterr().out
    assert "nunca iniciada" in out and "rodando: não" in out
    assert "token: presente" in out and token not in out


@pytest.mark.usefixtures("running")
def test_status_shows_the_state_the_queue_and_the_dead_letter(
    config: Path, home: Path, capsys: pytest.CaptureFixture[str], token: str
) -> None:
    (home / STATE_FILENAME).write_text(
        json.dumps(
            {
                "phase": "calculando",
                "name": "Lab do Prof. X",
                "cell_id": "cell-1",
                "admission": "passed",
                "workers": 2,
                "crest_threads": 4,
                "busy": [
                    {
                        "attempt_id": "a" * 32,
                        "mol_id": "cbs_00187",
                        "smiles": "[H][C]=C",
                        "started_at": "2026-10-20T12:00:00+00:00",
                    }
                ],
                "last_contact_at": "2026-10-20T12:01:00+00:00",
                "last_error": None,
            }
        )
    )
    box = Outbox(home)
    box.put("b" * 32, "cbs_00001", {"attempt_id": "b" * 32, "failure": "x"})
    box.put("c" * 32, "cbs_00002", {"attempt_id": "c" * 32, "failure": "x"})
    (dead,) = [e for e in box.pending() if e.attempt_id == "c" * 32]
    box.kill(dead, Reply(422, {"detail": "Assinatura diferente"}))

    assert cli.main(["--config", str(config), "status"]) == 0

    out = capsys.readouterr().out
    assert "rodando: sim" in out and "calculando" in out
    assert "Lab do Prof. X" in out and "cbs_00187" in out
    assert "resultados a enviar: 1" in out and "dead-letter: 1" in out
    assert "cbs_00002" in out and "422" in out and "Assinatura diferente" in out
    assert token not in out


def test_update_installs_the_wheel_with_this_python(
    config: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    wheel = tmp_path / "grimperium-0.3.0-py3-none-any.whl"
    wheel.write_bytes(b"PK")
    calls: list[list[str]] = []

    def pip(args: list[str]) -> int:
        calls.append(args)
        return 0

    assert cli.main(["--config", str(config), "update", str(wheel)], pip=pip) == 0

    assert calls == [
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-deps",
            "--force-reinstall",
            str(wheel),
        ]
    ]
    assert "celula start" in capsys.readouterr().out


@pytest.mark.usefixtures("running")
def test_update_refuses_while_the_cell_runs(config: Path, tmp_path: Path) -> None:
    wheel = tmp_path / "x-1.0-py3-none-any.whl"
    wheel.write_bytes(b"PK")

    assert cli.main(["--config", str(config), "update", str(wheel)], pip=_never) == 3


@pytest.mark.parametrize("name", ["missing.whl", "package.zip"])
def test_update_needs_an_existing_wheel(
    config: Path, tmp_path: Path, name: str
) -> None:
    if name.endswith(".zip"):
        (tmp_path / name).write_bytes(b"PK")

    code = cli.main(
        ["--config", str(config), "update", str(tmp_path / name)], pip=_never
    )

    assert code == 2


def test_a_failed_install_is_reported(config: Path, tmp_path: Path) -> None:
    wheel = tmp_path / "x-1.0-py3-none-any.whl"
    wheel.write_bytes(b"PK")

    code = cli.main(
        ["--config", str(config), "update", str(wheel)], pip=lambda _args: 1
    )

    assert code == 1


def _never(_args: list[str]) -> int:
    raise AssertionError("pip must not run")
