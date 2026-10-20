"""The cell's config file and its HTTP client."""

from __future__ import annotations

import json
import secrets
from pathlib import Path
from typing import Any

import httpx
import pytest

from semi_imperium.cell.client import CentralClient
from semi_imperium.cell.config import (
    CellConfig,
    ConfigError,
    load_config,
    worker_count,
)


def _write(path: Path, data: dict[str, Any], mode: int = 0o600) -> Path:
    path.write_text(json.dumps(data))
    path.chmod(mode)
    return path


def test_a_minimal_config_works_in_its_own_folder(tmp_path: Path) -> None:
    token = secrets.token_urlsafe(32)
    path = _write(
        tmp_path / "config.json",
        {"central_url": "http://100.64.0.1:8100/", "token": token},
    )

    config = load_config(path)

    assert config.central_url == "http://100.64.0.1:8100"
    assert config.root == tmp_path and config.crest_threads == 4
    assert config.work_dir == tmp_path / "work"
    assert token not in repr(config)


def test_a_config_others_can_read_is_refused(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "config.json",
        {"central_url": "http://c:8100", "token": "t"},
        mode=0o644,
    )
    with pytest.raises(ConfigError, match="chmod 600"):
        load_config(path)


@pytest.mark.parametrize(
    ("data", "message"),
    [
        ({"token": "t"}, "central_url"),
        ({"central_url": "http://c"}, "token"),
        ({"central_url": "c:8100", "token": "t"}, "http"),
        ({"central_url": "http://c", "token": "t", "crest_threads": 0}, ">= 1"),
        ({"central_url": "http://c", "token": "t", "workers": "x"}, "inválido"),
        ({"central_url": "http://c", "token": "t", "wrkers": 2}, "desconhecidas"),
        (
            {"central_url": "http://c", "token": "t", "root": "/mnt/c/celula"},
            "Windows",
        ),
    ],
)
def test_bad_configs_say_what_is_wrong(
    tmp_path: Path, data: dict[str, Any], message: str
) -> None:
    with pytest.raises(ConfigError, match=message):
        load_config(_write(tmp_path / "config.json", data))


def test_a_missing_config_is_a_config_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="não encontrada"):
        load_config(tmp_path / "config.json")


def test_workers_use_every_core_unless_set() -> None:
    config = CellConfig(central_url="http://c", token="t", crest_threads=4)
    assert worker_count(config, 16) == 4
    assert worker_count(config, 6) == 1
    assert worker_count(config, 2) == 1  # never zero
    fixed = CellConfig(central_url="http://c", token="t", workers=3)
    assert worker_count(fixed, 16) == 3


def test_the_client_sends_the_token_only_in_the_header() -> None:
    token = secrets.token_urlsafe(32)
    seen: list[httpx.Request] = []

    def answer(request: httpx.Request) -> httpx.Response:
        seen.append(request)
        return httpx.Response(200, json={"renewed": 0, "not_held": []})

    client = CentralClient(
        "http://c:8100", token, http=httpx.Client(transport=httpx.MockTransport(answer))
    )

    reply = client.heartbeat([])

    assert reply.status == 200 and reply.body["renewed"] == 0
    (request,) = seen
    assert request.headers["Authorization"] == f"Bearer {token}"
    assert token not in str(request.url) and token.encode() not in request.content
    assert token not in repr(client)


def test_no_answer_is_a_reply_without_status() -> None:
    def refuse(request: httpx.Request) -> httpx.Response:
        raise httpx.ConnectError("connection refused", request=request)

    client = CentralClient(
        "http://c:8100", "t", http=httpx.Client(transport=httpx.MockTransport(refuse))
    )

    reply = client.claim(2)

    assert reply.status is None and "ConnectError" in reply.detail


def test_a_body_that_is_not_json_still_gives_a_detail() -> None:
    client = CentralClient(
        "http://c:8100",
        "t",
        http=httpx.Client(
            transport=httpx.MockTransport(
                lambda _r: httpx.Response(502, text="Bad Gateway")
            )
        ),
    )

    reply = client.result({"attempt_id": "a", "failure": "x"})

    assert reply.status == 502 and reply.detail == "Bad Gateway"
