"""Tests for worker/client.py — WorkerClient HTTP calls."""

import json as _json
from collections.abc import Callable
from typing import Any

import httpx
import pytest

from grimperium.worker.client import ServerError, WorkerClient, WorkerClientConfig

# ── Minimal mock transport ─────────────────────────────────────────────────────

HandlerFn = Callable[[httpx.Request], httpx.Response]


class _MockTransport(httpx.BaseTransport):
    """Records requests and returns pre-programmed responses."""

    def __init__(self) -> None:
        self._routes: list[tuple[str, str, httpx.Response | HandlerFn]] = []
        self.requests: list[httpx.Request] = []

    def add(
        self,
        method: str,
        path: str,
        response: httpx.Response | HandlerFn,
    ) -> None:
        self._routes.append((method, path, response))

    def handle_request(self, request: httpx.Request) -> httpx.Response:
        self.requests.append(request)
        for method, path, resp in self._routes:
            if request.method == method and request.url.path == path:
                if callable(resp):
                    return resp(request)
                return resp
        return httpx.Response(
            500,
            json={"detail": f"no mock: {request.method} {request.url.path}"},
        )


def _json_response(
    status: int = 200, body: dict[str, Any] | None = None
) -> httpx.Response:
    return httpx.Response(status, json=body or {})


def _make_client(
    transport: _MockTransport,
    token: str = "",
    max_retries: int = 0,
    retry_delay_s: float = 0.0,
) -> WorkerClient:
    cfg = WorkerClientConfig(
        server_url="http://test-server",
        worker_id="w1",
        api_token=token,
        max_retries=max_retries,
        retry_delay_s=retry_delay_s,
    )
    return WorkerClient(cfg, _transport=transport)


# ── /register ─────────────────────────────────────────────────────────────────


class TestRegister:
    def test_register_succeeds(self) -> None:
        t = _MockTransport()
        t.add("POST", "/register", _json_response(200, {"status": "registered"}))
        client = _make_client(t)
        client.register("lab-01")  # must not raise

    def test_register_sends_worker_id(self) -> None:
        t = _MockTransport()
        t.add("POST", "/register", _json_response(200, {"status": "registered"}))
        client = _make_client(t)
        client.register("lab-01")
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"
        assert body["hostname"] == "lab-01"

    def test_register_uses_socket_hostname_when_not_given(self) -> None:
        import socket

        t = _MockTransport()
        t.add("POST", "/register", _json_response(200, {"status": "registered"}))
        client = _make_client(t)
        client.register()
        body = _json.loads(t.requests[-1].content)
        assert body["hostname"] == socket.gethostname()

    def test_register_raises_on_server_error(self) -> None:
        t = _MockTransport()
        t.add("POST", "/register", _json_response(500, {"detail": "boom"}))
        client = _make_client(t, max_retries=0)
        with pytest.raises(ServerError):
            client.register("lab-01")


# ── /claim ────────────────────────────────────────────────────────────────────


class TestClaim:
    def test_claim_returns_mol_id_and_smiles(self) -> None:
        t = _MockTransport()
        t.add(
            "POST",
            "/claim",
            _json_response(
                200, {"mol_id": "m1", "smiles": "CCO", "attempt_id": "att-1"}
            ),
        )
        client = _make_client(t)
        result = client.claim()
        assert result == ("m1", "CCO", "att-1")

    def test_claim_returns_none_when_queue_empty(self) -> None:
        t = _MockTransport()
        t.add("POST", "/claim", _json_response(200, {"mol_id": None, "smiles": None}))
        client = _make_client(t)
        result = client.claim()
        assert result is None

    def test_claim_sends_worker_id(self) -> None:
        t = _MockTransport()
        t.add("POST", "/claim", _json_response(200, {"mol_id": None, "smiles": None}))
        client = _make_client(t)
        client.claim()
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"


# ── /heartbeat ────────────────────────────────────────────────────────────────


class TestHeartbeat:
    def test_heartbeat_succeeds(self) -> None:
        t = _MockTransport()
        t.add("PUT", "/heartbeat/m1", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.heartbeat("m1", attempt_id="att-1")  # must not raise

    def test_heartbeat_silently_tolerates_404(self) -> None:
        t = _MockTransport()
        t.add("PUT", "/heartbeat/ghost", _json_response(404, {"detail": "not running"}))
        client = _make_client(t)
        client.heartbeat("ghost", attempt_id="att-1")  # must not raise

    def test_heartbeat_raises_on_5xx(self) -> None:
        t = _MockTransport()
        t.add("PUT", "/heartbeat/m1", _json_response(500, {"detail": "error"}))
        client = _make_client(t)
        with pytest.raises(ServerError):
            client.heartbeat("m1", attempt_id="att-1")

    def test_heartbeat_sends_worker_id(self) -> None:
        t = _MockTransport()
        t.add("PUT", "/heartbeat/m1", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.heartbeat("m1")
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"

    def test_heartbeat_sends_attempt_id_when_available(self) -> None:
        t = _MockTransport()
        t.add("PUT", "/heartbeat/m1", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.heartbeat("m1", attempt_id="att-42")
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"
        assert body["attempt_id"] == "att-42"


# ── /report/success ───────────────────────────────────────────────────────────


class TestReportSuccess:
    def test_report_success_succeeds(self) -> None:
        t = _MockTransport()
        t.add("POST", "/report/success", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.report_success("m1", {"H298_pm7": -42.0})

    def test_report_success_sends_payload(self) -> None:
        t = _MockTransport()
        t.add("POST", "/report/success", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.report_success("m1", {"H298_pm7": -42.0})
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"
        assert body["mol_id"] == "m1"
        assert body["result_update"]["H298_pm7"] == -42.0

    def test_report_success_raises_on_404(self) -> None:
        t = _MockTransport()
        t.add("POST", "/report/success", _json_response(404, {"detail": "not running"}))
        client = _make_client(t)
        with pytest.raises(ServerError):
            client.report_success("m1", {})


# ── /report/failure ───────────────────────────────────────────────────────────


class TestReportFailure:
    def test_report_failure_succeeds(self) -> None:
        t = _MockTransport()
        t.add("POST", "/report/failure", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.report_failure("m1", "CREST timeout")

    def test_report_failure_sends_payload(self) -> None:
        t = _MockTransport()
        t.add("POST", "/report/failure", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.report_failure("m1", "CREST timeout", force_skip=True)
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"
        assert body["mol_id"] == "m1"
        assert body["error"] == "CREST timeout"
        assert body["force_skip"] is True

    def test_report_failure_force_skip_defaults_false(self) -> None:
        t = _MockTransport()
        t.add("POST", "/report/failure", _json_response(200, {"status": "ok"}))
        client = _make_client(t)
        client.report_failure("m1", "err")
        body = _json.loads(t.requests[-1].content)
        assert body["force_skip"] is False


# ── /sync_results ─────────────────────────────────────────────────────────────


class TestSyncResults:
    def test_sync_results_returns_response_dict(self) -> None:
        t = _MockTransport()
        t.add(
            "POST",
            "/sync_results",
            _json_response(
                200,
                {
                    "accepted": 3,
                    "rejected": 1,
                    "items": [
                        {
                            "result_id": "a",
                            "mol_id": "m1",
                            "status": "applied",
                        }
                    ],
                },
            ),
        )
        client = _make_client(t)
        data = client.sync_results([])
        assert data["accepted"] == 3
        assert data["rejected"] == 1
        assert data["items"][0]["status"] == "applied"

    def test_sync_results_sends_worker_id(self) -> None:
        t = _MockTransport()
        t.add(
            "POST", "/sync_results", _json_response(200, {"accepted": 0, "rejected": 0})
        )
        client = _make_client(t)
        client.sync_results(
            [
                {
                    "mol_id": "m1",
                    "success": True,
                    "result_update": {},
                    "error": None,
                    "completed_at": "2024-01-01T00:00:00Z",
                }
            ]
        )
        body = _json.loads(t.requests[-1].content)
        assert body["worker_id"] == "w1"
        assert len(body["results"]) == 1


# ── Auth header ───────────────────────────────────────────────────────────────


class TestAuthHeader:
    def test_token_sent_in_header(self) -> None:
        t = _MockTransport()
        t.add("POST", "/register", _json_response(200, {"status": "registered"}))
        client = _make_client(t, token="secret-token")  # noqa: S106
        client.register("lab-01")
        assert t.requests[-1].headers.get("x-token") == "secret-token"

    def test_no_token_header_when_empty(self) -> None:
        t = _MockTransport()
        t.add("POST", "/register", _json_response(200, {"status": "registered"}))
        client = _make_client(t, token="")
        client.register("lab-01")
        assert "x-token" not in t.requests[-1].headers


# ── Retry logic ───────────────────────────────────────────────────────────────


class TestRetryLogic:
    def test_retries_on_5xx_and_succeeds(self) -> None:
        call_count = 0

        def handler(req: httpx.Request) -> httpx.Response:
            nonlocal call_count
            call_count += 1
            if call_count < 3:
                return _json_response(503, {"detail": "unavailable"})
            return _json_response(200, {"status": "registered"})

        t = _MockTransport()
        t.add("POST", "/register", handler)
        client = _make_client(t, max_retries=3, retry_delay_s=0.0)
        client.register("lab-01")
        assert call_count == 3

    def test_raises_after_all_retries_exhausted(self) -> None:
        t = _MockTransport()
        t.add("POST", "/register", _json_response(503, {"detail": "unavailable"}))
        client = _make_client(t, max_retries=2, retry_delay_s=0.0)
        with pytest.raises(ServerError):
            client.register("lab-01")

    def test_no_retry_on_4xx(self) -> None:
        call_count = 0

        def handler(req: httpx.Request) -> httpx.Response:
            nonlocal call_count
            call_count += 1
            return _json_response(404, {"detail": "not found"})

        t = _MockTransport()
        t.add("POST", "/report/success", handler)
        client = _make_client(t, max_retries=3, retry_delay_s=0.0)
        with pytest.raises(ServerError):
            client.report_success("m1", {})
        assert call_count == 1
