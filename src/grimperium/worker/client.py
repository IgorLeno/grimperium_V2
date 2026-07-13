"""HTTP client for the Grimperium distributed processing server."""

import logging
import socket
import time
from dataclasses import dataclass
from typing import Any

import httpx

LOG = logging.getLogger(__name__)


class ServerError(Exception):
    """Raised when the server returns an unexpected or error response."""


@dataclass
class WorkerClientConfig:
    """Configuration for WorkerClient."""

    server_url: str
    worker_id: str
    api_token: str = ""
    timeout_s: float = 30.0
    max_retries: int = 3
    retry_delay_s: float = 1.0


class WorkerClient:
    """HTTP client wrapping all Grimperium server API calls.

    Retries on 5xx up to ``config.max_retries`` times.
    Raises ``ServerError`` on 4xx or exhausted retries.
    Authentication via ``X-Token`` header when ``api_token`` is set.
    """

    def __init__(
        self,
        config: WorkerClientConfig,
        _transport: httpx.BaseTransport | None = None,
    ) -> None:
        self._config = config
        headers: dict[str, str] = {}
        if config.api_token:
            headers["X-Token"] = config.api_token
        self._http = httpx.Client(
            base_url=config.server_url,
            headers=headers,
            timeout=config.timeout_s,
            transport=_transport,
        )

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _post(self, path: str, body: dict[str, Any]) -> dict[str, Any]:
        """POST with automatic retry on 5xx. Raises ServerError on 4xx."""
        last_exc: Exception | None = None
        for attempt in range(self._config.max_retries + 1):
            try:
                r = self._http.post(path, json=body)
                if r.status_code == 200:
                    return dict(r.json())
                if r.status_code >= 500 and attempt < self._config.max_retries:
                    LOG.warning(
                        "POST %s returned %d — retrying (%d/%d)",
                        path,
                        r.status_code,
                        attempt + 1,
                        self._config.max_retries,
                    )
                    time.sleep(self._config.retry_delay_s)
                    continue
                raise ServerError(f"POST {path} → {r.status_code}: {r.text}")
            except httpx.TransportError as exc:
                last_exc = exc
                if attempt < self._config.max_retries:
                    LOG.warning("POST %s transport error — retrying: %s", path, exc)
                    time.sleep(self._config.retry_delay_s)
        raise ServerError(
            f"POST {path} failed after {self._config.max_retries + 1} attempts"
        ) from last_exc

    # ── Public API ────────────────────────────────────────────────────────────

    def register(self, hostname: str | None = None) -> dict[str, Any]:
        """Register with the server.

        Returns:
            Server response dict containing worker config (crest_timeout_minutes,
            mopac_timeout_minutes, batch_size, profile_name).
        """
        hn = hostname if hostname is not None else socket.gethostname()
        return self._post(
            "/register", {"worker_id": self._config.worker_id, "hostname": hn}
        )

    def get_config(self) -> dict[str, Any]:
        """Fetch current server-side config for this worker.

        Used by WorkerRunner to reconfigure before each claim().

        Returns:
            Config dict from GET /configure/{worker_id}.
        """
        r = self._http.get(f"/configure/{self._config.worker_id}")
        if r.status_code == 200:
            return dict(r.json())
        LOG.warning("GET /configure/%s → %d", self._config.worker_id, r.status_code)
        return {}

    def claim(self) -> tuple[str, str, str | None] | None:
        data = self._post("/claim", {"worker_id": self._config.worker_id})
        mol_id = data.get("mol_id")
        if mol_id is None:
            return None
        smiles = data.get("smiles") or ""
        attempt_id = data.get("attempt_id")
        return mol_id, smiles, str(attempt_id) if attempt_id else None

    def heartbeat(self, mol_id: str) -> None:
        r = self._http.put(
            f"/heartbeat/{mol_id}",
            json={"worker_id": self._config.worker_id},
        )
        if r.status_code not in (200, 404):
            raise ServerError(f"PUT /heartbeat/{mol_id} → {r.status_code}: {r.text}")

    def report_success(self, mol_id: str, result_update: dict[str, Any]) -> None:
        """Legado: preferir sync_results. Mantido para clientes antigos/testes."""
        self._post(
            "/report/success",
            {
                "worker_id": self._config.worker_id,
                "mol_id": mol_id,
                "result_update": result_update,
            },
        )

    def report_failure(self, mol_id: str, error: str, force_skip: bool = False) -> None:
        """Legado: preferir sync_results. Mantido para clientes antigos/testes."""
        self._post(
            "/report/failure",
            {
                "worker_id": self._config.worker_id,
                "mol_id": mol_id,
                "error": error,
                "force_skip": force_skip,
            },
        )

    def sync_results(self, results: list[dict[str, Any]]) -> dict[str, Any]:
        """Enviar resultados pelo protocolo canônico idempotente.

        Returns:
            Resposta completa do servidor, incluindo ``items`` por result_id.
        """
        data = self._post(
            "/sync_results",
            {"worker_id": self._config.worker_id, "results": results},
        )
        return dict(data)
