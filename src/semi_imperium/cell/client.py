"""HTTP calls from the cell to the central (``grimperium.webapi.cells``).

Every call starts here: the central never connects to a cell. A call never
raises for a network problem or an error status; it returns a ``Reply`` and
the caller decides (``outbox.classify``). The token goes only in the
``Authorization`` header of each request and is never logged or echoed.
"""

from __future__ import annotations

import logging
from typing import Any

import httpx

from semi_imperium.cell.outbox import Reply

log = logging.getLogger(__name__)

DEFAULT_TIMEOUT_S = 30.0


class CentralClient:
    def __init__(
        self,
        central_url: str,
        token: str,
        *,
        http: httpx.Client | None = None,
        timeout_s: float = DEFAULT_TIMEOUT_S,
    ) -> None:
        self._base = central_url.rstrip("/")
        self._headers = {"Authorization": f"Bearer {token}"}
        self._http = http or httpx.Client(timeout=timeout_s)

    def __repr__(self) -> str:
        return f"CentralClient({self._base!r})"

    def close(self) -> None:
        self._http.close()

    def _post(self, path: str, body: dict[str, Any]) -> Reply:
        try:
            response = self._http.post(
                f"{self._base}{path}", json=body, headers=self._headers
            )
        except httpx.HTTPError as exc:
            # The message names the URL at most, never the request headers.
            return Reply(None, {}, f"{type(exc).__name__}: {exc}")
        try:
            data = response.json()
        except ValueError:
            data = {"detail": response.text[:500]}
        if not isinstance(data, dict):
            data = {"detail": str(data)[:500]}
        return Reply(response.status_code, data)

    def register(self, body: dict[str, Any]) -> Reply:
        return self._post("/api/cell/register", body)

    def claim(self, count: int) -> Reply:
        return self._post("/api/cell/claim", {"count": count})

    def heartbeat(self, attempt_ids: list[str]) -> Reply:
        return self._post("/api/cell/heartbeat", {"attempt_ids": attempt_ids})

    def result(self, request: dict[str, Any]) -> Reply:
        return self._post("/api/cell/result", request)
