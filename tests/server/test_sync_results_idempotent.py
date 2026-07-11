from pathlib import Path

import pandas as pd
import pytest
from httpx import ASGITransport, AsyncClient

from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.server.app import create_app
from grimperium.server.config import ServerConfig


def _write_csv(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            }
        ]
    ).to_csv(path, index=False)


async def _client(csv_path: Path) -> AsyncClient:
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    return AsyncClient(transport=ASGITransport(app=app), base_url="http://test")


async def _claim(client: AsyncClient) -> str:
    await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
    await client.post("/dispatch/start")
    response = await client.post("/claim", json={"worker_id": "w1"})
    mol_id = response.json()["mol_id"]
    assert mol_id is not None
    return str(mol_id)


@pytest.mark.anyio
async def test_sync_results_duplicate_result_id_is_not_applied_twice(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "result-1",
                    "mol_id": mol_id,
                    "success": False,
                    "result_update": None,
                    "error": "MOPAC timeout",
                    "completed_at": "2026-04-21T10:00:00Z",
                }
            ],
        }

        first = await client.post("/sync_results", json=payload)
        second = await client.post("/sync_results", json=payload)

    assert first.status_code == 200
    assert first.json()["accepted"] == 1
    assert second.status_code == 200
    assert second.json()["duplicate"] is True
    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1


@pytest.mark.anyio
async def test_sync_results_conflicting_result_id_returns_409(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id = await _claim(client)
        base_result = {
            "result_id": "result-1",
            "mol_id": mol_id,
            "success": False,
            "result_update": None,
            "error": "MOPAC timeout",
            "completed_at": "2026-04-21T10:00:00Z",
        }
        first = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [base_result]}
        )
        conflicting = dict(base_result)
        conflicting["success"] = True
        conflicting["result_update"] = {"H298_pm7": -55.0}
        conflicting["error"] = None
        second = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [conflicting]}
        )

    assert first.status_code == 200
    assert second.status_code == 409
