"""The cell's configuration: ``~/grimperium-celula/config.json``.

The file holds the cell's token, so it must be readable by its owner only
(``chmod 600``); the token never appears in ``repr``, logs or the state file.
The working folder lives in the Linux home of WSL, never on a Windows drive
(``/mnt/c``): CREST writes thousands of small files and the 9p bridge to
Windows is both slow and lossy on a power cut.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

CONFIG_FILENAME = "config.json"


def default_home() -> Path:
    return Path.home() / "grimperium-celula"


class ConfigError(ValueError):
    """The configuration cannot be used; the message says why."""


@dataclass(frozen=True)
class CellConfig:
    central_url: str
    token: str = field(repr=False)
    root: Path = field(default_factory=default_home)
    crest_threads: int = 4
    #: ``None``: maximum use of the CPU, cores ÷ crest_threads.
    workers: int | None = None
    #: Folder with crest, xtb and mopac; ``None``: the fixed environment.
    bin_dir: Path | None = None
    crest_timeout_s: float = 3600.0
    mopac_timeout_s: float = 900.0
    keep_artifacts: bool = False

    def __post_init__(self) -> None:
        if not self.central_url.startswith(("http://", "https://")):
            raise ConfigError(f"central_url sem http(s)://: {self.central_url!r}")
        if not self.token.strip():
            raise ConfigError("token vazio")
        if self.crest_threads < 1:
            raise ConfigError(f"crest_threads deve ser >= 1: {self.crest_threads}")
        if self.workers is not None and self.workers < 1:
            raise ConfigError(f"workers deve ser >= 1: {self.workers}")
        if str(self.root.resolve()).startswith("/mnt/"):
            raise ConfigError(
                f"pasta de trabalho num disco do Windows ({self.root}); "
                "use uma pasta do Linux, como ~/grimperium-celula"
            )

    @property
    def store_root(self) -> Path:
        return self.root / "store"

    @property
    def work_dir(self) -> Path:
        return self.root / "work"


def cpu_cores() -> int:
    """Cores this process may use (WSL2 sees what ``.wslconfig`` grants)."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:  # not Linux
        return os.cpu_count() or 1


def worker_count(config: CellConfig, cores: int) -> int:
    """Full use of the CPU: one CREST of ``crest_threads`` per worker."""
    if config.workers is not None:
        return config.workers
    return max(1, cores // config.crest_threads)


def _optional_path(value: Any) -> Path | None:
    return None if value in (None, "") else Path(str(value)).expanduser()


def load_config(path: Path) -> CellConfig:
    """Read ``path``; refuse it when others can read the token."""
    try:
        mode = path.stat().st_mode
    except FileNotFoundError:
        raise ConfigError(f"configuração não encontrada: {path}") from None
    if mode & 0o077:
        raise ConfigError(
            f"{path} pode ser lido por outros usuários e contém o token; "
            f"rode: chmod 600 {path}"
        )
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ConfigError(f"{path} não é JSON válido: {exc}") from None
    if not isinstance(data, dict):
        raise ConfigError(f"{path} deve conter um objeto JSON")
    missing = [key for key in ("central_url", "token") if not data.get(key)]
    if missing:
        raise ConfigError(f"{path} sem: {', '.join(missing)}")
    known = {
        "central_url",
        "token",
        "root",
        "crest_threads",
        "workers",
        "bin_dir",
        "crest_timeout_s",
        "mopac_timeout_s",
        "keep_artifacts",
    }
    unknown = sorted(set(data) - known)
    if unknown:
        raise ConfigError(f"{path}: chaves desconhecidas: {', '.join(unknown)}")
    try:
        return CellConfig(
            central_url=str(data["central_url"]).rstrip("/"),
            token=str(data["token"]),
            root=_optional_path(data.get("root")) or path.parent,
            crest_threads=int(data.get("crest_threads", 4)),
            workers=None if data.get("workers") is None else int(data["workers"]),
            bin_dir=_optional_path(data.get("bin_dir")),
            crest_timeout_s=float(data.get("crest_timeout_s", 3600.0)),
            mopac_timeout_s=float(data.get("mopac_timeout_s", 900.0)),
            keep_artifacts=bool(data.get("keep_artifacts", False)),
        )
    except ConfigError:
        raise
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{path}: valor inválido: {exc}") from None
