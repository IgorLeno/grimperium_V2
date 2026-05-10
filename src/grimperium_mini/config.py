"""Configuration for grimperium_mini."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


def _env_int(name: str, default: int) -> int:
    value = os.getenv(name)
    if value is None or value == "":
        return default
    return int(value)


def _env_float(name: str, default: float) -> float:
    value = os.getenv(name)
    if value is None or value == "":
        return default
    return float(value)


def _env_path(name: str, default: Path) -> Path:
    value = os.getenv(name)
    if value is None or value == "":
        return default
    return Path(value)


@dataclass
class MiniConfig:
    """Runtime settings for the distilled TCC pipeline."""

    crest_executable: str = field(
        default_factory=lambda: os.getenv("GRIMPERIUM_MINI_CREST", "crest")
    )
    mopac_executable: str = field(
        default_factory=lambda: os.getenv("GRIMPERIUM_MINI_MOPAC", "mopac")
    )
    xtb_executable: str = field(
        default_factory=lambda: os.getenv("GRIMPERIUM_MINI_XTB", "xtb")
    )
    xtb_enabled: bool = field(
        default_factory=lambda: os.getenv("GRIMPERIUM_MINI_XTB_ENABLED", "true").lower()
        != "false"
    )
    threads: int = field(default_factory=lambda: _env_int("GRIMPERIUM_MINI_THREADS", 8))
    crest_method: str = field(
        default_factory=lambda: os.getenv("GRIMPERIUM_MINI_CREST_METHOD", "gfn2")
    )
    crest_ewin: float = field(
        default_factory=lambda: _env_float("GRIMPERIUM_MINI_CREST_EWIN", 6.0)
    )
    crest_rthr: float = field(
        default_factory=lambda: _env_float("GRIMPERIUM_MINI_CREST_RTHR", 0.125)
    )
    crest_quick_mode: str | None = field(
        default_factory=lambda: os.getenv("GRIMPERIUM_MINI_CREST_QUICK") or None
    )
    mopac_keywords_extra: list[str] = field(
        default_factory=lambda: ["PRECISE", "EF", "AUX"]
    )
    max_conformers_to_optimize: int | None = field(
        default_factory=lambda: _env_int("GRIMPERIUM_MINI_MAX_CONFORMERS", 10) or None
    )
    work_root: Path = field(
        default_factory=lambda: _env_path(
            "GRIMPERIUM_MINI_WORK_ROOT",
            Path("/home/ifernandes/Projetos/grimperium/runs"),
        )
    )
    results_dir: Path = field(
        default_factory=lambda: _env_path(
            "GRIMPERIUM_MINI_RESULTS_DIR",
            Path("/home/ifernandes/Projetos/grimperium/results"),
        )
    )
    timeout_xtb_s: int = field(
        default_factory=lambda: _env_int("GRIMPERIUM_MINI_TIMEOUT_XTB_S", 600)
    )
    timeout_crest_s: int = field(
        default_factory=lambda: _env_int("GRIMPERIUM_MINI_TIMEOUT_CREST_S", 82800)
    )
    timeout_mopac_s: int = field(
        default_factory=lambda: _env_int("GRIMPERIUM_MINI_TIMEOUT_MOPAC_S", 3600)
    )

    def __post_init__(self) -> None:
        self.crest_method = self.crest_method.lower().strip()
        if self.crest_quick_mode:
            self.crest_quick_mode = self.crest_quick_mode.lower().strip()
        if self.threads <= 0:
            raise ValueError("threads must be > 0")
        if self.crest_method not in {"gfn2", "gfnff", "gfn2//gfnff"}:
            raise ValueError("crest_method must be gfn2, gfnff, or gfn2//gfnff")
        if self.crest_quick_mode not in {None, "quick", "squick", "mquick"}:
            raise ValueError("crest_quick_mode must be quick, squick, mquick, or None")
