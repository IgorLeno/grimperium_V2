"""The cell package: one zip with everything a new cell needs to install.

Contents (under one folder, so "Extrair tudo" on Windows keeps it together):

- ``conda-lock.yml`` / ``environment.yml``: CREST, xtb and MOPAC exactly as
  on the central; the installer builds ``grimperium-quimica`` from the lock,
  so ``environment_sha256`` matches and registration is accepted.
- ``wheels/``: one wheel of the committed ``HEAD`` (``grimperium`` and
  ``semi_imperium`` together: the cell's code imports both).
- ``requirements-celula.txt`` + ``constraints-central.txt``: the Python
  libraries the cell imports, pinned to the versions the central runs.
- ``instalar.sh`` (Ubuntu/WSL2) and the Windows script it fills in.
- ``config.json``: central URL and the cell's token. **A secret**: the zip
  is generated on demand and never stored on the central.

The token exists only inside the zip. The central keeps just its hash, so a
package always comes with a *new* token (``CellRegistry.rotate_token``) and
the cell's previous token stops working. The wheel is built before the token
rotates: a failed build leaves the installed cell running.
"""

from __future__ import annotations

import importlib.metadata
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

#: What the cell imports (``semi_imperium.cell`` → ``grimperium``), by
#: distribution name; pip brings their dependencies under the constraints.
#: ``tests/webapi/test_cell_package.py`` fails when an import outgrows it.
CELL_DISTRIBUTIONS = (
    "httpx",
    "numpy",
    "pandas",
    "pydantic",
    "rdkit",
    "rich",
    "scikit-learn",
    "scipy",
)

#: Port of the cell API in the guide (``uvicorn … --port 8100``).
CELL_API_PORT = 8100

TEMPLATES = Path(__file__).with_name("celula")

# src/grimperium/webapi/cell_package.py -> repository root
_REPO_ROOT = Path(__file__).resolve().parents[3]
_SUBPROCESS_TIMEOUT_S = 600


class PackageError(RuntimeError):
    """The package cannot be built; the message says why, in Portuguese."""


@dataclass(frozen=True)
class WheelBuild:
    path: Path
    commit: str
    #: Tracked files differ from ``commit``; the wheel still holds ``commit``.
    dirty: bool


WheelBuilder = Callable[[Path], WheelBuild]


def _run(args: list[str], cwd: Path) -> subprocess.CompletedProcess[bytes]:
    try:
        return subprocess.run(  # noqa: S603  (fixed programs, no shell)
            args,
            cwd=cwd,
            capture_output=True,
            check=True,
            timeout=_SUBPROCESS_TIMEOUT_S,
        )
    except FileNotFoundError:
        raise PackageError(f"programa não encontrado: {args[0]}") from None
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.decode(errors="replace").strip().splitlines()[-1:]
        raise PackageError(
            f"{' '.join(args[:3])} falhou: {' '.join(detail) or exc.returncode}"
        ) from None
    except subprocess.TimeoutExpired:
        raise PackageError(f"{' '.join(args[:3])} demorou demais") from None


def build_wheel_from_head(out_dir: Path, repo_root: Path = _REPO_ROOT) -> WheelBuild:
    """Build the wheel of the committed ``HEAD``, never of uncommitted edits.

    ``git archive`` exports the tracked files of the commit into a scratch
    folder and ``poetry build`` packs them there, so what the cell runs is a
    commit anyone can check out again.
    """
    poetry = shutil.which("poetry")
    if poetry is None:
        raise PackageError("poetry não encontrado no PATH da API do dashboard")
    commit = _run(["git", "rev-parse", "HEAD"], repo_root).stdout.decode().strip()
    status = _run(
        ["git", "status", "--porcelain", "--untracked-files=no"], repo_root
    ).stdout
    archive = _run(["git", "archive", "--format=tar", commit], repo_root).stdout
    with tempfile.TemporaryDirectory(prefix="grimperium-wheel-") as scratch:
        source = Path(scratch) / "src"
        source.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
            tar.extractall(source, filter="data")
        _run([poetry, "build", "--format", "wheel", "--output", str(out_dir)], source)
    wheels = sorted(out_dir.glob("*.whl"))
    if len(wheels) != 1:
        raise PackageError(f"poetry build gerou {len(wheels)} wheels em {out_dir}")
    return WheelBuild(path=wheels[0], commit=commit, dirty=bool(status.strip()))


def python_series() -> str:
    """``3.14``: the cell gets the central's Python series from conda-forge."""
    return f"{sys.version_info.major}.{sys.version_info.minor}"


def pinned_requirements(names: tuple[str, ...] = CELL_DISTRIBUTIONS) -> str:
    """``name==version`` of each library the cell imports, as the central has it."""
    lines = [f"{name}=={importlib.metadata.version(name)}" for name in names]
    return "\n".join(lines) + "\n"


def central_constraints() -> str:
    """Every distribution installed on the central, pinned.

    As a pip constraints file it fixes the dependencies pip pulls in for
    ``pinned_requirements`` without forcing anything else onto the cell.
    The project itself (installed from the checkout) is left out.
    """
    pins: dict[str, str] = {}
    for dist in importlib.metadata.distributions():
        name = dist.metadata["Name"]
        if not name or name.lower() == "grimperium":
            continue
        pins.setdefault(name.lower(), f"{name}=={dist.version}")
    return "\n".join(pins[key] for key in sorted(pins)) + "\n"


def suggest_central_url(env: Mapping[str, str] | None = None) -> str:
    """Where cells should call: ``GRIMPERIUM_CELL_URL``, else this Tailscale IP.

    Empty when neither is known; the dashboard then asks for it.
    """
    environ = os.environ if env is None else env
    explicit = environ.get("GRIMPERIUM_CELL_URL", "").strip()
    if explicit:
        return explicit.rstrip("/")
    tailscale = shutil.which("tailscale")
    if tailscale is None:
        return ""
    try:
        completed = subprocess.run(  # noqa: S603
            [tailscale, "ip", "-4"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    address = completed.stdout.strip().splitlines()[:1]
    if completed.returncode != 0 or not address:
        return ""
    return f"http://{address[0]}:{CELL_API_PORT}"


def _slug(name: str) -> str:
    ascii_name = name.encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", "-", ascii_name).strip("-")[:40] or "celula"


def package_folder(cell_id: str, name: str) -> str:
    """Folder (and zip) name: readable, and unique per cell."""
    return f"celula-{_slug(name)}-{cell_id.removeprefix('cell-')}"


def _add(archive: zipfile.ZipFile, path: str, data: bytes, mode: int) -> None:
    info = zipfile.ZipInfo(path, date_time=datetime.now().timetuple()[:6])
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = (0o100000 | mode) << 16
    archive.writestr(info, data)


@dataclass(frozen=True)
class CellPackage:
    cell_id: str
    cell_name: str
    central_url: str
    token: str
    wheel: WheelBuild
    #: ``probe_programs`` of the central, without install paths: the
    #: installer compares the cell's programs to these before the first start.
    central_programs: Mapping[str, str]
    repo_root: Path = _REPO_ROOT

    def build(self) -> bytes:
        """The zip, in memory (it carries the token: never written to disk here)."""
        folder = package_folder(self.cell_id, self.cell_name)
        manifest = {
            "cell_id": self.cell_id,
            "cell_name": self.cell_name,
            "central_url": self.central_url,
            "commit": self.wheel.commit,
            "commit_dirty": self.wheel.dirty,
            "built_at": datetime.now(timezone.utc).isoformat(),
            "python": python_series(),
            "wheel": self.wheel.path.name,
            "central_programs": {
                k: v
                for k, v in self.central_programs.items()
                if not k.endswith("_path")
            },
        }
        shell_vars = {
            "PYTHON_VERSION": python_series(),
            "WHEEL": f"wheels/{self.wheel.path.name}",
            "CENTRAL_URL": self.central_url,
            "CELL_NAME": self.cell_name,
        }
        config = {"central_url": self.central_url, "token": self.token}
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as archive:

            def add(path: str, data: bytes, mode: int = 0o644) -> None:
                _add(archive, f"{folder}/{path}", data, mode)

            add("LEIA-ME.txt", (TEMPLATES / "LEIA-ME.txt").read_bytes())
            add("instalar.sh", (TEMPLATES / "instalar.sh").read_bytes(), 0o755)
            add(
                "windows/configurar-windows.ps1",
                (TEMPLATES / "configurar-windows.ps1").read_bytes(),
            )
            add("config.json", (json.dumps(config, indent=2) + "\n").encode(), 0o600)
            add(
                "pacote.json",
                (json.dumps(manifest, indent=2, ensure_ascii=False) + "\n").encode(),
            )
            add(
                "pacote.env",
                "".join(
                    f"{key}={shlex.quote(value)}\n" for key, value in shell_vars.items()
                ).encode(),
            )
            for lock in ("conda-lock.yml", "environment.yml"):
                add(lock, (self.repo_root / lock).read_bytes())
            add("requirements-celula.txt", pinned_requirements().encode())
            add("constraints-central.txt", central_constraints().encode())
            add(f"wheels/{self.wheel.path.name}", self.wheel.path.read_bytes())
        return buffer.getvalue()
