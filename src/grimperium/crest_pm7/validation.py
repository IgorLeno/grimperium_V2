"""Environment validation for CREST-PM7 Pipeline.

Validates that all required external tools are available and properly
configured.
"""

import shutil
import subprocess
from dataclasses import dataclass, field

from .config import PM7Config


@dataclass
class ValidationResult:
    """Result of environment validation.

    Attributes:
        valid: Whether the environment is valid
        crest_available: Whether CREST executable is found
        crest_version: CREST version string if available
        mopac_available: Whether MOPAC executable is found
        mopac_version: MOPAC version string if available
        errors: List of validation errors
        warnings: List of validation warnings
    """

    valid: bool = True
    crest_available: bool = False
    crest_version: str | None = None
    mopac_available: bool = False
    mopac_version: str | None = None
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def check_executable(
    name: str, version_flag: str = "--version"
) -> tuple[bool, str | None]:
    """Check if an executable is available and get its version.

    Args:
        name: Executable name or path
        version_flag: Flag to get version

    Returns:
        Tuple of (available, version_string)
    """
    path = shutil.which(name)
    if path is None:
        return False, None

    try:
        result = subprocess.run(
            [path, version_flag],
            capture_output=True,
            text=True,
            timeout=10,
        )
        version_str = result.stdout.strip() or result.stderr.strip()
        # Take first line only
        version: str | None = version_str.split("\n")[0] if version_str else None
        return True, version
    except (subprocess.TimeoutExpired, subprocess.SubprocessError):
        # Subprocess failure means the executable is not functioning correctly
        return False, None


def validate_environment(config: PM7Config) -> ValidationResult:
    """Validate the execution environment.

    Checks for required external tools: CREST and MOPAC.

    Args:
        config: Pipeline configuration

    Returns:
        ValidationResult with status and details
    """
    result = ValidationResult()

    # Check CREST
    result.crest_available, result.crest_version = check_executable(
        config.crest_executable, "--version"
    )
    if not result.crest_available:
        result.valid = False
        result.errors.append(f"CREST executable not found: {config.crest_executable}")

    # Check MOPAC
    result.mopac_available, result.mopac_version = check_executable(
        config.mopac_executable, "-v"
    )
    if not result.mopac_available:
        result.valid = False
        result.errors.append(f"MOPAC executable not found: {config.mopac_executable}")

    # Check directories
    try:
        config.temp_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        result.valid = False
        result.errors.append(f"Cannot create temp directory {config.temp_dir}: {e}")

    try:
        config.output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        result.valid = False
        result.errors.append(f"Cannot create output directory {config.output_dir}: {e}")

    return result


def validate_environment_strict(config: PM7Config) -> ValidationResult:
    """Strict environment validation requiring all tools.

    Requires the configured external tools and writable directories.

    Args:
        config: Pipeline configuration

    Returns:
        ValidationResult with status and details
    """
    return validate_environment(config)
