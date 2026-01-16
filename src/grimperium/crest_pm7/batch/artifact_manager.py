"""Artifact manager for preserving debug/audit files.

This module provides ArtifactManager for:
- Copying CREST output files (xyz, logs) to artifact directory
- Copying MOPAC files (input, output, archive) to artifact directory
- Creating structured artifact directories per molecule
- Managing artifact retention and cleanup

Thread-safety:
- SAFE for concurrent access to DIFFERENT mol_ids
- NOT safe for concurrent access to the SAME mol_id

Usage:
    artifact_manager = ArtifactManager(artifact_dir)
    artifact_paths = artifact_manager.save_artifacts(
        mol_id="mol_001",
        batch_id="batch_001",
        crest_work_dir=crest_output_path,
        mopac_work_dir=mopac_output_path,
    )
"""

from __future__ import annotations

import json
import logging
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

LOG = logging.getLogger("grimperium.crest_pm7.batch.artifact_manager")


@dataclass
class ArtifactPaths:
    """Paths to saved artifacts for a molecule."""

    mol_id: str
    batch_id: str
    artifact_root: Path

    # CREST artifacts
    crest_dir: Path | None = None
    crest_best_xyz: Path | None = None
    crest_conformers_xyz: Path | None = None
    crest_output_log: Path | None = None
    crest_energies: Path | None = None

    # MOPAC artifacts (per conformer)
    mopac_dir: Path | None = None
    mopac_inputs: list[Path] = field(default_factory=list)
    mopac_outputs: list[Path] = field(default_factory=list)
    mopac_archives: list[Path] = field(default_factory=list)

    # Metadata
    manifest_file: Path | None = None
    timestamp: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "mol_id": self.mol_id,
            "batch_id": self.batch_id,
            "artifact_root": str(self.artifact_root),
            "crest": {
                "dir": str(self.crest_dir) if self.crest_dir else None,
                "best_xyz": str(self.crest_best_xyz) if self.crest_best_xyz else None,
                "conformers_xyz": (
                    str(self.crest_conformers_xyz) if self.crest_conformers_xyz else None
                ),
                "output_log": (
                    str(self.crest_output_log) if self.crest_output_log else None
                ),
                "energies": str(self.crest_energies) if self.crest_energies else None,
            },
            "mopac": {
                "dir": str(self.mopac_dir) if self.mopac_dir else None,
                "inputs": [str(p) for p in self.mopac_inputs],
                "outputs": [str(p) for p in self.mopac_outputs],
                "archives": [str(p) for p in self.mopac_archives],
            },
            "manifest_file": str(self.manifest_file) if self.manifest_file else None,
            "timestamp": self.timestamp.isoformat(),
        }


class ArtifactManager:
    """Manages artifact preservation for debug and audit purposes.

    Creates structured directories for each molecule containing:
    - CREST output files (geometries, logs, energies)
    - MOPAC files (input, output, archive)
    - Manifest with file listing and metadata

    Directory Structure:
        artifact_dir/
        ├── batch_001/
        │   ├── mol_00001/
        │   │   ├── manifest.json
        │   │   ├── crest/
        │   │   │   ├── crest_best.xyz
        │   │   │   ├── crest_conformers.xyz
        │   │   │   ├── crest.output
        │   │   │   └── crest.energies
        │   │   └── mopac/
        │   │       ├── conf_0.mop
        │   │       ├── conf_0.out
        │   │       ├── conf_0.arc
        │   │       └── ...
        │   └── mol_00002/
        │       └── ...
        └── batch_002/
            └── ...

    Attributes:
        artifact_dir: Root directory for artifacts
        preserve_on_success: Whether to preserve artifacts for successful molecules
        preserve_on_failure: Whether to preserve artifacts for failed molecules
    """

    # CREST files to preserve
    CREST_FILES = [
        ("crest_best.xyz", "crest_best_xyz"),
        ("crest_conformers.xyz", "crest_conformers_xyz"),
        ("crest.output", "crest_output_log"),
        ("crest.engrad", "crest_energies"),
        ("crest.energies", "crest_energies"),  # Alternative name
    ]

    # MOPAC file patterns
    MOPAC_INPUT_PATTERN = "*.mop"
    MOPAC_OUTPUT_PATTERN = "*.out"
    MOPAC_ARCHIVE_PATTERN = "*.arc"

    def __init__(
        self,
        artifact_dir: Path,
        preserve_on_success: bool = True,
        preserve_on_failure: bool = True,
    ) -> None:
        """Initialize artifact manager.

        Args:
            artifact_dir: Root directory for artifacts (created if needed)
            preserve_on_success: Save artifacts for successful molecules
            preserve_on_failure: Save artifacts for failed molecules
        """
        self.artifact_dir = Path(artifact_dir)
        self.artifact_dir.mkdir(parents=True, exist_ok=True)
        self.preserve_on_success = preserve_on_success
        self.preserve_on_failure = preserve_on_failure
        LOG.debug(f"ArtifactManager initialized at {self.artifact_dir}")

    # Windows reserved names (case-insensitive)
    _WINDOWS_RESERVED = frozenset({
        "CON", "PRN", "AUX", "NUL",
        "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9",
        "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9",
    })

    def _sanitize_id(self, identifier: str) -> str:
        """Sanitize molecule or batch ID to prevent path injection attacks.

        Handles:
        - Directory traversal (../)
        - Null bytes (\\x00)
        - Path separators (/, \\)
        - Shell metacharacters (*, ?, <, >, |, :, ")
        - Windows reserved names (CON, PRN, etc.)
        - Filesystem length limits (255 chars)

        Args:
            identifier: mol_id or batch_id string

        Returns:
            Safe, filesystem-compatible string
        """
        if not identifier:
            return "_empty_"

        # Remove/replace dangerous characters: / \ : * ? " < > | and null byte
        safe = re.sub(r'[/\\:*?"<>|\x00]', "_", identifier)

        # Remove .. (prevent directory traversal)
        safe = safe.replace("..", "__")

        # Check for Windows reserved names (with or without extension)
        name_part = safe.split(".")[0].upper()
        if name_part in self._WINDOWS_RESERVED:
            safe = f"_{safe}"

        # Truncate to 255 chars (filesystem limit)
        if len(safe) > 255:
            safe = safe[:255]

        return safe

    def get_mol_artifact_dir(self, mol_id: str, batch_id: str) -> Path:
        """Get artifact directory path for a molecule.

        Args:
            mol_id: Molecule identifier
            batch_id: Batch identifier

        Returns:
            Path to molecule artifact directory
        """
        safe_mol_id = self._sanitize_id(mol_id)
        safe_batch_id = self._sanitize_id(batch_id)
        return self.artifact_dir / safe_batch_id / safe_mol_id

    def save_artifacts(
        self,
        mol_id: str,
        batch_id: str,
        crest_work_dir: Path | None = None,
        mopac_work_dir: Path | None = None,
        success: bool = True,
        extra_metadata: dict[str, Any] | None = None,
    ) -> ArtifactPaths | None:
        """Save artifacts for a processed molecule.

        Args:
            mol_id: Molecule identifier
            batch_id: Batch identifier
            crest_work_dir: CREST working directory (contains output files)
            mopac_work_dir: MOPAC working directory (contains conformer files)
            success: Whether processing succeeded
            extra_metadata: Additional metadata to include in manifest

        Returns:
            ArtifactPaths with saved file locations, or None if skipped
        """
        # Check if we should preserve based on success status
        if success and not self.preserve_on_success:
            LOG.debug(f"Skipping artifacts for {mol_id}: success=True, preserve=False")
            return None

        if not success and not self.preserve_on_failure:
            LOG.debug(f"Skipping artifacts for {mol_id}: success=False, preserve=False")
            return None

        # Create artifact directory
        mol_artifact_dir = self.get_mol_artifact_dir(mol_id, batch_id)
        mol_artifact_dir.mkdir(parents=True, exist_ok=True)

        artifact_paths = ArtifactPaths(
            mol_id=mol_id,
            batch_id=batch_id,
            artifact_root=mol_artifact_dir,
        )

        # Save CREST artifacts
        if crest_work_dir and crest_work_dir.exists():
            self._save_crest_artifacts(crest_work_dir, mol_artifact_dir, artifact_paths)

        # Save MOPAC artifacts
        if mopac_work_dir and mopac_work_dir.exists():
            self._save_mopac_artifacts(mopac_work_dir, mol_artifact_dir, artifact_paths)

        # Save manifest
        manifest_path = mol_artifact_dir / "manifest.json"
        manifest_data = artifact_paths.to_dict()
        manifest_data["success"] = success
        if extra_metadata:
            manifest_data["extra"] = extra_metadata

        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest_data, f, indent=2)

        artifact_paths.manifest_file = manifest_path

        LOG.info(f"Saved artifacts for {mol_id} to {mol_artifact_dir}")
        return artifact_paths

    def _save_crest_artifacts(
        self,
        crest_work_dir: Path,
        mol_artifact_dir: Path,
        artifact_paths: ArtifactPaths,
    ) -> None:
        """Copy CREST output files to artifact directory.

        Args:
            crest_work_dir: Source CREST working directory
            mol_artifact_dir: Destination artifact directory
            artifact_paths: ArtifactPaths to update
        """
        crest_dest = mol_artifact_dir / "crest"
        crest_dest.mkdir(exist_ok=True)
        artifact_paths.crest_dir = crest_dest

        for filename, attr_name in self.CREST_FILES:
            src_file = crest_work_dir / filename
            if src_file.exists():
                dest_file = crest_dest / filename
                try:
                    shutil.copy2(src_file, dest_file)
                    setattr(artifact_paths, attr_name, dest_file)
                    LOG.debug(f"Copied {filename} to {dest_file}")
                except OSError as e:
                    LOG.warning(f"Failed to copy {src_file}: {e}")

    def _save_mopac_artifacts(
        self,
        mopac_work_dir: Path,
        mol_artifact_dir: Path,
        artifact_paths: ArtifactPaths,
    ) -> None:
        """Copy MOPAC files to artifact directory.

        Args:
            mopac_work_dir: Source MOPAC working directory
            mol_artifact_dir: Destination artifact directory
            artifact_paths: ArtifactPaths to update
        """
        mopac_dest = mol_artifact_dir / "mopac"
        mopac_dest.mkdir(exist_ok=True)
        artifact_paths.mopac_dir = mopac_dest

        # Copy input files (.mop)
        for src_file in mopac_work_dir.glob(self.MOPAC_INPUT_PATTERN):
            dest_file = mopac_dest / src_file.name
            try:
                shutil.copy2(src_file, dest_file)
                artifact_paths.mopac_inputs.append(dest_file)
                LOG.debug(f"Copied {src_file.name} to {dest_file}")
            except OSError as e:
                LOG.warning(f"Failed to copy {src_file}: {e}")

        # Copy output files (.out)
        for src_file in mopac_work_dir.glob(self.MOPAC_OUTPUT_PATTERN):
            dest_file = mopac_dest / src_file.name
            try:
                shutil.copy2(src_file, dest_file)
                artifact_paths.mopac_outputs.append(dest_file)
                LOG.debug(f"Copied {src_file.name} to {dest_file}")
            except OSError as e:
                LOG.warning(f"Failed to copy {src_file}: {e}")

        # Copy archive files (.arc)
        for src_file in mopac_work_dir.glob(self.MOPAC_ARCHIVE_PATTERN):
            dest_file = mopac_dest / src_file.name
            try:
                shutil.copy2(src_file, dest_file)
                artifact_paths.mopac_archives.append(dest_file)
                LOG.debug(f"Copied {src_file.name} to {dest_file}")
            except OSError as e:
                LOG.warning(f"Failed to copy {src_file}: {e}")

    def load_manifest(self, mol_id: str, batch_id: str) -> dict[str, Any] | None:
        """Load manifest for a molecule.

        Args:
            mol_id: Molecule identifier
            batch_id: Batch identifier

        Returns:
            Manifest data dict or None if not found
        """
        manifest_path = self.get_mol_artifact_dir(mol_id, batch_id) / "manifest.json"
        if not manifest_path.exists():
            return None

        try:
            with open(manifest_path, encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            LOG.warning(f"Failed to load manifest for {mol_id}: {e}")
            return None

    def list_batch_artifacts(self, batch_id: str) -> list[str]:
        """List all mol_ids with artifacts for a batch.

        Args:
            batch_id: Batch identifier

        Returns:
            List of mol_ids
        """
        safe_batch_id = self._sanitize_id(batch_id)
        batch_dir = self.artifact_dir / safe_batch_id

        if not batch_dir.exists():
            return []

        return [d.name for d in batch_dir.iterdir() if d.is_dir()]

    def delete_artifacts(self, mol_id: str, batch_id: str) -> bool:
        """Delete artifacts for a molecule.

        Args:
            mol_id: Molecule identifier
            batch_id: Batch identifier

        Returns:
            True if deleted, False if not found
        """
        mol_artifact_dir = self.get_mol_artifact_dir(mol_id, batch_id)

        if not mol_artifact_dir.exists():
            return False

        try:
            shutil.rmtree(mol_artifact_dir)
            LOG.info(f"Deleted artifacts for {mol_id}")
            return True
        except OSError as e:
            LOG.error(f"Failed to delete artifacts for {mol_id}: {e}")
            return False

    def get_artifact_stats(self, batch_id: str | None = None) -> dict[str, Any]:
        """Get statistics about stored artifacts.

        Args:
            batch_id: Optional batch to filter by

        Returns:
            Dict with artifact counts and sizes
        """
        if batch_id:
            safe_batch_id = self._sanitize_id(batch_id)
            search_dirs = [self.artifact_dir / safe_batch_id]
        else:
            search_dirs = [
                d for d in self.artifact_dir.iterdir() if d.is_dir()
            ]

        total_molecules = 0
        total_files = 0
        total_size_bytes = 0
        crest_files = 0
        mopac_files = 0

        for batch_dir in search_dirs:
            if not batch_dir.exists():
                continue

            for mol_dir in batch_dir.iterdir():
                if not mol_dir.is_dir():
                    continue

                total_molecules += 1

                # Count CREST files
                crest_dir = mol_dir / "crest"
                if crest_dir.exists():
                    for f in crest_dir.iterdir():
                        if f.is_file():
                            crest_files += 1
                            total_files += 1
                            total_size_bytes += f.stat().st_size

                # Count MOPAC files
                mopac_dir = mol_dir / "mopac"
                if mopac_dir.exists():
                    for f in mopac_dir.iterdir():
                        if f.is_file():
                            mopac_files += 1
                            total_files += 1
                            total_size_bytes += f.stat().st_size

        return {
            "total_molecules": total_molecules,
            "total_files": total_files,
            "total_size_mb": round(total_size_bytes / (1024 * 1024), 2),
            "crest_files": crest_files,
            "mopac_files": mopac_files,
        }

    def cleanup_batch(self, batch_id: str) -> int:
        """Delete all artifacts for a batch.

        Args:
            batch_id: Batch identifier

        Returns:
            Number of molecule directories deleted
        """
        safe_batch_id = self._sanitize_id(batch_id)
        batch_dir = self.artifact_dir / safe_batch_id

        if not batch_dir.exists():
            return 0

        count = 0
        for mol_dir in list(batch_dir.iterdir()):
            if mol_dir.is_dir():
                try:
                    shutil.rmtree(mol_dir)
                    count += 1
                except OSError as e:
                    LOG.warning(f"Failed to delete {mol_dir}: {e}")

        # Try to remove the batch directory if empty
        try:
            batch_dir.rmdir()
        except OSError:
            pass  # Not empty or other error

        LOG.info(f"Cleaned up {count} molecule artifacts from batch {batch_id}")
        return count
