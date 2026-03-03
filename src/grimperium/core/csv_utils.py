"""Atomic CSV write utility.

Provides ``atomic_to_csv`` which writes a DataFrame to a CSV file using the
temp-file-then-rename pattern, preventing data loss on Ctrl-C or crashes.
A ``.bak`` backup is created automatically before overwriting an existing file.
"""

import logging
import os
import shutil
import tempfile
from pathlib import Path

import pandas as pd

LOG = logging.getLogger(__name__)


def atomic_to_csv(path: Path, df: pd.DataFrame) -> None:
    """Write *df* to *path* atomically with automatic backup.

    Algorithm:
        1. Write to a temp file in the same directory (same filesystem).
        2. ``flush`` + ``fsync`` to force data to disk.
        3. If *path* already exists and is non-empty, copy it to
           ``<path>.bak`` so it can be recovered later.
        4. ``os.replace`` (atomic on POSIX when same filesystem).
        5. On failure the temp file is removed and the original is untouched.

    Args:
        path: Destination CSV file path.
        df: DataFrame to persist.

    Raises:
        TypeError: If *path* or *df* is ``None``.
    """
    if path is None:
        raise TypeError("path must not be None")
    if df is None:
        raise TypeError("df must not be None")

    path.parent.mkdir(parents=True, exist_ok=True)

    fd, tmp_path_str = tempfile.mkstemp(
        dir=path.parent, suffix=".csv.tmp", prefix=".atomic_"
    )
    tmp_path = Path(tmp_path_str)

    try:
        # Write DataFrame to temp file with fsync
        f = None
        try:
            f = os.fdopen(fd, "w", newline="", encoding="utf-8")
        except Exception:
            os.close(fd)
            raise

        with f:
            df.to_csv(f, index=False)
            f.flush()
            os.fsync(f.fileno())

        # Backup existing file (if non-empty) before replacing
        if path.exists() and path.stat().st_size > 0:
            backup_path = path.with_name(path.name + ".bak")
            shutil.copy2(path, backup_path)

        # Atomic rename
        os.replace(tmp_path, path)
        LOG.debug("Atomically wrote %d rows to %s", len(df), path)

    except Exception:
        # Cleanup temp file on any failure
        try:
            tmp_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise
