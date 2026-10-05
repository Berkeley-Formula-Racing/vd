"""Application data paths and an operating-system-backed single-instance lock."""

from __future__ import annotations

import errno
import os
from pathlib import Path
from typing import Mapping

if os.name == "nt":
    import msvcrt
else:  # Allows the storage package's tests to run on non-Windows hosts too.
    import fcntl

from .models import DataPaths


class AlreadyRunningError(RuntimeError):
    """The selected data root is already held by another process."""


class DataRootLock:
    """Hold an exclusive OS lock on one data root for the app's lifetime."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._handle = None

    def acquire(self) -> None:
        if self._handle is not None:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+b")
        try:
            handle.seek(0, os.SEEK_END)
            if handle.tell() == 0:
                handle.write(b"\0")
                handle.flush()
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            handle.close()
            if _is_lock_contention(error):
                raise AlreadyRunningError(
                    f"The test log data folder is already open: {self.path.parent}"
                ) from error
            raise
        self._handle = handle

    def close(self) -> None:
        handle = self._handle
        if handle is None:
            return
        self._handle = None
        try:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        except OSError:
            # Closing the held descriptor releases OS locks even if explicit unlock fails.
            pass
        finally:
            handle.close()

    def __enter__(self) -> DataRootLock:
        self.acquire()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def _is_lock_contention(error: OSError) -> bool:
    if os.name == "nt":
        return getattr(error, "winerror", None) in (33, 36) or error.errno in (
            errno.EACCES,
            errno.EDEADLK,
        )
    return error.errno in (errno.EACCES, errno.EAGAIN, errno.EDEADLK)


def resolve_data_paths(
    data_dir: Path | None,
    environ: Mapping[str, str],
    home: Path,
) -> DataPaths:
    """Resolve the explicit, per-user, or home-fallback data root."""
    if data_dir is not None:
        root = Path(data_dir).expanduser()
    else:
        local_app_data = environ.get("LOCALAPPDATA")
        if local_app_data and str(local_app_data).strip():
            root = Path(local_app_data).expanduser() / "VDTestLog"
        else:
            root = Path(home).expanduser() / "VDTestLog"
    resolved_root = root.resolve(strict=False)
    return DataPaths(
        root=resolved_root,
        database=resolved_root / "test_log.sqlite3",
        attachments=resolved_root / "attachments",
        lock_file=resolved_root / ".test_log.lock",
        backup_root=resolved_root / "backups",
    )
