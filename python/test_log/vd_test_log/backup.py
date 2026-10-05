"""SQLite-consistent backups with snapshot-derived attachment manifests."""

from __future__ import annotations

from datetime import datetime, timezone
import shutil
import sqlite3
from pathlib import Path

from .attachments import AttachmentPathError, AttachmentStore
from .models import DataPaths, new_id


class BackupError(OSError):
    """A backup could not be completed without replacing earlier backups."""


class BackupManager:
    def __init__(self, paths: DataPaths):
        self.paths = paths

    def create_backup(self, connection: sqlite3.Connection, attachment_store: AttachmentStore) -> Path:
        backup_root = Path(self.paths.backup_root)
        temporary_path: Path | None = None
        temporary_name: str | None = None
        try:
            backup_root.mkdir(parents=True, exist_ok=True)
            backup_root = backup_root.resolve(strict=False)
            temporary_name = f".tmp-backup-{new_id()}"
            temporary_path = backup_root / temporary_name
            temporary_path.mkdir()
            (temporary_path / "attachments").mkdir()

            database_copy = temporary_path / Path(self.paths.database).name
            snapshot = sqlite3.connect(database_copy)
            try:
                connection.backup(snapshot)
            finally:
                snapshot.close()

            manifest = self._snapshot_manifest(database_copy)
            for relative_path in manifest:
                source = attachment_store.resolve_relative(relative_path, require_exists=True)
                destination = self._backup_attachment_path(temporary_path, relative_path)
                destination.parent.mkdir(parents=True, exist_ok=True)
                # Recheck after creating folders, immediately before writing the copy.
                destination = self._backup_attachment_path(temporary_path, relative_path)
                shutil.copyfile(source, destination)

            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
            final_path = backup_root / f"{stamp}-{new_id()[:8]}"
            temporary_path.rename(final_path)
            return final_path
        except BaseException as error:
            cleanup_error = self._remove_temporary_folder(
                backup_root, temporary_path, temporary_name
            )
            if cleanup_error is not None and hasattr(error, "add_note"):
                error.add_note(f"Temporary backup cleanup also failed: {cleanup_error}")
            if isinstance(error, BackupError):
                raise
            if isinstance(error, Exception):
                raise BackupError(f"Could not create backup: {error}") from error
            raise

    @staticmethod
    def _snapshot_manifest(database_path: Path) -> list[str]:
        snapshot = sqlite3.connect(database_path)
        try:
            rows = snapshot.execute(
                "SELECT relative_path FROM attachments ORDER BY relative_path"
            )
            return [row[0] for row in rows]
        finally:
            snapshot.close()

    @staticmethod
    def _backup_attachment_path(temporary_root: Path, relative_path: str) -> Path:
        relative = Path(relative_path)
        if relative.is_absolute() or relative.drive:
            raise AttachmentPathError("Backup attachment paths must be relative.")
        attachment_root = (temporary_root / "attachments").resolve(strict=False)
        candidate = (temporary_root / relative).resolve(strict=False)
        try:
            candidate.relative_to(attachment_root)
        except ValueError as error:
            raise AttachmentPathError(
                "A snapshot attachment path escapes the backup attachments folder."
            ) from error
        if candidate == attachment_root:
            raise AttachmentPathError("A snapshot attachment path must name a file.")
        return candidate

    @staticmethod
    def _remove_temporary_folder(
        backup_root: Path,
        temporary_path: Path | None,
        temporary_name: str | None,
    ) -> OSError | None:
        if temporary_path is None or temporary_name is None:
            return None
        try:
            resolved_root = Path(backup_root).resolve(strict=False)
            expected = resolved_root / temporary_name
            if (
                temporary_path.name != temporary_name
                or temporary_path.parent.resolve(strict=False) != resolved_root
                or temporary_path.resolve(strict=False) != expected
            ):
                return OSError("Temporary backup path no longer matches the generated folder.")
            if temporary_path.is_symlink():
                return OSError("Temporary backup folder was replaced by a symlink.")
            if temporary_path.exists():
                shutil.rmtree(temporary_path)
        except OSError as error:
            return error
        return None
