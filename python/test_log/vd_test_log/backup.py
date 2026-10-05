"""SQLite-consistent backups with snapshot-derived attachment manifests."""

from __future__ import annotations

from datetime import datetime, timezone
import os
import re
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

    def _validated_backup_root(self) -> Path:
        """Return the configured backup root only when it stays inside the data root."""
        selected_root = Path(os.path.abspath(os.fspath(self.paths.root)))
        backup_root = Path(os.path.abspath(os.fspath(self.paths.backup_root)))
        try:
            relative = backup_root.relative_to(selected_root)
        except ValueError as error:
            raise BackupError("The backup folder must be inside the selected data folder.") from error
        if not relative.parts:
            raise BackupError("The backup folder must be below the selected data folder.")

        resolved_root = Path(self.paths.root).resolve(strict=False)
        resolved_backup = Path(self.paths.backup_root).resolve(strict=False)
        try:
            resolved_backup.relative_to(resolved_root)
        except ValueError as error:
            raise BackupError("The resolved backup folder escapes the selected data folder.") from error
        if resolved_backup == resolved_root:
            raise BackupError("The backup folder must be below the selected data folder.")

        current = selected_root
        for part in relative.parts:
            current /= part
            if _is_link_or_junction(current):
                raise BackupError("The backup folder cannot pass through a link or junction.")
        return backup_root

    def _validated_backup_path(self, path: Path) -> Path:
        """Check a generated backup path against the currently validated root."""
        backup_root = self._validated_backup_root()
        candidate = Path(os.path.abspath(os.fspath(path)))
        try:
            relative = candidate.relative_to(backup_root)
        except ValueError as error:
            raise BackupError("A generated backup path escapes the backup folder.") from error

        current = backup_root
        for part in relative.parts:
            current /= part
            if _is_link_or_junction(current):
                raise BackupError("A generated backup path cannot pass through a link or junction.")

        resolved_root = Path(self.paths.root).resolve(strict=False)
        resolved_backup = backup_root.resolve(strict=False)
        resolved_candidate = candidate.resolve(strict=False)
        try:
            resolved_backup.relative_to(resolved_root)
            resolved_candidate.relative_to(resolved_backup)
        except ValueError as error:
            raise BackupError("A generated backup path resolves outside its allowed folder.") from error
        return candidate

    def create_backup(self, connection: sqlite3.Connection, attachment_store: AttachmentStore) -> Path:
        backup_root = self._validated_backup_root()
        temporary_path: Path | None = None
        temporary_name: str | None = None
        try:
            # Validate before mkdir so malformed DataPaths cannot create an outside folder.
            backup_root.mkdir(parents=True, exist_ok=True)
            backup_root = self._validated_backup_root()

            temporary_name = f".tmp-backup-{new_id()}"
            temporary_path = self._validated_backup_path(backup_root / temporary_name)
            temporary_path.mkdir()
            temporary_path = self._validated_backup_path(temporary_path)
            attachments_path = self._validated_backup_path(temporary_path / "attachments")
            attachments_path.mkdir()
            self._validated_backup_path(attachments_path)

            database_name = Path(self.paths.database).name
            if database_name in ("", ".", ".."):
                raise BackupError("The database path must include a valid file name.")
            database_copy = self._validated_backup_path(temporary_path / database_name)
            snapshot = sqlite3.connect(database_copy)
            try:
                connection.backup(snapshot)
            finally:
                snapshot.close()
            self._validated_backup_path(database_copy)

            manifest = self._snapshot_manifest(database_copy)
            for relative_path in manifest:
                source = attachment_store.resolve_relative(relative_path, require_exists=True)
                destination = self._backup_attachment_path(temporary_path, relative_path)
                destination = self._validated_backup_path(destination)
                destination.parent.mkdir(parents=True, exist_ok=True)
                self._validated_backup_path(destination.parent)
                destination = self._validated_backup_path(destination)
                shutil.copyfile(source, destination)
                self._validated_backup_path(destination)

            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
            final_name = f"{stamp}-{new_id()[:8]}"
            final_path = self._validated_backup_path(backup_root / final_name)
            if final_path.exists():
                raise BackupError("The generated backup destination already exists.")
            temporary_path = self._validated_backup_path(temporary_path)
            temporary_path.rename(final_path)
            return self._validated_backup_path(final_path)
        except BaseException as error:
            cleanup_error = self._remove_temporary_folder(temporary_path, temporary_name)
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
        candidate = temporary_root / relative
        resolved_candidate = candidate.resolve(strict=False)
        try:
            resolved_candidate.relative_to(attachment_root)
        except ValueError as error:
            raise AttachmentPathError(
                "A snapshot attachment path escapes the backup attachments folder."
            ) from error
        if resolved_candidate == attachment_root:
            raise AttachmentPathError("A snapshot attachment path must name a file.")
        return candidate

    def _remove_temporary_folder(
        self,
        temporary_path: Path | None,
        temporary_name: str | None,
    ) -> OSError | None:
        if temporary_path is None or temporary_name is None:
            return None
        try:
            if not re.fullmatch(r"\.tmp-backup-[0-9a-f]{32}", temporary_name):
                return OSError("Temporary backup name is not a generated folder name.")
            backup_root = self._validated_backup_root()
            expected = backup_root / temporary_name
            if Path(os.path.abspath(os.fspath(temporary_path))) != expected:
                return OSError("Temporary backup path no longer matches the generated folder.")
            target = self._validated_backup_path(expected)
            if target.exists():
                if not target.is_dir():
                    return OSError("Temporary backup path is not a folder.")
                # Recheck confinement and link status immediately before recursive cleanup.
                target = self._validated_backup_path(target)
                shutil.rmtree(target)
        except OSError as error:
            return error
        return None


def _is_link_or_junction(path: Path) -> bool:
    if path.is_symlink():
        return True
    is_junction = getattr(os.path, "isjunction", None)
    return bool(is_junction and is_junction(path))
