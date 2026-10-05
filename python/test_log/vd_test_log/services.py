"""UI-facing service facade and lifecycle ownership for the local test log."""

from __future__ import annotations

import os
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from threading import Condition, local
from typing import Iterator, Sequence

from .attachments import AttachmentStore
from .backup import BackupManager
from .csv_export import export_setup as write_setup_csv
from .data_folder import DataRootLock
from .models import (
    Attachment,
    AttachmentRole,
    DataPaths,
    EventLayout,
    Lap,
    Setup,
    StagedAttachment,
    TestDay,
    new_id,
    utc_now_iso,
)
from .repository import SQLiteRepository
from .validation import ValidationError


class UnsafeCsvDestinationError(ValueError):
    """A CSV export would overwrite data owned by the test log."""


def _resolved_path(path: Path) -> Path:
    return Path(path).expanduser().resolve(strict=False)


def _same_existing_file(first: Path, second: Path) -> bool:
    try:
        return first.exists() and second.exists() and os.path.samefile(first, second)
    except OSError:
        return False


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _guard_csv_destination(paths: DataPaths, destination: Path) -> None:
    resolved_destination = _resolved_path(destination)
    database = Path(paths.database)
    protected_files = (
        database,
        Path(f"{database}-wal"),
        Path(f"{database}-shm"),
        Path(f"{database}-journal"),
        Path(paths.lock_file),
    )
    for protected_file in protected_files:
        resolved_protected = _resolved_path(protected_file)
        if resolved_destination == resolved_protected or _same_existing_file(
            resolved_destination, resolved_protected
        ):
            raise UnsafeCsvDestinationError(
                "Choose a CSV destination outside the test log database, sidecars, and lock file."
            )

    for owned_directory in (Path(paths.attachments), Path(paths.backup_root)):
        resolved_directory = _resolved_path(owned_directory)
        if _is_within(resolved_destination, resolved_directory):
            raise UnsafeCsvDestinationError(
                "Choose a CSV destination outside the test log attachments and backup folders."
            )


class BackupInProgressError(RuntimeError):
    """A mutation was requested while a consistent backup is being written."""


class AttachmentOwnershipError(ValueError):
    """A staged attachment was not created and still owned by this facade."""


class ServicesClosedError(RuntimeError):
    """The service facade has already released its storage resources."""


class TestLogServices:
    """Own the root lock, SQLite connection, and staged attachment copies."""

    def __init__(
        self,
        paths: DataPaths,
        lock: DataRootLock,
        repository: SQLiteRepository,
    ):
        self.paths = paths
        self._lock = lock
        self._repository = repository
        self._attachments = AttachmentStore(paths)
        self._backup_manager = BackupManager(paths)
        self._staged: dict[str, StagedAttachment] = {}
        self._condition = Condition()
        self._active_mutations = 0
        self._backup_active = False
        self._closed = False
        self._mutation_local = local()

    @classmethod
    def open(cls, paths: DataPaths) -> TestLogServices:
        lock = DataRootLock(paths.lock_file)
        repository: SQLiteRepository | None = None
        try:
            lock.acquire()
            Path(paths.root).mkdir(parents=True, exist_ok=True)
            Path(paths.attachments).mkdir(parents=True, exist_ok=True)
            repository = SQLiteRepository.open(paths.database)
            return cls(paths, lock, repository)
        except BaseException:
            if repository is not None:
                repository.close()
            lock.close()
            raise

    def close(self) -> None:
        with self._condition:
            while self._backup_active or self._active_mutations:
                self._condition.wait()
            if self._closed:
                return
            self._closed = True
        try:
            self._repository.close()
        finally:
            self._lock.close()

    @contextmanager
    def _mutation(self) -> Iterator[None]:
        depth = getattr(self._mutation_local, "depth", 0)
        if depth:
            self._mutation_local.depth = depth + 1
            try:
                yield
            finally:
                self._mutation_local.depth -= 1
            return

        with self._condition:
            self._ensure_open_locked()
            if self._backup_active:
                raise BackupInProgressError("Edits are paused while a backup is being written.")
            self._active_mutations += 1
        self._mutation_local.depth = 1
        try:
            yield
        finally:
            self._mutation_local.depth = 0
            with self._condition:
                self._active_mutations -= 1
                self._condition.notify_all()

    def _ensure_open_locked(self) -> None:
        if self._closed:
            raise ServicesClosedError("The test log service is closed.")

    def _ensure_open(self) -> None:
        with self._condition:
            self._ensure_open_locked()

    def _owned_staged(self, staged: Sequence[StagedAttachment]) -> tuple[StagedAttachment, ...]:
        values = tuple(staged)
        seen: set[str] = set()
        for attachment in values:
            if not isinstance(attachment, StagedAttachment):
                raise AttachmentOwnershipError("Pass staged attachment records created by this service.")
            current = self._staged.get(attachment.relative_path)
            if current != attachment:
                raise AttachmentOwnershipError("This service does not own that staged attachment copy.")
            if attachment.relative_path in seen:
                raise AttachmentOwnershipError("A staged attachment cannot be attached more than once.")
            seen.add(attachment.relative_path)
        return values

    def _discard_owned(self, staged: Sequence[StagedAttachment]) -> list[Path]:
        failures: list[Path] = []
        for attachment in dict.fromkeys(staged):
            if self._staged.get(attachment.relative_path) != attachment:
                continue
            current_failures = self._attachments.remove_staged((attachment,))
            failures.extend(current_failures)
            if not current_failures:
                self._staged.pop(attachment.relative_path, None)
        return failures

    @staticmethod
    def _add_cleanup_note(error: BaseException, failures: Sequence[Path]) -> None:
        if failures and hasattr(error, "add_note"):
            paths = ", ".join(str(path) for path in failures)
            error.add_note(f"Could not remove staged attachment copies: {paths}")

    def _save(
        self,
        method,
        record,
        staged: Sequence[StagedAttachment],
    ):
        with self._mutation():
            values = self._owned_staged(staged)
            try:
                for attachment in values:
                    self._attachments.resolve_staged(attachment)
                saved = method(record, attachments=values)
            except BaseException as error:
                self._add_cleanup_note(error, self._discard_owned(values))
                raise
            for attachment in values:
                self._staged.pop(attachment.relative_path, None)
            return saved

    def list_days(self) -> list[TestDay]:
        self._ensure_open()
        return self._repository.list_days()

    def list_setups(self, day_id: str) -> list[Setup]:
        self._ensure_open()
        return self._repository.list_setups(day_id)

    def list_laps(self, setup_id: str) -> list[Lap]:
        self._ensure_open()
        return self._repository.list_laps(setup_id)

    def list_event_layouts(self, include_archived: bool = False) -> list[EventLayout]:
        self._ensure_open()
        return self._repository.list_event_layouts(include_archived)

    def get_day(self, identifier: str) -> TestDay | None:
        self._ensure_open()
        return self._repository.get_day(identifier)

    def get_setup(self, identifier: str) -> Setup | None:
        self._ensure_open()
        return self._repository.get_setup(identifier)

    def get_lap(self, identifier: str) -> Lap | None:
        self._ensure_open()
        return self._repository.get_lap(identifier)

    def get_event_layout(self, identifier: str) -> EventLayout | None:
        self._ensure_open()
        return self._repository.get_event_layout(identifier)

    def list_attachments(self, owner_type: str, owner_id: str) -> list[Attachment]:
        self._ensure_open()
        return self._repository.list_attachments(owner_type, owner_id)

    def save_day(
        self,
        record: TestDay,
        staged: Sequence[StagedAttachment] = (),
    ) -> TestDay:
        return self._save(self._repository.save_day, record, staged)

    def save_setup(
        self,
        record: Setup,
        staged: Sequence[StagedAttachment] = (),
    ) -> Setup:
        return self._save(self._repository.save_setup, record, staged)

    def save_event_layout(
        self,
        record: EventLayout,
        staged: Sequence[StagedAttachment] = (),
    ) -> EventLayout:
        return self._save(self._repository.save_event_layout, record, staged)

    def save_lap(
        self,
        record: Lap,
        staged: Sequence[StagedAttachment] = (),
        *,
        lap_order: Sequence[str] | None = None,
    ) -> Lap:
        if lap_order is None:
            return self._save(self._repository.save_lap, record, staged)
        def save_with_order(lap, *, attachments):
            return self._repository.save_lap(lap, attachments, lap_order=lap_order)

        return self._save(save_with_order, record, staged)

    def stage_attachment(self, source: Path, role: AttachmentRole) -> StagedAttachment:
        with self._mutation():
            staged = self._attachments.stage(source, role)
            self._staged[staged.relative_path] = staged
            return staged

    def discard_staged(self, staged: Sequence[StagedAttachment]) -> list[Path]:
        with self._mutation():
            return self._discard_owned(tuple(staged))

    def resolve_attachment(self, attachment: Attachment) -> Path:
        self._ensure_open()
        return self._attachments.resolve(attachment)

    def _cleanup_deleted(self, attachments: Sequence[Attachment]) -> list[Path]:
        return self._attachments.remove_committed(attachments)

    def delete_day(self, identifier: str) -> list[Path]:
        with self._mutation():
            deleted = self._repository.delete_day(identifier)
            return self._cleanup_deleted(deleted)

    def delete_setup(self, identifier: str) -> list[Path]:
        with self._mutation():
            deleted = self._repository.delete_setup(identifier)
            return self._cleanup_deleted(deleted)

    def delete_lap(self, identifier: str) -> list[Path]:
        with self._mutation():
            deleted = self._repository.delete_lap(identifier)
            return self._cleanup_deleted(deleted)

    def event_is_referenced(self, identifier: str) -> bool:
        self._ensure_open()
        return self._repository.event_is_referenced(identifier)

    def delete_event_layout(self, identifier: str) -> list[Path]:
        with self._mutation():
            deleted = self._repository.delete_event_layout(identifier)
            return self._cleanup_deleted(deleted)

    def archive_event_layout(
        self,
        identifier: str,
        archived: bool = True,
    ) -> EventLayout:
        with self._mutation():
            return self._repository.archive_event_layout(identifier, archived)

    def reorder_laps(self, setup_id: str, lap_ids: Sequence[str]) -> list[Lap]:
        with self._mutation():
            return self._repository.reorder_laps(setup_id, lap_ids)

    def best_laps_by_event(self, setup_id: str) -> dict[str, Lap]:
        self._ensure_open()
        return self._repository.best_laps_by_event(setup_id)

    def duplicate_setup(self, identifier: str) -> Setup:
        with self._mutation():
            original = self._repository.get_setup(identifier)
            if original is None:
                raise ValidationError("setup_id", "The setup no longer exists.")
            copied_attachments: list[StagedAttachment] = []
            try:
                for attachment in self._repository.list_attachments("setup", identifier):
                    source = self._attachments.resolve(attachment)
                    copied = self._attachments.stage(source, attachment.role)
                    copied = replace(copied, original_name=attachment.original_name)
                    self._staged[copied.relative_path] = copied
                    copied_attachments.append(copied)
                setups = self._repository.list_setups(original.test_day_id)
                next_order = max((setup.order for setup in setups), default=0) + 1
                duplicate = replace(
                    original,
                    id=new_id(),
                    name=f"{original.name} (copy)",
                    setup_code=None,
                    order=next_order,
                    created_at=utc_now_iso(),
                )
                return self._save(self._repository.save_setup, duplicate, copied_attachments)
            except BaseException as error:
                self._add_cleanup_note(error, self._discard_owned(copied_attachments))
                raise

    def duplicate_event_layout(self, identifier: str) -> EventLayout:
        with self._mutation():
            original = self._repository.get_event_layout(identifier)
            if original is None:
                raise ValidationError("event_layout_id", "The event/layout no longer exists.")
            copied_attachments: list[StagedAttachment] = []
            try:
                for attachment in self._repository.list_attachments("event_layout", identifier):
                    source = self._attachments.resolve(attachment)
                    copied = self._attachments.stage(source, attachment.role)
                    copied = replace(copied, original_name=attachment.original_name)
                    self._staged[copied.relative_path] = copied
                    copied_attachments.append(copied)
                duplicate = replace(
                    original,
                    id=new_id(),
                    archived=False,
                    created_at=utc_now_iso(),
                )
                return self._save(
                    self._repository.save_event_layout,
                    duplicate,
                    copied_attachments,
                )
            except BaseException as error:
                self._add_cleanup_note(error, self._discard_owned(copied_attachments))
                raise

    def create_backup(self) -> Path:
        if getattr(self._mutation_local, "depth", 0):
            raise BackupInProgressError("A backup cannot start inside a write operation.")
        with self._condition:
            self._ensure_open_locked()
            if self._backup_active:
                raise BackupInProgressError("A backup is already being written.")
            self._backup_active = True
            while self._active_mutations:
                self._condition.wait()
        try:
            return self._backup_manager.create_backup(
                self._repository._connection,
                self._attachments,
            )
        finally:
            with self._condition:
                self._backup_active = False
                self._condition.notify_all()

    def export_setup(self, setup_id: str, destination: Path) -> Path:
        self._ensure_open()
        _guard_csv_destination(self.paths, destination)
        return write_setup_csv(self._repository, setup_id, destination)
