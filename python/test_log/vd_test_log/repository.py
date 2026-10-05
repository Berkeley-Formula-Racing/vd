"""Transactional SQLite storage for offline test log records."""

from __future__ import annotations

import sqlite3
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from .models import (
    Attachment,
    EventLayout,
    Lap,
    Setup,
    StagedAttachment,
    TestDay,
    new_id,
    utc_now_iso,
)
from .validation import (
    ValidationError,
    validate_day,
    validate_event_layout,
    validate_lap,
    validate_setup,
    validate_staged_attachment,
)


SCHEMA_VERSION = 1

_OWNER_COLUMNS = {
    "day": "day_owner_id",
    "setup": "setup_owner_id",
    "lap": "lap_owner_id",
    "event_layout": "event_layout_owner_id",
}

_CREATE_STATEMENTS = (
    """
    CREATE TABLE test_days (
        id TEXT PRIMARY KEY,
        date TEXT NOT NULL,
        location TEXT NOT NULL,
        weather TEXT,
        notes TEXT,
        created_at TEXT NOT NULL
    )
    """,
    """
    CREATE TABLE setups (
        id TEXT PRIMARY KEY,
        test_day_id TEXT NOT NULL REFERENCES test_days(id) ON DELETE CASCADE,
        name TEXT NOT NULL,
        setup_code TEXT,
        settings_text TEXT NOT NULL,
        notes TEXT,
        sort_order INTEGER NOT NULL CHECK (sort_order > 0),
        created_at TEXT NOT NULL,
        UNIQUE (test_day_id, sort_order)
    )
    """,
    """
    CREATE TABLE event_layouts (
        id TEXT PRIMARY KEY,
        track_name TEXT NOT NULL,
        layout_name TEXT NOT NULL,
        event_name TEXT,
        event_type TEXT,
        length_m REAL CHECK (length_m IS NULL OR length_m > 0),
        notes TEXT,
        archived INTEGER NOT NULL CHECK (archived IN (0, 1)),
        created_at TEXT NOT NULL
    )
    """,
    """
    CREATE TABLE laps (
        id TEXT PRIMARY KEY,
        setup_id TEXT NOT NULL REFERENCES setups(id) ON DELETE CASCADE,
        event_layout_id TEXT NOT NULL REFERENCES event_layouts(id) ON DELETE RESTRICT,
        sequence INTEGER NOT NULL CHECK (sequence > 0),
        time_ms INTEGER NOT NULL CHECK (time_ms > 0),
        status TEXT NOT NULL CHECK (status IN ('valid', 'invalid')),
        driver TEXT,
        time_of_day TEXT,
        notes TEXT,
        created_at TEXT NOT NULL,
        UNIQUE (setup_id, sequence)
    )
    """,
    """
    CREATE TABLE attachments (
        id TEXT PRIMARY KEY,
        owner_type TEXT NOT NULL CHECK (owner_type IN ('day', 'setup', 'lap', 'event_layout')),
        day_owner_id TEXT REFERENCES test_days(id) ON DELETE CASCADE,
        setup_owner_id TEXT REFERENCES setups(id) ON DELETE CASCADE,
        lap_owner_id TEXT REFERENCES laps(id) ON DELETE CASCADE,
        event_layout_owner_id TEXT REFERENCES event_layouts(id) ON DELETE CASCADE,
        role TEXT NOT NULL CHECK (role IN ('file', 'map')),
        original_name TEXT NOT NULL,
        relative_path TEXT NOT NULL UNIQUE,
        created_at TEXT NOT NULL,
        CHECK (
            (owner_type = 'day' AND day_owner_id IS NOT NULL
                AND setup_owner_id IS NULL AND lap_owner_id IS NULL
                AND event_layout_owner_id IS NULL)
            OR
            (owner_type = 'setup' AND day_owner_id IS NULL
                AND setup_owner_id IS NOT NULL AND lap_owner_id IS NULL
                AND event_layout_owner_id IS NULL)
            OR
            (owner_type = 'lap' AND day_owner_id IS NULL
                AND setup_owner_id IS NULL AND lap_owner_id IS NOT NULL
                AND event_layout_owner_id IS NULL)
            OR
            (owner_type = 'event_layout' AND day_owner_id IS NULL
                AND setup_owner_id IS NULL AND lap_owner_id IS NULL
                AND event_layout_owner_id IS NOT NULL)
        ),
        CHECK (role = 'file' OR (role = 'map' AND owner_type = 'event_layout'))
    )
    """,
    "CREATE INDEX setups_by_day_order ON setups(test_day_id, sort_order)",
    "CREATE INDEX laps_by_setup_sequence ON laps(setup_id, sequence)",
    "CREATE INDEX laps_by_event_layout ON laps(event_layout_id)",
    "CREATE INDEX attachments_by_day ON attachments(day_owner_id)",
    "CREATE INDEX attachments_by_setup ON attachments(setup_owner_id)",
    "CREATE INDEX attachments_by_lap ON attachments(lap_owner_id)",
    "CREATE INDEX attachments_by_event_layout ON attachments(event_layout_owner_id)",
    """
    CREATE TRIGGER referenced_event_fields_immutable
    BEFORE UPDATE OF track_name, layout_name, event_name, event_type, length_m, notes
    ON event_layouts
    WHEN EXISTS (SELECT 1 FROM laps WHERE event_layout_id = OLD.id)
      AND (
          OLD.track_name IS NOT NEW.track_name
          OR OLD.layout_name IS NOT NEW.layout_name
          OR OLD.event_name IS NOT NEW.event_name
          OR OLD.event_type IS NOT NEW.event_type
          OR OLD.length_m IS NOT NEW.length_m
          OR OLD.notes IS NOT NEW.notes
      )
    BEGIN
        SELECT RAISE(ABORT, 'referenced event/layout details are immutable');
    END
    """,
    """
    CREATE TRIGGER referenced_event_attachment_insert_immutable
    BEFORE INSERT ON attachments
    WHEN NEW.event_layout_owner_id IS NOT NULL
      AND EXISTS (SELECT 1 FROM laps WHERE event_layout_id = NEW.event_layout_owner_id)
    BEGIN
        SELECT RAISE(ABORT, 'referenced event/layout attachments are immutable');
    END
    """,
    """
    CREATE TRIGGER referenced_event_attachment_update_immutable
    BEFORE UPDATE ON attachments
    WHEN (
        OLD.event_layout_owner_id IS NOT NULL
        AND EXISTS (SELECT 1 FROM laps WHERE event_layout_id = OLD.event_layout_owner_id)
    ) OR (
        NEW.event_layout_owner_id IS NOT NULL
        AND EXISTS (SELECT 1 FROM laps WHERE event_layout_id = NEW.event_layout_owner_id)
    )
    BEGIN
        SELECT RAISE(ABORT, 'referenced event/layout attachments are immutable');
    END
    """,
    """
    CREATE TRIGGER referenced_event_attachment_delete_immutable
    BEFORE DELETE ON attachments
    WHEN OLD.event_layout_owner_id IS NOT NULL
      AND EXISTS (SELECT 1 FROM laps WHERE event_layout_id = OLD.event_layout_owner_id)
    BEGIN
        SELECT RAISE(ABORT, 'referenced event/layout attachments are immutable');
    END
    """,
)

class SQLiteRepository:
    """Owns one SQLite connection and commits each write operation atomically."""

    def __init__(self, connection: sqlite3.Connection):
        self._connection = connection
        self._closed = False

    @classmethod
    def open(cls, database_path: Path) -> SQLiteRepository:
        path = Path(database_path)
        existed_before_open = path.exists()
        original_size = path.stat().st_size if existed_before_open else 0
        connection = sqlite3.connect(path)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        try:
            enabled = connection.execute("PRAGMA foreign_keys").fetchone()[0]
            if enabled != 1:
                raise ValidationError("database", "SQLite foreign-key enforcement is unavailable.")
            cls._initialize_or_migrate(connection, existed_before_open, original_size)
        except Exception:
            connection.close()
            raise
        return cls(connection)

    @staticmethod
    def _initialize_or_migrate(
        connection: sqlite3.Connection,
        existed_before_open: bool,
        original_size: int,
    ) -> None:
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        table_names = {
            row["name"]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        if version == SCHEMA_VERSION:
            required = {"test_days", "setups", "event_layouts", "laps", "attachments"}
            if not required.issubset(table_names):
                raise ValidationError(
                    "database",
                    "The version 1 test log database is incomplete and was not changed.",
                )
            return
        if version > SCHEMA_VERSION:
            raise ValidationError(
                "database",
                f"Database schema version {version} is newer than this app supports.",
            )
        if version != 0 or table_names or (existed_before_open and original_size > 0):
            raise ValidationError(
                "database",
                f"Database schema version {version} is unknown; the existing file was not changed.",
            )

        connection.execute("BEGIN IMMEDIATE")
        try:
            for statement in _CREATE_STATEMENTS:
                connection.execute(statement)
            connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
            connection.commit()
        except Exception:
            connection.rollback()
            raise

    def close(self) -> None:
        if not self._closed:
            self._connection.close()
            self._closed = True

    @staticmethod
    def _day_from_row(row: sqlite3.Row | None) -> TestDay | None:
        if row is None:
            return None
        return TestDay(
            id=row["id"],
            date=row["date"],
            location=row["location"],
            weather=row["weather"],
            notes=row["notes"],
            created_at=row["created_at"],
        )

    @staticmethod
    def _setup_from_row(row: sqlite3.Row | None) -> Setup | None:
        if row is None:
            return None
        return Setup(
            id=row["id"],
            test_day_id=row["test_day_id"],
            name=row["name"],
            setup_code=row["setup_code"],
            settings_text=row["settings_text"],
            notes=row["notes"],
            order=row["sort_order"],
            created_at=row["created_at"],
        )

    @staticmethod
    def _event_from_row(row: sqlite3.Row | None) -> EventLayout | None:
        if row is None:
            return None
        length = row["length_m"]
        return EventLayout(
            id=row["id"],
            track_name=row["track_name"],
            layout_name=row["layout_name"],
            event_name=row["event_name"],
            event_type=row["event_type"],
            length_m=float(length) if length is not None else None,
            notes=row["notes"],
            archived=bool(row["archived"]),
            created_at=row["created_at"],
        )

    @staticmethod
    def _lap_from_row(row: sqlite3.Row | None) -> Lap | None:
        if row is None:
            return None
        return Lap(
            id=row["id"],
            setup_id=row["setup_id"],
            event_layout_id=row["event_layout_id"],
            sequence=row["sequence"],
            time_ms=row["time_ms"],
            status=row["status"],
            driver=row["driver"],
            time_of_day=row["time_of_day"],
            notes=row["notes"],
            created_at=row["created_at"],
        )

    @staticmethod
    def _attachment_from_row(row: sqlite3.Row) -> Attachment:
        return Attachment(
            id=row["id"],
            owner_type=row["owner_type"],
            owner_id=row["owner_id"],
            role=row["role"],
            original_name=row["original_name"],
            relative_path=row["relative_path"],
            created_at=row["created_at"],
        )

    def list_days(self) -> list[TestDay]:
        rows = self._connection.execute(
            "SELECT * FROM test_days ORDER BY date DESC, created_at DESC, id"
        )
        return [self._day_from_row(row) for row in rows]

    def list_setups(self, day_id: str) -> list[Setup]:
        rows = self._connection.execute(
            "SELECT * FROM setups WHERE test_day_id = ? ORDER BY sort_order, id",
            (day_id,),
        )
        return [self._setup_from_row(row) for row in rows]

    def list_laps(self, setup_id: str) -> list[Lap]:
        rows = self._connection.execute(
            "SELECT * FROM laps WHERE setup_id = ? ORDER BY sequence, id",
            (setup_id,),
        )
        return [self._lap_from_row(row) for row in rows]

    def list_event_layouts(self, include_archived: bool = False) -> list[EventLayout]:
        if include_archived:
            rows = self._connection.execute(
                """
                SELECT * FROM event_layouts
                ORDER BY track_name, layout_name, event_name, created_at, id
                """
            )
        else:
            rows = self._connection.execute(
                """
                SELECT * FROM event_layouts WHERE archived = 0
                ORDER BY track_name, layout_name, event_name, created_at, id
                """
            )
        return [self._event_from_row(row) for row in rows]

    def get_day(self, identifier: str) -> TestDay | None:
        row = self._connection.execute(
            "SELECT * FROM test_days WHERE id = ?", (identifier,)
        ).fetchone()
        return self._day_from_row(row)

    def get_setup(self, identifier: str) -> Setup | None:
        row = self._connection.execute(
            "SELECT * FROM setups WHERE id = ?", (identifier,)
        ).fetchone()
        return self._setup_from_row(row)

    def get_lap(self, identifier: str) -> Lap | None:
        row = self._connection.execute(
            "SELECT * FROM laps WHERE id = ?", (identifier,)
        ).fetchone()
        return self._lap_from_row(row)

    def get_event_layout(self, identifier: str) -> EventLayout | None:
        row = self._connection.execute(
            "SELECT * FROM event_layouts WHERE id = ?", (identifier,)
        ).fetchone()
        return self._event_from_row(row)

    def list_attachments(self, owner_type: str, owner_id: str) -> list[Attachment]:
        if owner_type not in _OWNER_COLUMNS:
            raise ValidationError("owner_type", "Choose a supported attachment owner.")
        if not isinstance(owner_id, str) or not owner_id.strip():
            raise ValidationError("owner_id", "An attachment owner ID is required.")
        owner_column = _OWNER_COLUMNS[owner_type]
        rows = self._connection.execute(
            f"""
            SELECT id, owner_type,
                   COALESCE(day_owner_id, setup_owner_id, lap_owner_id,
                            event_layout_owner_id) AS owner_id,
                   role, original_name, relative_path, created_at
            FROM attachments
            WHERE {owner_column} = ?
            ORDER BY created_at, id
            """,
            (owner_id,),
        )
        return [self._attachment_from_row(row) for row in rows]

    def save_day(
        self,
        record: TestDay,
        attachments: Sequence[StagedAttachment] = (),
    ) -> TestDay:
        staged = self._prepare_attachments(attachments, "day")
        with self._connection:
            previous = self.get_day(record.id)
            saved = replace(record, created_at=previous.created_at) if previous else record
            validate_day(saved)
            if previous is None:
                self._connection.execute(
                    """
                    INSERT INTO test_days(id, date, location, weather, notes, created_at)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (saved.id, saved.date, saved.location, saved.weather, saved.notes, saved.created_at),
                )
            else:
                self._connection.execute(
                    """
                    UPDATE test_days
                    SET date = ?, location = ?, weather = ?, notes = ?
                    WHERE id = ?
                    """,
                    (saved.date, saved.location, saved.weather, saved.notes, saved.id),
                )
            self._insert_attachments("day", saved.id, staged)
        return self.get_day(saved.id)

    def save_setup(
        self,
        record: Setup,
        attachments: Sequence[StagedAttachment] = (),
    ) -> Setup:
        staged = self._prepare_attachments(attachments, "setup")
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            previous = self.get_setup(record.id)
            if previous is not None and previous.test_day_id != record.test_day_id:
                raise ValidationError("test_day_id", "A saved setup cannot be moved to another test day.")
            saved = replace(record, created_at=previous.created_at) if previous else record
            validate_setup(saved)
            if self.get_day(saved.test_day_id) is None:
                raise ValidationError("test_day_id", "Select an existing test day.")
            rows = list(
                self._connection.execute(
                    "SELECT id, sort_order FROM setups WHERE test_day_id = ? ORDER BY sort_order, id",
                    (saved.test_day_id,),
                )
            )
            ordered_ids = [row["id"] for row in rows if row["id"] != saved.id]
            insertion_index = min(saved.order - 1, len(ordered_ids))
            ordered_ids.insert(insertion_index, saved.id)
            max_order = max((row["sort_order"] for row in rows), default=0)
            offset = max_order + len(rows) + 1
            if rows:
                self._connection.execute(
                    "UPDATE setups SET sort_order = sort_order + ? WHERE test_day_id = ?",
                    (offset, saved.test_day_id),
                )
            if previous is None:
                temporary_order = max_order + offset + len(rows) + 1
                self._connection.execute(
                    """
                    INSERT INTO setups(
                        id, test_day_id, name, setup_code, settings_text, notes,
                        sort_order, created_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        saved.id,
                        saved.test_day_id,
                        saved.name,
                        saved.setup_code,
                        saved.settings_text,
                        saved.notes,
                        temporary_order,
                        saved.created_at,
                    ),
                )
            else:
                self._connection.execute(
                    """
                    UPDATE setups
                    SET name = ?, setup_code = ?, settings_text = ?, notes = ?
                    WHERE id = ?
                    """,
                    (saved.name, saved.setup_code, saved.settings_text, saved.notes, saved.id),
                )
            for order_value, setup_id in enumerate(ordered_ids, start=1):
                self._connection.execute(
                    "UPDATE setups SET sort_order = ? WHERE id = ?",
                    (order_value, setup_id),
                )
            self._insert_attachments("setup", saved.id, staged)
        return self.get_setup(saved.id)

    def _prepare_attachments(
        self,
        attachments: Sequence[StagedAttachment],
        owner_type: str,
    ) -> tuple[StagedAttachment, ...]:
        staged = tuple(attachments)
        for attachment in staged:
            if not isinstance(attachment, StagedAttachment):
                raise ValidationError(
                    "attachments",
                    "Pass staged attachment records when saving attachments.",
                )
            validate_staged_attachment(attachment, owner_type)
        return staged

    def _insert_attachments(
        self,
        owner_type: str,
        owner_id: str,
        attachments: Sequence[StagedAttachment],
    ) -> None:
        owner_column = _OWNER_COLUMNS[owner_type]
        for attachment in attachments:
            try:
                self._connection.execute(
                    f"""
                    INSERT INTO attachments(
                        id, owner_type, {owner_column}, role, original_name,
                        relative_path, created_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        new_id(),
                        owner_type,
                        owner_id,
                        attachment.role,
                        attachment.original_name,
                        attachment.relative_path,
                        utc_now_iso(),
                    ),
                )
            except sqlite3.IntegrityError as error:
                raise ValidationError(
                    "attachments",
                    "An attachment could not be saved for this record.",
                ) from error

    def save_event_layout(
        self,
        record: EventLayout,
        attachments: Sequence[StagedAttachment] = (),
    ) -> EventLayout:
        staged = self._prepare_attachments(attachments, "event_layout")
        with self._connection:
            previous = self.get_event_layout(record.id)
            saved = replace(record, created_at=previous.created_at) if previous else record
            validate_event_layout(saved)
            referenced = self.event_is_referenced(saved.id)
            if referenced and previous is not None:
                immutable_fields = (
                    "track_name",
                    "layout_name",
                    "event_name",
                    "event_type",
                    "length_m",
                    "notes",
                )
                if any(getattr(previous, field) != getattr(saved, field) for field in immutable_fields):
                    raise ValidationError(
                        "event_layout",
                        "Referenced event/layout details are read-only; duplicate the record to change them.",
                    )
            if referenced and staged:
                raise ValidationError(
                    "attachments",
                    "Attachments on a referenced event/layout are read-only.",
                )
            if previous is None:
                self._connection.execute(
                    """
                    INSERT INTO event_layouts(
                        id, track_name, layout_name, event_name, event_type,
                        length_m, notes, archived, created_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        saved.id,
                        saved.track_name,
                        saved.layout_name,
                        saved.event_name,
                        saved.event_type,
                        saved.length_m,
                        saved.notes,
                        int(saved.archived),
                        saved.created_at,
                    ),
                )
            else:
                self._connection.execute(
                    """
                    UPDATE event_layouts
                    SET track_name = ?, layout_name = ?, event_name = ?, event_type = ?,
                        length_m = ?, notes = ?, archived = ?
                    WHERE id = ?
                    """,
                    (
                        saved.track_name,
                        saved.layout_name,
                        saved.event_name,
                        saved.event_type,
                        saved.length_m,
                        saved.notes,
                        int(saved.archived),
                        saved.id,
                    ),
                )
            self._insert_attachments("event_layout", saved.id, staged)
        return self.get_event_layout(saved.id)

    def save_lap(
        self,
        record: Lap,
        attachments: Sequence[StagedAttachment] = (),
        *,
        lap_order: Sequence[str] | None = None,
    ) -> Lap:
        staged = self._prepare_attachments(attachments, "lap")
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            previous = self.get_lap(record.id)
            if previous is not None and previous.setup_id != record.setup_id:
                raise ValidationError("setup_id", "A saved lap cannot be moved to another setup.")
            saved = replace(record, created_at=previous.created_at) if previous else record
            validate_lap(saved)
            if self.get_setup(saved.setup_id) is None:
                raise ValidationError("setup_id", "Select an existing setup.")
            event = self.get_event_layout(saved.event_layout_id)
            if event is None:
                raise ValidationError("event_layout_id", "Select an existing event/layout.")
            if event.archived and (
                previous is None or previous.event_layout_id != saved.event_layout_id
            ):
                raise ValidationError(
                    "event_layout_id",
                    "Archived event/layouts can remain on historical laps but cannot be assigned to a lap.",
                )
            rows = list(
                self._connection.execute(
                    "SELECT id, sequence FROM laps WHERE setup_id = ? ORDER BY sequence, id",
                    (saved.setup_id,),
                )
            )
            ordered_ids = [row["id"] for row in rows if row["id"] != saved.id]
            insertion_index = min(saved.sequence - 1, len(ordered_ids))
            ordered_ids.insert(insertion_index, saved.id)
            max_sequence = max((row["sequence"] for row in rows), default=0)
            offset = max_sequence + len(rows) + 1
            if rows:
                self._connection.execute(
                    "UPDATE laps SET sequence = sequence + ? WHERE setup_id = ?",
                    (offset, saved.setup_id),
                )
            if previous is None:
                temporary_sequence = max_sequence + offset + len(rows) + 1
                self._connection.execute(
                    """
                    INSERT INTO laps(
                        id, setup_id, event_layout_id, sequence, time_ms, status,
                        driver, time_of_day, notes, created_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        saved.id,
                        saved.setup_id,
                        saved.event_layout_id,
                        temporary_sequence,
                        saved.time_ms,
                        saved.status,
                        saved.driver,
                        saved.time_of_day,
                        saved.notes,
                        saved.created_at,
                    ),
                )
            else:
                self._connection.execute(
                    """
                    UPDATE laps
                    SET event_layout_id = ?, time_ms = ?, status = ?, driver = ?,
                        time_of_day = ?, notes = ?
                    WHERE id = ?
                    """,
                    (
                        saved.event_layout_id,
                        saved.time_ms,
                        saved.status,
                        saved.driver,
                        saved.time_of_day,
                        saved.notes,
                        saved.id,
                    ),
                )
            for sequence, lap_id in enumerate(ordered_ids, start=1):
                self._connection.execute(
                    "UPDATE laps SET sequence = ? WHERE id = ?",
                    (sequence, lap_id),
                )
            self._insert_attachments("lap", saved.id, staged)
            if lap_order is not None:
                self._reorder_laps_in_transaction(saved.setup_id, lap_order)
        return self.get_lap(saved.id)

    def delete_day(self, identifier: str) -> list[Attachment]:
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            attachments = self._attachments_for_query(
                """
                day_owner_id = ?
                OR setup_owner_id IN (
                    SELECT id FROM setups WHERE test_day_id = ?
                )
                OR lap_owner_id IN (
                    SELECT laps.id FROM laps
                    JOIN setups ON setups.id = laps.setup_id
                    WHERE setups.test_day_id = ?
                )
                """,
                (identifier, identifier, identifier),
            )
            self._connection.execute("DELETE FROM test_days WHERE id = ?", (identifier,))
        return attachments

    def delete_setup(self, identifier: str) -> list[Attachment]:
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            setup = self.get_setup(identifier)
            if setup is None:
                return []
            attachments = self._attachments_for_query(
                """
                setup_owner_id = ?
                OR lap_owner_id IN (SELECT id FROM laps WHERE setup_id = ?)
                """,
                (identifier, identifier),
            )
            self._connection.execute("DELETE FROM setups WHERE id = ?", (identifier,))
            self._renumber_setups(setup.test_day_id)
        return attachments

    def delete_lap(self, identifier: str) -> list[Attachment]:
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            lap = self.get_lap(identifier)
            if lap is None:
                return []
            attachments = self._attachments_for_query(
                "lap_owner_id = ?",
                (identifier,),
            )
            self._connection.execute("DELETE FROM laps WHERE id = ?", (identifier,))
            self._renumber_laps(lap.setup_id)
        return attachments

    def event_is_referenced(self, identifier: str) -> bool:
        row = self._connection.execute(
            "SELECT 1 FROM laps WHERE event_layout_id = ? LIMIT 1",
            (identifier,),
        ).fetchone()
        return row is not None

    def delete_event_layout(self, identifier: str) -> list[Attachment]:
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            event = self.get_event_layout(identifier)
            if event is None:
                return []
            if self.event_is_referenced(identifier):
                raise ValidationError(
                    "event_layout",
                    "A referenced event/layout cannot be deleted; archive it instead.",
                )
            attachments = self._attachments_for_query(
                "event_layout_owner_id = ?",
                (identifier,),
            )
            self._connection.execute(
                "DELETE FROM event_layouts WHERE id = ?",
                (identifier,),
            )
        return attachments

    def archive_event_layout(
        self,
        identifier: str,
        archived: bool = True,
    ) -> EventLayout:
        if not isinstance(archived, bool):
            raise ValidationError("archived", "Archived must be true or false.")
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            record = self.get_event_layout(identifier)
            if record is None:
                raise ValidationError("event_layout", "The event/layout no longer exists.")
            self._connection.execute(
                "UPDATE event_layouts SET archived = ? WHERE id = ?",
                (int(archived), identifier),
            )
        return self.get_event_layout(identifier)

    def _reorder_laps_in_transaction(self, setup_id: str, lap_ids: Sequence[str]) -> list[Lap]:
        if isinstance(lap_ids, (str, bytes)):
            raise ValidationError("lap_ids", "Provide the lap IDs in the desired order.")
        order = tuple(lap_ids)
        if self.get_setup(setup_id) is None:
            raise ValidationError("setup_id", "Select an existing setup.")
        rows = list(
            self._connection.execute(
                "SELECT id, sequence FROM laps WHERE setup_id = ? ORDER BY sequence, id",
                (setup_id,),
            )
        )
        existing_ids = [row["id"] for row in rows]
        if len(order) != len(existing_ids) or set(order) != set(existing_ids):
            raise ValidationError(
                "lap_ids",
                "Lap order must contain every lap in this setup exactly once.",
            )
        if not order:
            return []
        offset = max(row["sequence"] for row in rows) + len(rows) + 1
        self._connection.execute(
            "UPDATE laps SET sequence = sequence + ? WHERE setup_id = ?",
            (offset, setup_id),
        )
        for sequence, lap_id in enumerate(order, start=1):
            self._connection.execute(
                "UPDATE laps SET sequence = ? WHERE id = ? AND setup_id = ?",
                (sequence, lap_id, setup_id),
            )
        return self.list_laps(setup_id)

    def reorder_laps(self, setup_id: str, lap_ids: Sequence[str]) -> list[Lap]:
        with self._connection:
            self._connection.execute("BEGIN IMMEDIATE")
            return self._reorder_laps_in_transaction(setup_id, lap_ids)

    def best_laps_by_event(self, setup_id: str) -> dict[str, Lap]:
        rows = self._connection.execute(
            """
            SELECT * FROM laps
            WHERE setup_id = ? AND status = 'valid'
            ORDER BY event_layout_id, time_ms, sequence, id
            """,
            (setup_id,),
        )
        best: dict[str, Lap] = {}
        for row in rows:
            event_id = row["event_layout_id"]
            if event_id not in best:
                best[event_id] = self._lap_from_row(row)
        return best

    def _attachments_for_query(
        self,
        predicate: str,
        parameters: tuple[object, ...],
    ) -> list[Attachment]:
        rows = self._connection.execute(
            f"""
            SELECT id, owner_type,
                   COALESCE(day_owner_id, setup_owner_id, lap_owner_id,
                            event_layout_owner_id) AS owner_id,
                   role, original_name, relative_path, created_at
            FROM attachments
            WHERE {predicate}
            ORDER BY created_at, id
            """,
            parameters,
        )
        return [self._attachment_from_row(row) for row in rows]

    def _renumber_setups(self, day_id: str) -> None:
        rows = list(
            self._connection.execute(
                "SELECT id, sort_order FROM setups WHERE test_day_id = ? ORDER BY sort_order, id",
                (day_id,),
            )
        )
        if not rows:
            return
        offset = max(row["sort_order"] for row in rows) + 1
        self._connection.execute(
            "UPDATE setups SET sort_order = sort_order + ? WHERE test_day_id = ?",
            (offset, day_id),
        )
        for order_value, row in enumerate(rows, start=1):
            self._connection.execute(
                "UPDATE setups SET sort_order = ? WHERE id = ?",
                (order_value, row["id"]),
            )

    def _renumber_laps(self, setup_id: str) -> None:
        rows = list(
            self._connection.execute(
                "SELECT id, sequence FROM laps WHERE setup_id = ? ORDER BY sequence, id",
                (setup_id,),
            )
        )
        if not rows:
            return
        offset = max(row["sequence"] for row in rows) + 1
        self._connection.execute(
            "UPDATE laps SET sequence = sequence + ? WHERE setup_id = ?",
            (offset, setup_id),
        )
        for sequence, row in enumerate(rows, start=1):
            self._connection.execute(
                "UPDATE laps SET sequence = ? WHERE id = ?",
                (sequence, row["id"]),
            )
