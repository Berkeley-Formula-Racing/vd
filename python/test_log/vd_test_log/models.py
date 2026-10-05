"""Immutable records shared by the repository, services, and UI."""

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal
from uuid import uuid4


OwnerType = Literal["day", "setup", "lap", "event_layout"]
AttachmentRole = Literal["file", "map"]
LapStatus = Literal["valid", "invalid"]


def new_id() -> str:
    """Return a stable, compact UUID suitable for a record primary key."""
    return uuid4().hex


def utc_now_iso() -> str:
    """Return the current timezone-aware UTC timestamp in ISO 8601 form."""
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


@dataclass(frozen=True, slots=True)
class TestDay:
    id: str
    date: str
    location: str
    weather: str | None
    notes: str | None
    created_at: str


@dataclass(frozen=True, slots=True)
class Setup:
    id: str
    test_day_id: str
    name: str
    setup_code: str | None
    settings_text: str
    notes: str | None
    order: int
    created_at: str


@dataclass(frozen=True, slots=True)
class EventLayout:
    id: str
    track_name: str
    layout_name: str
    event_name: str | None
    event_type: str | None
    length_m: float | None
    notes: str | None
    archived: bool
    created_at: str


@dataclass(frozen=True, slots=True)
class Lap:
    id: str
    setup_id: str
    event_layout_id: str
    sequence: int
    time_ms: int
    status: LapStatus
    driver: str | None
    time_of_day: str | None
    notes: str | None
    created_at: str


@dataclass(frozen=True, slots=True)
class Attachment:
    id: str
    owner_type: OwnerType
    owner_id: str
    role: AttachmentRole
    original_name: str
    relative_path: str
    created_at: str


@dataclass(frozen=True, slots=True)
class StagedAttachment:
    role: AttachmentRole
    original_name: str
    relative_path: str


@dataclass(frozen=True, slots=True)
class DataPaths:
    root: Path
    database: Path
    attachments: Path
    lock_file: Path
    backup_root: Path
