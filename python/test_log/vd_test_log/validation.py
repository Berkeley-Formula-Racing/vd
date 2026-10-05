"""Field validation shared by input parsing and SQLite persistence."""

from __future__ import annotations

import math
import re
from datetime import date, datetime
from typing import Any

from .models import EventLayout, Lap, Setup, StagedAttachment, TestDay


class ValidationError(ValueError):
    """A user-correctable input error with its associated field name."""

    def __init__(self, field: str, message: str):
        self.field = field
        self.message = message
        super().__init__(message)


_DATE_PATTERN = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}\Z")


def _required_text(value: Any, field: str, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(field, f"Enter a {label}.")


def _optional_text(value: Any, field: str) -> None:
    if value is not None and not isinstance(value, str):
        raise ValidationError(field, "Enter text or leave this field blank.")


def _aware_timestamp(value: Any, field: str = "created_at") -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(field, "A creation timestamp is required.")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise ValidationError(field, "Use an ISO 8601 timestamp.") from error
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValidationError(field, "Use a timezone-aware ISO 8601 timestamp.")


def validate_day(record: TestDay) -> None:
    _required_text(record.id, "id", "record ID")
    if not isinstance(record.date, str) or not _DATE_PATTERN.fullmatch(record.date):
        raise ValidationError("date", "Enter a date in YYYY-MM-DD format.")
    try:
        parsed_date = date.fromisoformat(record.date)
    except ValueError as error:
        raise ValidationError("date", "Enter a real calendar date.") from error
    if parsed_date.isoformat() != record.date:
        raise ValidationError("date", "Enter a real calendar date.")
    _required_text(record.location, "location", "location")
    _optional_text(record.weather, "weather")
    _optional_text(record.notes, "notes")
    _aware_timestamp(record.created_at)


def validate_setup(record: Setup) -> None:
    _required_text(record.id, "id", "record ID")
    _required_text(record.test_day_id, "test_day_id", "test day")
    _required_text(record.name, "name", "setup name")
    _optional_text(record.setup_code, "setup_code")
    if not isinstance(record.settings_text, str):
        raise ValidationError("settings_text", "Enter setup settings as text.")
    _optional_text(record.notes, "notes")
    if isinstance(record.order, bool) or not isinstance(record.order, int) or record.order <= 0:
        raise ValidationError("order", "Setup order must be a positive whole number.")
    _aware_timestamp(record.created_at)


def validate_event_layout(record: EventLayout) -> None:
    _required_text(record.id, "id", "record ID")
    _required_text(record.track_name, "track_name", "track or venue name")
    _required_text(record.layout_name, "layout_name", "layout name")
    _optional_text(record.event_name, "event_name")
    _optional_text(record.event_type, "event_type")
    _optional_text(record.notes, "notes")
    if record.length_m is not None:
        if isinstance(record.length_m, bool) or not isinstance(record.length_m, (int, float)):
            raise ValidationError("length_m", "Length must be a finite positive number of metres.")
        if not math.isfinite(record.length_m) or record.length_m <= 0:
            raise ValidationError("length_m", "Length must be a finite positive number of metres.")
    if not isinstance(record.archived, bool):
        raise ValidationError("archived", "Archived must be true or false.")
    _aware_timestamp(record.created_at)


def validate_lap(record: Lap) -> None:
    _required_text(record.id, "id", "record ID")
    _required_text(record.setup_id, "setup_id", "setup")
    _required_text(record.event_layout_id, "event_layout_id", "event and layout")
    if isinstance(record.sequence, bool) or not isinstance(record.sequence, int) or record.sequence <= 0:
        raise ValidationError("sequence", "Lap sequence must be a positive whole number.")
    if isinstance(record.time_ms, bool) or not isinstance(record.time_ms, int) or record.time_ms <= 0:
        raise ValidationError("time_ms", "Lap time must be a positive number of milliseconds.")
    if record.status not in ("valid", "invalid"):
        raise ValidationError("status", "Choose valid or invalid lap status.")
    _optional_text(record.driver, "driver")
    _optional_text(record.time_of_day, "time_of_day")
    _optional_text(record.notes, "notes")
    _aware_timestamp(record.created_at)


def validate_staged_attachment(attachment: StagedAttachment, owner_type: str) -> None:
    if attachment.role not in ("file", "map"):
        raise ValidationError("role", "Attachment role must be file or map.")
    if attachment.role == "map" and owner_type != "event_layout":
        raise ValidationError("role", "Map attachments belong to event/layout records.")
    _required_text(attachment.original_name, "original_name", "original filename")
    _required_text(attachment.relative_path, "relative_path", "stored relative path")
