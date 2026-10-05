"""CSV export for one setup, preserving all lap statuses and notes."""

from __future__ import annotations

import csv
from pathlib import Path

from .repository import SQLiteRepository
from .time_value import format_lap_time
from .validation import ValidationError


CSV_COLUMNS = (
    "test_date",
    "test_location",
    "setup_name",
    "setup_code",
    "lap_sequence",
    "lap_time",
    "status",
    "driver",
    "time_of_day",
    "notes",
    "track_name",
    "event_name",
    "layout_name",
    "event_type",
    "length_m",
)


def export_setup(
    repository: SQLiteRepository,
    setup_id: str,
    destination: Path,
) -> Path:
    setup = repository.get_setup(setup_id)
    if setup is None:
        raise ValidationError("setup_id", "Select an existing setup to export.")
    day = repository.get_day(setup.test_day_id)
    if day is None:
        raise ValidationError("setup_id", "The setup's test day no longer exists.")

    rows: list[dict[str, object]] = []
    for lap in repository.list_laps(setup_id):
        event = repository.get_event_layout(lap.event_layout_id)
        if event is None:
            raise ValidationError("event_layout_id", "A lap's event/layout no longer exists.")
        rows.append(
            {
                "test_date": day.date,
                "test_location": day.location,
                "setup_name": setup.name,
                "setup_code": setup.setup_code,
                "lap_sequence": lap.sequence,
                "lap_time": format_lap_time(lap.time_ms),
                "status": lap.status,
                "driver": lap.driver,
                "time_of_day": lap.time_of_day,
                "notes": lap.notes,
                "track_name": event.track_name,
                "event_name": event.event_name,
                "layout_name": event.layout_name,
                "event_type": event.event_type,
                "length_m": event.length_m,
            }
        )

    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    return output
