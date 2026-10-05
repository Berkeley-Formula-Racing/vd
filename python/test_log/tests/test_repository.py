import dataclasses
from contextlib import closing
import math
import sqlite3
import tempfile
from threading import Thread
import unittest
from datetime import datetime
from pathlib import Path

from vd_test_log.models import (
    EventLayout,
    Lap,
    Setup,
    StagedAttachment,
    TestDay,
)
from vd_test_log.repository import SQLiteRepository
from vd_test_log.validation import ValidationError


CREATED_AT = "2026-10-05T10:15:30+00:00"


def make_day(identifier="day-1", **changes):
    return dataclasses.replace(
        TestDay(
            id=identifier,
            date="2026-10-05",
            location="Bruntingthorpe",
            weather=None,
            notes=None,
            created_at=CREATED_AT,
        ),
        **changes,
    )


def make_setup(identifier="setup-1", day_id="day-1", order=1, **changes):
    return dataclasses.replace(
        Setup(
            id=identifier,
            test_day_id=day_id,
            name="Baseline",
            setup_code="A01",
            settings_text="Cold tyre pressures: 20 psi",
            notes=None,
            order=order,
            created_at=CREATED_AT,
        ),
        **changes,
    )


def make_event(identifier="event-1", **changes):
    return dataclasses.replace(
        EventLayout(
            id=identifier,
            track_name="Bruntingthorpe",
            layout_name="Short Loop",
            event_name="Handling Day",
            event_type="test",
            length_m=1234.5,
            notes="Dry surface",
            archived=False,
            created_at=CREATED_AT,
        ),
        **changes,
    )


def make_lap(
    identifier="lap-1",
    setup_id="setup-1",
    event_id="event-1",
    sequence=1,
    time_ms=42_318,
    status="valid",
    **changes,
):
    return dataclasses.replace(
        Lap(
            id=identifier,
            setup_id=setup_id,
            event_layout_id=event_id,
            sequence=sequence,
            time_ms=time_ms,
            status=status,
            driver="Driver 1",
            time_of_day="10:15",
            notes=None,
            created_at=CREATED_AT,
        ),
        **changes,
    )


class RepositoryTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.database_path = Path(self.temporary_directory.name) / "test_log.sqlite3"
        self.repository = SQLiteRepository.open(self.database_path)

    def tearDown(self):
        self.repository.close()
        self.temporary_directory.cleanup()

    def add_context(self):
        day = make_day()
        setup = make_setup()
        event = make_event()
        self.repository.save_day(day)
        self.repository.save_setup(setup)
        self.repository.save_event_layout(event)
        return day, setup, event

    def test_round_trips_records_and_attachments_after_reopen(self):
        day = make_day()
        setup = make_setup()
        event = make_event()
        lap = make_lap()
        day_attachment = StagedAttachment(
            role="file",
            original_name="run notes.txt",
            relative_path="days/day-1/notes-1.txt",
        )
        self.repository.save_day(day, attachments=(day_attachment,))
        self.repository.save_setup(setup)
        self.repository.save_event_layout(event)
        self.repository.save_lap(lap)

        self.repository.close()
        self.repository = SQLiteRepository.open(self.database_path)

        self.assertEqual(self.repository.get_day(day.id), day)
        self.assertEqual(self.repository.get_setup(setup.id), setup)
        self.assertEqual(self.repository.get_event_layout(event.id), event)
        self.assertEqual(self.repository.get_lap(lap.id), lap)
        attachments = self.repository.list_attachments("day", day.id)
        self.assertEqual(len(attachments), 1)
        self.assertEqual(attachments[0].owner_type, "day")
        self.assertEqual(attachments[0].owner_id, day.id)
        self.assertEqual(attachments[0].original_name, "run notes.txt")
        self.assertEqual(attachments[0].relative_path, "days/day-1/notes-1.txt")
        parsed_created_at = datetime.fromisoformat(attachments[0].created_at)
        self.assertIsNotNone(parsed_created_at.tzinfo)

    def test_new_database_uses_schema_version_one(self):
        with closing(sqlite3.connect(self.database_path)) as connection:
            version = connection.execute("PRAGMA user_version").fetchone()[0]
        self.assertEqual(version, 1)

        self.repository.close()
        self.repository = SQLiteRepository.open(self.database_path)
        self.assertEqual(self.repository.list_days(), [])

    def test_rejects_impossible_dates_and_blank_required_text(self):
        invalid_records = (
            (make_day(date="2026-02-30"), "date"),
            (make_day(location="  \t "), "location"),
            (make_setup(name=" \n "), "name"),
            (make_event(track_name=" "), "track_name"),
            (make_event(layout_name="\t"), "layout_name"),
        )
        for record, field in invalid_records:
            with self.subTest(field=field):
                with self.assertRaises(ValidationError) as raised:
                    if isinstance(record, TestDay):
                        self.repository.save_day(record)
                    elif isinstance(record, Setup):
                        self.repository.save_day(make_day())
                        self.repository.save_setup(record)
                    else:
                        self.repository.save_event_layout(record)
                self.assertEqual(raised.exception.field, field)

        self.assertEqual(self.repository.list_days(), [make_day()])

    def test_rejects_nonpositive_or_nonfinite_event_lengths(self):
        for length in (0, -1.0, math.nan, math.inf, -math.inf):
            with self.subTest(length=length):
                with self.assertRaises(ValidationError) as raised:
                    self.repository.save_event_layout(make_event(length_m=length))
                self.assertEqual(raised.exception.field, "length_m")
        self.assertEqual(self.repository.list_event_layouts(), [])

    def test_rejects_map_attachments_on_non_event_owners(self):
        day = make_day()
        invalid_map = StagedAttachment(
            role="map",
            original_name="course.png",
            relative_path="days/day-1/course.png",
        )
        with self.assertRaises(ValidationError) as raised:
            self.repository.save_day(day, attachments=(invalid_map,))
        self.assertEqual(raised.exception.field, "role")
        self.assertIsNone(self.repository.get_day(day.id))
        self.assertEqual(self.repository.list_attachments("day", day.id), [])

    def test_delete_day_returns_cascaded_attachment_rows(self):
        day, setup, event = self.add_context()
        self.repository.save_lap(make_lap())
        self.repository.save_day(
            day,
            attachments=(
                StagedAttachment("file", "day.txt", "day/day.txt"),
            ),
        )
        self.repository.save_setup(
            setup,
            attachments=(
                StagedAttachment("file", "setup.txt", "setup/setup.txt"),
            ),
        )
        self.repository.save_lap(
            make_lap(),
            attachments=(
                StagedAttachment("file", "lap.txt", "lap/lap.txt"),
            ),
        )

        removed = self.repository.delete_day(day.id)

        self.assertEqual(
            {attachment.owner_type for attachment in removed},
            {"day", "setup", "lap"},
        )
        self.assertIsNone(self.repository.get_day(day.id))
        self.assertIsNone(self.repository.get_setup(setup.id))
        self.assertIsNone(self.repository.get_lap("lap-1"))
        self.assertIsNotNone(self.repository.get_event_layout(event.id))

    def test_setup_order_edits_preserve_unique_positions_and_created_at(self):
        day = make_day()
        self.repository.save_day(day)
        first = make_setup("setup-1", order=1)
        second = make_setup("setup-2", order=2, name="Alternate")
        self.repository.save_setup(first)
        self.repository.save_setup(second)

        self.repository.save_setup(dataclasses.replace(second, order=1))
        ordered = self.repository.list_setups(day.id)
        self.assertEqual([(item.id, item.order) for item in ordered], [
            ("setup-2", 1),
            ("setup-1", 2),
        ])

        edited = dataclasses.replace(
            self.repository.get_setup(first.id),
            name="Revised",
            created_at="2000-01-01T00:00:00+00:00",
        )
        self.repository.save_setup(edited)
        self.assertEqual(self.repository.get_setup(first.id).created_at, CREATED_AT)
        self.assertEqual(
            [(item.id, item.order) for item in self.repository.list_setups(day.id)],
            [("setup-2", 1), ("setup-1", 2)],
        )

    def test_lap_reorder_requires_exact_permutation_and_renumbers_atomically(self):
        self.add_context()
        first = make_lap("lap-1", sequence=1, time_ms=45_000)
        second = make_lap("lap-2", sequence=2, time_ms=44_000)
        third = make_lap("lap-3", sequence=3, time_ms=43_000)
        for lap in (first, second, third):
            self.repository.save_lap(lap)

        with self.assertRaises(ValidationError):
            self.repository.reorder_laps("setup-1", ("lap-1", "lap-1", "lap-3"))
        with self.assertRaises(ValidationError):
            self.repository.reorder_laps("setup-1", ("lap-1", "lap-3"))

        self.repository.reorder_laps("setup-1", ("lap-3", "lap-1", "lap-2"))

        laps = self.repository.list_laps("setup-1")
        self.assertEqual([lap.id for lap in laps], ["lap-3", "lap-1", "lap-2"])
        self.assertEqual([lap.sequence for lap in laps], [1, 2, 3])

    def test_best_laps_are_valid_and_grouped_by_event(self):
        self.add_context()
        second_event = make_event(
            "event-2",
            track_name="Second Venue",
            layout_name="Long Loop",
        )
        self.repository.save_event_layout(second_event)
        for lap in (
            make_lap("lap-1", sequence=1, time_ms=50_000),
            make_lap("lap-2", sequence=2, time_ms=45_000),
            make_lap("lap-3", sequence=3, time_ms=40_000, status="invalid"),
            make_lap("lap-4", event_id="event-2", sequence=4, time_ms=60_000),
        ):
            self.repository.save_lap(lap)

        best = self.repository.best_laps_by_event("setup-1")

        self.assertEqual({event_id: lap.id for event_id, lap in best.items()}, {
            "event-1": "lap-2",
            "event-2": "lap-4",
        })

    def test_referenced_event_is_immutable_but_archiving_is_allowed(self):
        self.add_context()
        lap = make_lap()
        self.repository.save_lap(lap)

        with self.assertRaises(ValidationError) as raised:
            self.repository.save_event_layout(make_event(track_name="Changed Venue"))
        self.assertEqual(raised.exception.field, "event_layout")

        archived = self.repository.archive_event_layout("event-1")
        self.assertTrue(archived.archived)
        self.assertEqual(self.repository.get_event_layout("event-1"), archived)

        edited_historical_lap = dataclasses.replace(
            self.repository.get_lap(lap.id),
            status="invalid",
            notes="Flagged after the run",
        )
        self.repository.save_lap(edited_historical_lap)
        self.assertEqual(self.repository.get_lap(lap.id).status, "invalid")

        with self.assertRaises(ValidationError):
            self.repository.save_lap(make_lap("lap-2", sequence=2))

    def test_archived_event_cannot_be_assigned_to_new_or_reassigned_lap(self):
        self.add_context()
        second_event = make_event("event-2", layout_name="Long Loop")
        self.repository.save_event_layout(second_event)
        existing = make_lap("lap-2", event_id="event-2", sequence=1)
        self.repository.save_lap(existing)
        self.repository.archive_event_layout("event-1")

        with self.assertRaises(ValidationError):
            self.repository.save_lap(make_lap("lap-3", sequence=2))
        with self.assertRaises(ValidationError):
            self.repository.save_lap(
                dataclasses.replace(existing, event_layout_id="event-1")
            )

    def test_referenced_event_rejects_new_attachments(self):
        _, _, event = self.add_context()
        self.repository.save_lap(make_lap())
        for role in ("file", "map"):
            with self.subTest(role=role):
                attachment = StagedAttachment(
                    role=role,
                    original_name="course image",
                    relative_path=f"events/event-1/{role}.bin",
                )
                with self.assertRaises(ValidationError) as raised:
                    self.repository.save_event_layout(
                        self.repository.get_event_layout(event.id),
                        attachments=(attachment,),
                    )
                self.assertEqual(raised.exception.field, "attachments")
        self.assertEqual(self.repository.list_attachments("event_layout", event.id), [])

    def test_save_lap_holds_write_lock_across_archived_event_check(self):
        self.add_context()
        writer_results = []

        def archive_event_from_second_connection():
            connection = sqlite3.connect(self.database_path, timeout=0)
            try:
                connection.execute(
                    "UPDATE event_layouts SET archived = 1 WHERE id = ?",
                    ("event-1",),
                )
                connection.commit()
                writer_results.append("committed")
            except sqlite3.OperationalError as error:
                connection.rollback()
                writer_results.append(error.sqlite_errorname)
            finally:
                connection.close()

        original_get_event_layout = self.repository.get_event_layout
        archive_attempted = False

        def read_event_then_attempt_archive(identifier):
            nonlocal archive_attempted
            event = original_get_event_layout(identifier)
            if identifier == "event-1" and not archive_attempted:
                archive_attempted = True
                writer = Thread(target=archive_event_from_second_connection)
                writer.start()
                writer.join(timeout=3)
                self.assertFalse(writer.is_alive(), "the competing writer must finish promptly")
            return event

        self.repository.get_event_layout = read_event_then_attempt_archive
        try:
            saved = self.repository.save_lap(make_lap())
        finally:
            self.repository.get_event_layout = original_get_event_layout

        self.assertEqual(writer_results, ["SQLITE_BUSY"])
        self.assertFalse(self.repository.get_event_layout("event-1").archived)
        self.assertEqual(saved.id, "lap-1")
        self.assertEqual(self.repository.get_lap("lap-1").event_layout_id, "event-1")

    def test_reorder_laps_holds_write_lock_across_exact_membership_check(self):
        self.add_context()
        self.repository.save_lap(make_lap("lap-1", sequence=1))
        self.repository.save_lap(make_lap("lap-2", sequence=2))
        writer_results = []

        def insert_lap_from_second_connection():
            connection = sqlite3.connect(self.database_path, timeout=0)
            try:
                connection.execute(
                    """
                    INSERT INTO laps(
                        id, setup_id, event_layout_id, sequence, time_ms, status,
                        driver, time_of_day, notes, created_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        "concurrent-lap",
                        "setup-1",
                        "event-1",
                        3,
                        41_000,
                        "valid",
                        "Driver 2",
                        "10:16",
                        None,
                        CREATED_AT,
                    ),
                )
                connection.commit()
                writer_results.append("committed")
            except sqlite3.OperationalError as error:
                connection.rollback()
                writer_results.append(error.sqlite_errorname)
            finally:
                connection.close()

        membership_read_interleaved = False

        def attempt_insert_before_membership_read(statement):
            nonlocal membership_read_interleaved
            if (
                not membership_read_interleaved
                and statement.startswith("SELECT id, sequence FROM laps WHERE setup_id =")
            ):
                membership_read_interleaved = True
                writer = Thread(target=insert_lap_from_second_connection)
                writer.start()
                writer.join(timeout=3)
                self.assertFalse(writer.is_alive(), "the competing writer must finish promptly")

        self.repository._connection.set_trace_callback(attempt_insert_before_membership_read)
        try:
            reordered = self.repository.reorder_laps("setup-1", ("lap-2", "lap-1"))
        finally:
            self.repository._connection.set_trace_callback(None)

        self.assertEqual(writer_results, ["SQLITE_BUSY"])
        self.assertEqual(
            [(lap.id, lap.sequence) for lap in reordered],
            [("lap-2", 1), ("lap-1", 2)],
        )
        self.assertIsNone(self.repository.get_lap("concurrent-lap"))

    def test_delete_lap_serializes_attachment_capture_with_concurrent_writer(self):
        self.add_context()
        lap = make_lap()
        initial_attachment = StagedAttachment(
            "file",
            "before.txt",
            "laps/before.txt",
        )
        self.repository.save_lap(lap, attachments=(initial_attachment,))
        writer_results = []

        def concurrent_writer():
            connection = sqlite3.connect(self.database_path, timeout=0)
            try:
                connection.execute("PRAGMA foreign_keys = ON")
                connection.execute(
                    """
                    INSERT INTO attachments(
                        id, owner_type, lap_owner_id, role, original_name,
                        relative_path, created_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        "interleaved-attachment",
                        "lap",
                        lap.id,
                        "file",
                        "during.txt",
                        "laps/during.txt",
                        CREATED_AT,
                    ),
                )
                connection.commit()
                writer_results.append("committed")
            except sqlite3.OperationalError as error:
                connection.rollback()
                writer_results.append(error.sqlite_errorname)
            finally:
                connection.close()

        original_capture = self.repository._attachments_for_query

        def capture_then_attempt_write(predicate, parameters):
            captured = original_capture(predicate, parameters)
            writer = Thread(target=concurrent_writer)
            writer.start()
            writer.join(timeout=3)
            self.assertFalse(writer.is_alive(), "the competing writer must finish promptly")
            return captured

        self.repository._attachments_for_query = capture_then_attempt_write
        removed = self.repository.delete_lap(lap.id)

        self.assertEqual([item.relative_path for item in removed], ["laps/before.txt"])
        self.assertEqual(writer_results, ["SQLITE_BUSY"])
        self.assertIsNone(self.repository.get_lap(lap.id))
        self.assertEqual(self.repository.list_attachments("lap", lap.id), [])

    def test_delete_event_returns_attachments_and_refuses_referenced_event(self):
        self.add_context()
        map_attachment = StagedAttachment(
            "map",
            "short loop.png",
            "events/event-1/map.png",
        )
        self.repository.save_event_layout(make_event(), attachments=(map_attachment,))
        self.repository.save_lap(make_lap())

        with self.assertRaises(ValidationError):
            self.repository.delete_event_layout("event-1")
        self.assertIsNotNone(self.repository.get_event_layout("event-1"))

        self.repository.delete_lap("lap-1")
        removed = self.repository.delete_event_layout("event-1")
        self.assertEqual(len(removed), 1)
        self.assertEqual(removed[0].role, "map")
        self.assertIsNone(self.repository.get_event_layout("event-1"))

    def test_referenced_history_edits_keep_creation_timestamps(self):
        day, setup, event = self.add_context()
        lap = make_lap()
        self.repository.save_lap(lap)
        replacement_time = "2000-01-01T00:00:00+00:00"

        self.repository.save_day(dataclasses.replace(day, notes="Edited", created_at=replacement_time))
        self.repository.save_setup(dataclasses.replace(setup, notes="Edited", created_at=replacement_time))
        self.repository.save_event_layout(
            dataclasses.replace(event, archived=True, created_at=replacement_time)
        )
        self.repository.save_lap(
            dataclasses.replace(lap, notes="Edited", created_at=replacement_time)
        )

        self.assertEqual(self.repository.get_day(day.id).created_at, CREATED_AT)
        self.assertEqual(self.repository.get_setup(setup.id).created_at, CREATED_AT)
        self.assertEqual(self.repository.get_event_layout(event.id).created_at, CREATED_AT)
        self.assertEqual(self.repository.get_lap(lap.id).created_at, CREATED_AT)

    def test_refuses_nonempty_unknown_or_newer_schema_without_changing_bytes(self):
        for version in (0, 999):
            with self.subTest(version=version):
                path = Path(self.temporary_directory.name) / f"unknown-{version}.sqlite3"
                with closing(sqlite3.connect(path)) as connection:
                    connection.execute("CREATE TABLE sentinel(value TEXT)")
                    connection.execute(
                        "INSERT INTO sentinel(value) VALUES (?)",
                        ("preserve me",),
                    )
                    connection.execute(f"PRAGMA user_version = {version}")
                    connection.commit()
                before = path.read_bytes()

                with self.assertRaises(ValidationError):
                    SQLiteRepository.open(path)

                self.assertEqual(path.read_bytes(), before)
                with closing(sqlite3.connect(path)) as connection:
                    self.assertEqual(
                        connection.execute("SELECT value FROM sentinel").fetchone()[0],
                        "preserve me",
                    )


if __name__ == "__main__":
    unittest.main()
