import csv
import tempfile
import unittest
from pathlib import Path

from vd_test_log.models import DataPaths, EventLayout, Lap, Setup, TestDay, new_id, utc_now_iso
from vd_test_log.services import TestLogServices


class CsvExportTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.paths = DataPaths(
            root=self.root / "data",
            database=self.root / "data" / "test_log.sqlite3",
            attachments=self.root / "data" / "attachments",
            lock_file=self.root / "data" / ".test_log.lock",
            backup_root=self.root / "data" / "backups",
        )
        self.services = TestLogServices.open(self.paths)

    def tearDown(self):
        self.services.close()
        self.temporary_directory.cleanup()

    def test_export_uses_canonical_snake_case_columns_and_quotes_invalid_rows(self):
        created_at = utc_now_iso()
        day = TestDay(
            id=new_id(),
            date="2026-10-05",
            location="Bruntingthorpe",
            weather=None,
            notes=None,
            created_at=created_at,
        )
        setup = Setup(
            id=new_id(),
            test_day_id=day.id,
            name="Baseline",
            setup_code="A01",
            settings_text="pressures",
            notes=None,
            order=1,
            created_at=created_at,
        )
        event = EventLayout(
            id=new_id(),
            track_name="Bruntingthorpe",
            layout_name="Short Loop",
            event_name="Handling Day",
            event_type="test",
            length_m=1234.5,
            notes=None,
            archived=False,
            created_at=created_at,
        )
        lap = Lap(
            id=new_id(),
            setup_id=setup.id,
            event_layout_id=event.id,
            sequence=1,
            time_ms=123_456,
            status="invalid",
            driver="Driver 1",
            time_of_day="10:15",
            notes='off line, said "driver"\ncheck pressure',
            created_at=created_at,
        )
        self.services.save_day(day)
        self.services.save_setup(setup)
        self.services.save_event_layout(event)
        self.services.save_lap(lap)
        destination = self.root / "exports" / "laps.csv"

        result = self.services.export_setup(setup.id, destination)

        self.assertEqual(result, destination)
        with destination.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            self.assertEqual(
                reader.fieldnames,
                [
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
                ],
            )
            rows = list(reader)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["test_date"], "2026-10-05")
        self.assertEqual(rows[0]["setup_name"], "Baseline")
        self.assertEqual(rows[0]["setup_code"], "A01")
        self.assertEqual(rows[0]["lap_time"], "2:03.456")
        self.assertEqual(rows[0]["status"], "invalid")
        self.assertEqual(rows[0]["notes"], 'off line, said "driver"\ncheck pressure')
        self.assertEqual(rows[0]["track_name"], "Bruntingthorpe")


if __name__ == "__main__":
    unittest.main()
