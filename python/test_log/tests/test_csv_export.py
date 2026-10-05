import csv
import os
import subprocess
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


    def _new_export_fixture(self, name):
        root = self.root / name
        paths = DataPaths(
            root=root,
            database=root / "test_log.sqlite3",
            attachments=root / "attachments",
            lock_file=root / ".test_log.lock",
            backup_root=root / "backups",
        )
        services = TestLogServices.open(paths)
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
        services.save_day(day)
        services.save_setup(setup)
        return services, paths, setup

    @staticmethod
    def _create_junction(link, target):
        quote = lambda value: str(value).replace("'", "''")
        command = (
            f"New-Item -ItemType Junction -Path '{quote(link)}' "
            f"-Target '{quote(target)}' | Out-Null"
        )
        result = subprocess.run(
            ["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", command],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode:
            raise OSError(result.stderr or result.stdout or "Could not create junction.")

    def _create_setup(self):
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
        self.services.save_day(day)
        self.services.save_setup(setup)
        return setup

    def test_rejects_resolved_database_path_before_export(self):
        setup = self._create_setup()
        database_bytes = self.paths.database.read_bytes()
        resolved_alias = self.paths.root / "not-created" / ".." / self.paths.database.name

        with self.assertRaises(ValueError):
            self.services.export_setup(setup.id, resolved_alias)

        self.assertEqual(self.paths.database.read_bytes(), database_bytes)
        self.assertFalse((self.paths.root / "not-created").exists())

    def test_rejects_hard_link_to_database_before_export(self):
        setup = self._create_setup()
        database_bytes = self.paths.database.read_bytes()
        hard_link = self.paths.root / "database-hard-link.csv"
        os.link(self.paths.database, hard_link)

        with self.assertRaises(ValueError):
            self.services.export_setup(setup.id, hard_link)

        self.assertEqual(self.paths.database.read_bytes(), database_bytes)
        self.assertEqual(hard_link.read_bytes(), database_bytes)

    def test_rejects_owned_sidecars_attachments_and_backup_files(self):
        target_factories = (
            lambda paths: Path(f"{paths.database}-wal"),
            lambda paths: Path(f"{paths.database}-shm"),
            lambda paths: Path(f"{paths.database}-journal"),
            lambda paths: paths.attachments / "owned-attachment.bin",
            lambda paths: paths.backup_root / "existing-backup" / "attachments" / "owned.bin",
        )
        for index, target_factory in enumerate(target_factories):
            with self.subTest(target_index=index):
                services, paths, setup = self._new_export_fixture(f"protected-{index}")
                try:
                    destination = target_factory(paths)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(f"preserve-{index}".encode("ascii"))
                    before = destination.read_bytes()

                    error = None
                    try:
                        services.export_setup(setup.id, destination)
                    except Exception as caught:
                        error = caught

                    self.assertEqual(destination.read_bytes(), before)
                    self.assertIsInstance(error, ValueError)
                finally:
                    services.close()

    def test_rejects_resolved_attachment_junction_alias(self):
        setup = self._create_setup()
        attachment = self.paths.attachments / "owned-through-junction.bin"
        attachment.write_bytes(b"preserve attachment")
        junction = self.root / "attachments alias"
        self._create_junction(junction, self.paths.attachments)
        destination = junction / attachment.name

        with self.assertRaises(ValueError):
            self.services.export_setup(setup.id, destination)

        self.assertEqual(attachment.read_bytes(), b"preserve attachment")

    def test_rejects_data_root_lock_before_opening_destination(self):
        setup = self._create_setup()
        lock_file = self.paths.lock_file
        original_size = lock_file.stat().st_size

        try:
            self.services.export_setup(setup.id, lock_file)
        except Exception as error:
            self.assertIsInstance(error, ValueError)
        else:
            self.fail("export must reject the active data-root lock")

        self.assertEqual(lock_file.stat().st_size, original_size)
        self.services.close()
        self.assertEqual(lock_file.read_bytes(), b"\x00")

    def test_allows_ordinary_csv_inside_data_root(self):
        setup = self._create_setup()
        destination = self.paths.root / "exports" / "laps.csv"

        result = self.services.export_setup(setup.id, destination)

        self.assertEqual(result, destination)
        self.assertTrue(destination.is_file())
        self.assertTrue(destination.read_text(encoding="utf-8").startswith("test_date,test_location"))


if __name__ == "__main__":
    unittest.main()
