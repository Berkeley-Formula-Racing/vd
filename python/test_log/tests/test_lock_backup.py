import sqlite3
import tempfile
import unittest
from pathlib import Path

from vd_test_log.backup import BackupError
from vd_test_log.data_folder import DataRootLock
from vd_test_log.models import DataPaths, TestDay, new_id, utc_now_iso
from vd_test_log.repository import SQLiteRepository
from vd_test_log.services import BackupInProgressError, TestLogServices


class LockBackupTests(unittest.TestCase):
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

    def make_day(self):
        return TestDay(
            id=new_id(),
            date="2026-10-05",
            location="Bruntingthorpe",
            weather=None,
            notes=None,
            created_at=utc_now_iso(),
        )

    def create_junction(self, link, target):
        if not str(link).startswith("\\\\") and not Path(link).drive:
            self.fail("The junction test requires an absolute Windows path.")
        command = (
            f"New-Item -ItemType Junction -Path '{str(link).replace(chr(39), chr(39) * 2)}' "
            f"-Target '{str(target).replace(chr(39), chr(39) * 2)}' | Out-Null"
        )
        result = __import__("subprocess").run(
            ["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", command],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

    def test_backup_is_reopenable_and_copies_only_snapshot_manifest(self):
        source = self.root / "referenced.txt"
        source.write_bytes(b"referenced contents")
        referenced_stage = self.services.stage_attachment(source, "file")
        day = self.make_day()
        self.services.save_day(day, (referenced_stage,))
        referenced = self.services.list_attachments("day", day.id)[0]

        orphan_source = self.root / "orphan.txt"
        orphan_source.write_bytes(b"unreferenced contents")
        orphan = self.services.stage_attachment(orphan_source, "file")

        backup = self.services.create_backup()
        snapshot = SQLiteRepository.open(backup / "test_log.sqlite3")
        try:
            self.assertEqual(snapshot.list_days(), [day])
            self.assertEqual(snapshot.list_attachments("day", day.id), [referenced])
        finally:
            snapshot.close()
        self.assertEqual((backup / referenced.relative_path).read_bytes(), b"referenced contents")
        self.assertFalse((backup / orphan.relative_path).exists())

    def test_empty_backup_still_has_an_attachments_folder(self):
        backup = self.services.create_backup()

        self.assertTrue((backup / "attachments").is_dir())

    def test_failed_backup_keeps_previous_backup_and_removes_its_temp_folder(self):
        source = self.root / "referenced.txt"
        source.write_bytes(b"will go missing")
        staged = self.services.stage_attachment(source, "file")
        day = self.make_day()
        self.services.save_day(day, (staged,))
        referenced = self.services.list_attachments("day", day.id)[0]
        self.services.resolve_attachment(referenced).unlink()

        previous = self.paths.backup_root / "previous-backup"
        previous.mkdir(parents=True)
        marker = previous / "keep.txt"
        marker.write_text("valid earlier backup", encoding="utf-8")

        with self.assertRaises(BackupError):
            self.services.create_backup()

        self.assertEqual(marker.read_text(encoding="utf-8"), "valid earlier backup")
        self.assertEqual(list(self.paths.backup_root.iterdir()), [previous])

    def test_backup_rejects_backup_root_outside_data_root_without_touching_it(self):
        selected_root = self.root / "malformed-data"
        outside = self.root / "outside-backups"
        previous = outside / "prior-backup"
        previous.mkdir(parents=True)
        marker = previous / "keep.txt"
        marker.write_text("prior backup", encoding="utf-8")
        source = self.root / "original-source.txt"
        source.write_bytes(b"original bytes")
        paths = DataPaths(
            root=selected_root,
            database=selected_root / "test_log.sqlite3",
            attachments=selected_root / "attachments",
            lock_file=selected_root / ".test_log.lock",
            backup_root=outside,
        )
        services = TestLogServices.open(paths)
        try:
            day = self.make_day()
            staged = services.stage_attachment(source, "file")
            services.save_day(day, (staged,))

            with self.assertRaises(BackupError):
                services.create_backup()

            self.assertEqual(marker.read_text(encoding="utf-8"), "prior backup")
            self.assertEqual(list(outside.iterdir()), [previous])
            self.assertEqual(source.read_bytes(), b"original bytes")
            self.assertEqual(services.list_days(), [day])
        finally:
            services.close()

    def test_open_does_not_create_malformed_backup_root_outside_selected_data_root(self):
        selected_root = self.root / "uncreated-data"
        outside = self.root / "must-not-be-created"
        paths = DataPaths(
            root=selected_root,
            database=selected_root / "test_log.sqlite3",
            attachments=selected_root / "attachments",
            lock_file=selected_root / ".test_log.lock",
            backup_root=outside,
        )
        services = TestLogServices.open(paths)
        try:
            self.assertFalse(outside.exists())
            with self.assertRaises(BackupError):
                services.create_backup()
            self.assertFalse(outside.exists())
        finally:
            services.close()

    def test_backup_rejects_backup_junction_without_writing_through_it(self):
        outside = self.root / "junction-target"
        previous = outside / "prior-backup"
        previous.mkdir(parents=True)
        marker = previous / "keep.txt"
        marker.write_text("prior backup", encoding="utf-8")
        source = self.root / "original-source.txt"
        source.write_bytes(b"original bytes")
        day = self.make_day()
        staged = self.services.stage_attachment(source, "file")
        self.services.save_day(day, (staged,))
        self.create_junction(self.paths.backup_root, outside)

        try:
            with self.assertRaises(BackupError):
                self.services.create_backup()

            self.assertEqual(marker.read_text(encoding="utf-8"), "prior backup")
            self.assertEqual(list(outside.iterdir()), [previous])
            self.assertEqual(source.read_bytes(), b"original bytes")
            self.assertEqual(self.services.list_days(), [day])
        finally:
            self.paths.backup_root.rmdir()

    def test_failed_service_startup_releases_data_root_lock(self):
        self.services.close()
        self.paths.database.write_bytes(b"not a sqlite database")

        with self.assertRaises(sqlite3.DatabaseError):
            TestLogServices.open(self.paths)

        lock = DataRootLock(self.paths.lock_file)
        lock.acquire()
        lock.close()

    def test_facade_mutation_is_rejected_during_backup(self):
        manager = self.services._backup_manager
        real_create_backup = manager.create_backup
        attempted_day = self.make_day()
        observed = []

        def attempt_during_backup(connection, attachment_store):
            try:
                self.services.save_day(attempted_day)
            except BackupInProgressError:
                observed.append("blocked")
            return real_create_backup(connection, attachment_store)

        manager.create_backup = attempt_during_backup
        self.services.create_backup()

        self.assertEqual(observed, ["blocked"])
        self.assertIsNone(self.services.get_day(attempted_day.id))


if __name__ == "__main__":
    unittest.main()
