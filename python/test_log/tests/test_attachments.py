import subprocess
import tempfile
import unittest
from pathlib import Path

from vd_test_log.attachments import (
    AttachmentMissingError,
    AttachmentPathError,
    AttachmentStore,
)
from vd_test_log.models import Attachment, DataPaths


class AttachmentStoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.data_root = self.root / "data"
        self.paths = DataPaths(
            root=self.data_root,
            database=self.data_root / "test_log.sqlite3",
            attachments=self.data_root / "attachments",
            lock_file=self.data_root / ".test_log.lock",
            backup_root=self.data_root / "backups",
        )
        self.store = AttachmentStore(self.paths)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_stage_copies_to_a_generated_path_and_preserves_source(self):
        source = self.root / "setup notes.txt"
        source.write_bytes(b'front="2", rear="3"')

        staged = self.store.stage(source, "file")

        stored = self.data_root / staged.relative_path
        self.assertEqual(source.read_bytes(), b'front="2", rear="3"')
        self.assertEqual(stored.read_bytes(), source.read_bytes())
        self.assertEqual(Path(staged.relative_path).parts[0], "attachments")
        self.assertNotEqual(stored, source)
        self.assertEqual(staged.original_name, source.name)

    def test_resolve_rejects_traversal_before_returning_a_path(self):
        attachment = Attachment(
            id="attachment-1",
            owner_type="setup",
            owner_id="setup-1",
            role="file",
            original_name="outside.txt",
            relative_path="../outside.txt",
            created_at="2026-10-05T10:15:30+00:00",
        )

        with self.assertRaises(AttachmentPathError):
            self.store.resolve(attachment)

    def test_resolve_rejects_a_junction_that_escapes_attachment_root(self):
        outside = self.root / "external"
        outside.mkdir()
        outside_file = outside / "escape.txt"
        outside_file.write_text("private", encoding="utf-8")
        self.paths.attachments.mkdir(parents=True)
        link = self.paths.attachments / "escape"
        junction_command = (
            f"New-Item -ItemType Junction -Path '{str(link).replace(chr(39), chr(39) * 2)}' "
            f"-Target '{str(outside).replace(chr(39), chr(39) * 2)}' | Out-Null"
        )
        result = subprocess.run(
            ["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", junction_command],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            self.skipTest(f"junction creation is unavailable: {result.stderr or result.stdout}")
        attachment = Attachment(
            id="attachment-2",
            owner_type="setup",
            owner_id="setup-1",
            role="file",
            original_name="escape.txt",
            relative_path="attachments/escape/escape.txt",
            created_at="2026-10-05T10:15:30+00:00",
        )

        try:
            with self.assertRaises(AttachmentPathError):
                self.store.resolve(attachment)
            self.assertEqual(outside_file.read_text(encoding="utf-8"), "private")
        finally:
            link.rmdir()

    def test_missing_stored_file_remains_a_typed_open_error(self):
        attachment = Attachment(
            id="attachment-3",
            owner_type="lap",
            owner_id="lap-1",
            role="file",
            original_name="removed.txt",
            relative_path="attachments/removed.txt",
            created_at="2026-10-05T10:15:30+00:00",
        )

        with self.assertRaises(AttachmentMissingError):
            self.store.resolve(attachment)

    def test_stage_reports_a_missing_source_without_creating_a_copy(self):
        with self.assertRaises(OSError):
            self.store.stage(self.root / "missing.txt", "file")

        self.assertEqual(list(self.paths.attachments.glob("*")), [])


if __name__ == "__main__":
    unittest.main()
