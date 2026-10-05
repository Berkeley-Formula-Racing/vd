import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from vd_test_log.data_folder import AlreadyRunningError, DataRootLock, resolve_data_paths


class DataFolderTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_explicit_data_directory_is_resolved_without_using_profile(self):
        selected = self.root / "custom data"
        paths = resolve_data_paths(selected, {}, self.root / "home")

        self.assertEqual(paths.root, selected.resolve())
        self.assertEqual(paths.database, selected.resolve() / "test_log.sqlite3")
        self.assertEqual(paths.attachments, selected.resolve() / "attachments")
        self.assertEqual(paths.lock_file, selected.resolve() / ".test_log.lock")
        self.assertEqual(paths.backup_root, selected.resolve() / "backups")

    def test_local_app_data_is_used_when_available(self):
        local = self.root / "local"

        paths = resolve_data_paths(None, {"LOCALAPPDATA": str(local)}, self.root / "home")

        self.assertEqual(paths.root, (local / "VDTestLog").resolve())

    def test_home_fallback_is_used_when_local_app_data_is_empty(self):
        home = self.root / "home"

        paths = resolve_data_paths(None, {"LOCALAPPDATA": ""}, home)

        self.assertEqual(paths.root, (home / "VDTestLog").resolve())

    def test_lock_blocks_another_process_and_releases_on_close(self):
        lock = DataRootLock(self.root / ".test_log.lock")
        lock.acquire()
        code = """
from pathlib import Path
import sys
from vd_test_log.data_folder import AlreadyRunningError, DataRootLock
lock = DataRootLock(Path(sys.argv[1]))
try:
    lock.acquire()
except AlreadyRunningError:
    raise SystemExit(23)
lock.close()
"""
        try:
            child = subprocess.run(
                [sys.executable, "-c", code, str(lock.path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(child.returncode, 23, child.stderr)
        finally:
            lock.close()

        second = DataRootLock(lock.path)
        second.acquire()
        second.close()

    def test_lock_close_is_idempotent(self):
        lock = DataRootLock(self.root / ".test_log.lock")
        lock.acquire()

        lock.close()
        lock.close()


if __name__ == "__main__":
    unittest.main()
