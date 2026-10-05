import io
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from contextlib import redirect_stderr, redirect_stdout
from unittest.mock import patch
import tkinter as tk

from vd_test_log.app import run_app
from vd_test_log.data_folder import DataRootLock, resolve_data_paths


class StartupTests(unittest.TestCase):
    def test_run_app_creates_database_and_releases_lock_after_tk_exits(self):
        with tempfile.TemporaryDirectory() as folder:
            data_root = Path(folder) / "data root"
            paths = resolve_data_paths(data_root, {}, Path(folder))

            with patch.object(tk.Tk, "mainloop", return_value=None):
                self.assertEqual(run_app(data_root), 0)

            self.assertTrue(paths.database.is_file())
            second_lock = DataRootLock(paths.lock_file)
            second_lock.acquire()
            second_lock.close()

    def test_ui_build_failure_releases_data_root_lock(self):
        with tempfile.TemporaryDirectory() as folder:
            data_root = Path(folder) / "broken data"
            paths = resolve_data_paths(data_root, {}, Path(folder))
            error_output = io.StringIO()

            with patch(
                "vd_test_log.app.TestLogWindow.build",
                side_effect=RuntimeError("simulated widget build failure"),
            ):
                with redirect_stderr(error_output):
                    self.assertNotEqual(run_app(data_root), 0)

            self.assertIn("simulated widget build failure", error_output.getvalue())
            self.assertTrue(paths.database.is_file())
            self.assertTrue(paths.lock_file.is_file())
            second_lock = DataRootLock(paths.lock_file)
            second_lock.acquire()
            second_lock.close()

    def test_ui_build_failure_leaves_no_theme_error_for_next_tk_root(self):
        script = textwrap.dedent(
            """
            import tempfile
            import tkinter as tk
            from pathlib import Path
            from unittest.mock import patch
            from vd_test_log.app import run_app

            with tempfile.TemporaryDirectory() as folder:
                with patch(
                    "vd_test_log.app.TestLogWindow.build",
                    side_effect=RuntimeError("simulated widget build failure"),
                ):
                    if run_app(Path(folder) / "failed startup data") != 1:
                        raise AssertionError("the simulated startup failure must return 1")

            root = tk.Tk()
            root.withdraw()
            root.update()
            root.destroy()
            """
        )
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=Path(__file__).resolve().parents[1],
            capture_output=True,
            check=False,
            text=True,
            timeout=30,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn('can\'t invoke "event" command', result.stderr)

    def test_help_does_not_create_requested_data_root(self):
        from vd_test_log.__main__ import main

        with tempfile.TemporaryDirectory() as folder:
            data_root = Path(folder) / "help must stay offline"

            with redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    main(["--help", "--data-dir", str(data_root)])

            self.assertEqual(raised.exception.code, 0)
            self.assertFalse(data_root.exists())


if __name__ == "__main__":
    unittest.main()
