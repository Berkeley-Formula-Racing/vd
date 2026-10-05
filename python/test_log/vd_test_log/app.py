"""Application startup for the offline vehicle test log."""

from __future__ import annotations

import argparse
import os
import sys
import tkinter as tk
from collections.abc import Sequence
from pathlib import Path

from .data_folder import AlreadyRunningError, resolve_data_paths
from .services import TestLogServices
from .ui import TestLogWindow


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vd_test_log",
        description="Record offline vehicle test days, setups, laps, and files.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        metavar="PATH",
        help="use PATH for the local database, attachments, and backups",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    options = _argument_parser().parse_args(argv)
    return run_app(options.data_dir)


def run_app(data_dir: Path | None = None) -> int:
    """Open the data root, run Tk, and release storage on every exit path."""
    services: TestLogServices | None = None
    root: tk.Tk | None = None
    exit_code = 0

    try:
        paths = resolve_data_paths(data_dir, os.environ, Path.home())
        services = TestLogServices.open(paths)
        root = tk.Tk()
        root.withdraw()
        root.update_idletasks()
        window = TestLogWindow(root, services)
        window.build()
        root.deiconify()
        root.mainloop()
    except Exception as error:
        if isinstance(error, AlreadyRunningError):
            message = str(error)
        else:
            message = f"Could not start the test log: {error}"
        print(message, file=sys.stderr)
        exit_code = 1
    finally:
        if root is not None:
            try:
                root.destroy()
            except tk.TclError:
                pass
        if services is not None:
            try:
                services.close()
            except Exception as error:
                print(f"Could not close the test log data safely: {error}", file=sys.stderr)
                exit_code = 1

    return exit_code
