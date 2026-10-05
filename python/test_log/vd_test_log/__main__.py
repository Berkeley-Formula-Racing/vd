"""Command-line entry point for the offline vehicle test log."""

import argparse
from pathlib import Path

from .app import run_app


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="vd_test_log")
    parser.add_argument("--data-dir", type=Path, metavar="PATH")
    options = parser.parse_args(argv)
    return run_app(options.data_dir)


if __name__ == "__main__":
    raise SystemExit(main())
