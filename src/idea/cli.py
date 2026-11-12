"""Command-line wrappers for the current IDEA legacy scripts."""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


LEGACY_DIR = Path(__file__).resolve().parent / "legacy"


COMMANDS = {
    "generate-features": "generate_feature.py",
    "select": "select_snp_no_genome.py",
    "run": "uldsuite_acc.py",
    "run-legacy": "uldsuite.py",
    "accuracy": "calculate_true_accuracy.py",
}


def _run_legacy(script_name: str, argv: list[str]) -> None:
    script_path = LEGACY_DIR / script_name
    if not script_path.exists():
        raise SystemExit(f"Legacy script not found: {script_path}")

    old_argv = sys.argv[:]
    try:
        sys.argv = [str(script_path), *argv]
        runpy.run_path(str(script_path), run_name="__main__")
    finally:
        sys.argv = old_argv


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="idea-design",
        description="IDEA low-density SNP array design and imputation workflow.",
    )
    parser.add_argument(
        "command",
        choices=sorted(COMMANDS),
        help="Workflow command to run.",
    )
    parser.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="Arguments passed to the selected workflow command.",
    )
    parsed = parser.parse_args(argv)
    _run_legacy(COMMANDS[parsed.command], parsed.args)


if __name__ == "__main__":
    main()
