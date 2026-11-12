"""Command-line entry point for IDEA-GP, formerly HAGSSuite."""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


LEGACY_DIR = Path(__file__).resolve().parent / "legacy"
HAGS_RUNNER = LEGACY_DIR / "optimize_gs.py"


def _run_hags(argv: list[str]) -> None:
    if not HAGS_RUNNER.exists():
        raise SystemExit(f"HAGSSuite runner not found: {HAGS_RUNNER}")
    old_argv = sys.argv[:]
    try:
        sys.argv = [str(HAGS_RUNNER), *argv]
        runpy.run_path(str(HAGS_RUNNER), run_name="__main__")
    finally:
        sys.argv = old_argv


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="idea-gp",
        description="IDEA-GP genomic prediction workflow, formerly HAGSSuite.",
    )
    parser.add_argument("command", choices=["run", "status"])
    parser.add_argument("args", nargs=argparse.REMAINDER)

    parsed = parser.parse_args(argv)
    if parsed.command == "status":
        if parsed.args:
            parser.error("status does not accept additional arguments")
        print("idea-gp is the renamed HAGSSuite module.")
        print(f"Runner: {HAGS_RUNNER}")
        print(f"Runner attached: {HAGS_RUNNER.exists()}")
        return

    if parsed.command == "run":
        # A positional command consumes the wrapper's arguments before the
        # runner flags. Accept an optional separator without passing it to the
        # runner, where it would prematurely terminate option parsing.
        runner_args = parsed.args[1:] if parsed.args[:1] == ["--"] else parsed.args
        _run_hags(runner_args)


if __name__ == "__main__":
    main()
