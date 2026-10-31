"""Command line of the orbital-projected band plotter; stage 1 offers only "inspect"."""

from __future__ import annotations

import argparse
import sys

from .parser import VasprunError, parse_vasprun


# Build the argument parser with the stage 1 subcommand.
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="vasp-orbital-bands",
                                     description="Orbital-projected band structures from vasprun.xml (stage 1: inspect only).")
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser("inspect", help="print what the parser reads from a vasprun.xml")
    inspect.add_argument("vasprun", help="path to vasprun.xml or vasprun.xml.gz")
    inspect.add_argument("--allow-missing-projections", action="store_true",
                         help="print the summary even when the file has no <projected> block")
    return parser


# Entry point of the vasp-orbital-bands command; returns the process exit code.
def main(argv: list[str] | None = None) -> int:
    # "inspect" is the only subcommand, and add_subparsers(required=True) exits on anything else.
    args = build_parser().parse_args(argv)
    try:
        data = parse_vasprun(args.vasprun, require_projections=not args.allow_missing_projections)
    except (VasprunError, FileNotFoundError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print("\n".join(data.summary_lines()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
