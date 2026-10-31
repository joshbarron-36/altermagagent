"""Command line of the orbital-projected band plotter: "inspect", "plot", and plot flags directly on the main command.

Examples (the first one plots from vasprun.xml in the current folder):
  vasp-orbital-bands --element Mn --orbitals d --output mn_d.png
  vasp-orbital-bands plot --vasprun run/vasprun.xml --atoms 1,3 --orbitals dxy,dxz --mode color --output bands.png
  vasp-orbital-bands --list-atoms
  vasp-orbital-bands inspect --vasprun run/vasprun.xml --list-orbitals
Exit codes: 0 for success, 2 for a user error, printed as one "error:" line.
Plotting, k-path and Matplotlib modules are imported only when a figure is drawn, so listings and inspect start quickly.
"""

from __future__ import annotations

import argparse
import dataclasses
import re
import sys
from pathlib import Path
from typing import NoReturn, Sequence

from altermagagent.bands import slug

from .options import DEFAULT_WINDOW, FORMATS, MODES, KPathError, OptionsError, PlotOptions
from .parser import VasprunError, parse_vasprun
from .projections import SelectionError, atom_lines, build_selection, orbital_lines, select_component, select_weights

DEFAULTS = PlotOptions()


class UsageError(Exception):
    """A command-line argument is missing or malformed."""


class _Parser(argparse.ArgumentParser):
    """ArgumentParser that raises UsageError instead of exiting, and accepts each option only under its full name."""

    # Turn off prefix matching, so "--ato" or "--elem" is refused rather than completed to a longer option.
    def __init__(self, **kwargs: object) -> None:
        super().__init__(allow_abbrev=False, **kwargs)

    # Raise the message so main() can print it as a single "error:" line.
    def error(self, message: str) -> NoReturn:
        raise UsageError(message)


# Add --vasprun, --list-atoms and --list-orbitals, shared by "inspect" and "plot".
def _add_file_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--vasprun", default="vasprun.xml", help="path to vasprun.xml or .xml.gz (default: vasprun.xml)")
    parser.add_argument("--list-atoms", action="store_true", help="print every ion with its 1-based VASP number and stop")
    parser.add_argument("--list-orbitals", action="store_true", help="print the file's orbital labels and groups and stop")


# Parser of the "inspect" subcommand.
def inspect_parser() -> _Parser:
    parser = _Parser(prog="vasp-orbital-bands inspect", description="Print what the parser reads from a vasprun.xml.")
    _add_file_flags(parser)
    parser.add_argument("--allow-missing-projections", action="store_true",
                        help="print the summary even when the file has no <projected> block")
    return parser


# Parse "6,4.5" or "6x4.5" into a (width, height) pair in inches.
def _figsize(text: str) -> tuple[float, float]:
    parts = re.split(r"[,x]", text.strip())
    try:
        width, height = (float(p) for p in parts)
    except ValueError:
        raise argparse.ArgumentTypeError(f"figure size {text!r} must be WIDTH,HEIGHT in inches, for example 6,4.5") from None
    return width, height


# Parser of the "plot" subcommand, also used when the main command gets plot flags directly; defaults come from PlotOptions.
def plot_parser(prog: str = "vasp-orbital-bands") -> _Parser:
    d = DEFAULTS
    parser = _Parser(prog=prog, description="Orbital-projected band structure from vasprun.xml. "
                     "Subcommands: inspect, plot; without one, the flags below plot directly.")
    _add_file_flags(parser)
    pick = parser.add_argument_group("selection (ions use VASP's 1-based numbers)")
    pick.add_argument("--atoms", "--atom", help="ion numbers and ranges, for example 1,3-4 (--atom is accepted too)")
    pick.add_argument("--element", help="every ion of these elements, for example Mn or Mn,Te")
    pick.add_argument("--orbitals", "--orbital",
                      help="orbitals or groups (s, p, d, f, t2g, eg, all), for example dxy,eg (--orbital is accepted too)")
    pick.add_argument("--select", action="append", metavar="ATOMS:ORBITALS",
                      help="one atom and orbital combination, repeatable, for example --select Mn:d --select 3:pz")
    pick.add_argument("--normalize", action="store_true",
                      help="divide by the sum over all ions and orbitals of each state (default: raw weights)")
    look = parser.add_argument_group("plot")
    look.add_argument("--mode", choices=MODES, default=d.mode, help="fatband (marker size) or color (colour and colour bar)")
    look.add_argument("--emin", type=float, default=d.emin,
                      help=f"lower energy limit in eV, relative to the Fermi level unless --absolute (default -{DEFAULT_WINDOW:g})")
    look.add_argument("--emax", type=float, default=d.emax,
                      help=f"upper energy limit in eV, relative to the Fermi level unless --absolute (default {DEFAULT_WINDOW:g})")
    look.add_argument("--absolute", action="store_true", help="plot absolute energies instead of E - E_F")
    look.add_argument("--efermi", type=float, help="Fermi energy in eV to use instead of the file's value")
    look.add_argument("--high-symmetry", metavar="N:NAME,...",
                      help='high-symmetry labels at 1-based k-point numbers, for example "1:G,41:M,81:K,121:G"')
    look.add_argument("--output", help=f"output file ({', '.join(FORMATS)}); default: a name built from the selection")
    look.add_argument("--dpi", type=int, default=d.dpi, help=f"resolution of raster output (default {d.dpi})")
    look.add_argument("--figsize", type=_figsize, default=d.figsize, help="WIDTH,HEIGHT in inches (default %(default)s)")
    look.add_argument("--font-size", type=float, default=d.font_size, help="base font size in points (default %(default)s)")
    look.add_argument("--linewidth", type=float, default=d.linewidth, help="band line width in points (default %(default)s)")
    look.add_argument("--marker-scale", type=float, default=d.marker_scale,
                      help="fat-band marker area in points squared per unit weight (default %(default)s)")
    look.add_argument("--alpha", type=float, default=d.alpha, help="fat-band marker opacity, 0 to 1 (default %(default)s)")
    look.add_argument("--cmap", default=d.cmap, help="Matplotlib colour map of the color mode (default %(default)s)")
    look.add_argument("--cmin", type=float, default=d.cmin, help="lower colour limit (default 0)")
    look.add_argument("--cmax", type=float, default=d.cmax, help="upper colour limit (default: largest weight in the window)")
    look.add_argument("--title", default=d.title, help="figure title (default: the selection)")
    look.add_argument("--no-title", action="store_true", help="draw no title")
    look.add_argument("--legend", action="store_true", help="draw a legend naming the selection")
    return parser


# PlotOptions from the parsed flags: every PlotOptions field has a flag of the same name.
def options_from_args(args: argparse.Namespace) -> PlotOptions:
    return PlotOptions(**{f.name: getattr(args, f.name) for f in dataclasses.fields(PlotOptions)})


# Output path: --output as given, or a default name in the current folder that must not be the vasprun's folder.
def _output_path(args: argparse.Namespace, label: str, vasprun: Path) -> Path:
    if args.output:
        return Path(args.output).expanduser()
    out = Path.cwd() / f"bands_{slug(label)}_{args.mode}{'_norm' if args.normalize else ''}.png"
    if out.parent.resolve() == vasprun.resolve().parent:
        raise UsageError(f"the default output {out.name} would be written next to {vasprun.name}; "
                         f"pass --output to choose the file, for example --output {out.name}")
    return out


# Print the requested --list-atoms and --list-orbitals listings; the file is parsed without requiring a <projected> block.
def _print_lists(args: argparse.Namespace, vasprun: str | Path) -> None:
    data = parse_vasprun(vasprun, require_projections=False, with_projections=False)
    if args.list_atoms:
        print("\n".join(atom_lines(data)))
    if args.list_orbitals:
        if not data.orbital_labels:
            raise SelectionError("this file has no <projected> block, so it has no orbital labels; "
                                 "rerun the band calculation with LORBIT = 11 to get them")
        print("\n".join(orbital_lines(data)))


# Run "inspect": print the summary, or only the requested listings; the projection numbers are never read.
def run_inspect(args: argparse.Namespace) -> int:
    if args.list_atoms or args.list_orbitals:
        _print_lists(args, args.vasprun)
        return 0
    data = parse_vasprun(args.vasprun, require_projections=not args.allow_missing_projections, with_projections=False)
    print("\n".join(data.summary_lines()))
    return 0


# Run a plot: parse, select, build the k-path, draw and save; print the output path and every note.
def run_plot(args: argparse.Namespace) -> int:
    vasprun = Path(args.vasprun).expanduser()
    if args.list_atoms or args.list_orbitals:
        _print_lists(args, vasprun)
        return 0
    from .kpath import build_kpath
    from .plotting import plot_bands, resolve_window
    options = options_from_args(args)
    data = parse_vasprun(vasprun)
    if args.efermi is not None:
        data = dataclasses.replace(data, efermi=args.efermi)
    selection = build_selection(data, atoms=args.atoms, elements=args.element, orbitals=args.orbitals, select=args.select)
    component, component_label, component_note = select_component(data)
    weights = select_weights(data, selection, component, normalize=options.normalize)
    efermi, emin, emax = resolve_window(data, options.emin, options.emax, options.absolute)
    kpath = build_kpath(data, labels=args.high_symmetry)
    output = _output_path(args, selection.label, vasprun)
    title = None if args.no_title else (options.title or f"{selection.label} ({component_label}"
                                        f"{', normalised' if options.normalize else ''})")
    options = dataclasses.replace(options, emin=emin, emax=emax, title=title)
    notes = [component_note] if component_note else []
    path = plot_bands(data.eigenvalues[0], weights, kpath, output, options, efermi=efermi, label=selection.label,
                      notes=notes)
    nions, norbitals = len({ion for ion, _ in selection.pairs}), len({orb for _, orb in selection.pairs})
    print(f"wrote {path}")
    print(f"selection: {selection.label}; {nions} ions, {norbitals} orbitals, block {component_label}")
    print(f"Fermi energy {efermi:.6f} eV; k-path segments from {kpath.segment_source}; labels from {kpath.label_source}")
    for note in notes + kpath.notes:
        print(f"note: {note}")
    return 0


# Entry point of the vasp-orbital-bands command; returns the process exit code.
def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        if argv[:1] == ["inspect"]:
            return run_inspect(inspect_parser().parse_args(argv[1:]))
        if argv[:1] == ["plot"]:
            return run_plot(plot_parser("vasp-orbital-bands plot").parse_args(argv[1:]))
        return run_plot(plot_parser().parse_args(argv))
    except (UsageError, argparse.ArgumentError, VasprunError, SelectionError, KPathError, OptionsError,
            FileNotFoundError, IsADirectoryError, PermissionError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
