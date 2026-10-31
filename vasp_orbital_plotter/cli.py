"""Command line of the orbital-projected band plotter: "inspect", "plot", and plot flags directly on the main command.

Examples (the first one plots from vasprun.xml in the current folder):
  vasp-orbital-bands --element Mn --orbitals d --output mn_d.png
  vasp-orbital-bands plot --vasprun run/vasprun.xml --atoms 1,3 --orbitals dxy,dxz --mode color --output bands.png
  vasp-orbital-bands --list-atoms
  vasp-orbital-bands inspect --vasprun run/vasprun.xml --list-orbitals
  vasp-orbital-bands --vasprun cr/vasprun.xml --element Cr --orbitals d --spin both --output cr_d_spins.png
  vasp-orbital-bands --element Mn --orbitals d --component mz --mode color --output mn_d_mz.png
  vasp-orbital-bands --select Mn:d --select Te:p --overlay --style talk --output overlay.pdf
Exit codes: 0 for success, 2 for a user error, printed as one "error:" line.
Plotting, k-path and Matplotlib modules are imported only when a figure is drawn, so listings and inspect start quickly.
"""

from __future__ import annotations

import argparse
import dataclasses
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING, NoReturn, Sequence

import numpy as np

from altermagagent.bands import slug

from .options import (DEFAULT_WINDOW, FORMATS, MAX_SERIES, MODES, PANEL_SIZE, SERIES_COLOURS, STYLES, KPathError, OptionsError,
                      PlotOptions)
from .parser import VasprunData, VasprunError, parse_vasprun
from .projections import (COMPONENTS, SPINS, BlockChoice, SelectionError, Selection, atom_lines, bands_reaching,
                          build_selection, build_series, choose_blocks, flag_refusal, orbital_lines, select_weights,
                          state_totals)

if TYPE_CHECKING:
    from .plotting import Panel

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
    pick.add_argument("--overlay", action="store_true",
                      help="draw each --select term as its own series with its own colour and legend entry, instead of their sum")
    pick.add_argument("--label-series", action="append", metavar="NAME",
                      help="legend name of one overlay series, repeatable, in the order of the --select terms")
    spin = parser.add_argument_group("spin and magnetisation")
    spin.add_argument("--spin", choices=SPINS,
                      help="collinear ISPIN = 2 files only: one spin channel, or both as two panels (default both)")
    spin.add_argument("--component", choices=COMPONENTS,
                      help="noncollinear files only: total charge (default) or the signed magnetisation block mx, my or mz")
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
    look.add_argument("--label-map", metavar="OLD=NEW,...",
                      help='rename path labels, mathtext allowed, for example "G=Γ,C0=C$_0$"')
    look.add_argument("--no-labels", action="store_true", help="draw the high-symmetry lines without their names")
    look.add_argument("--output", help=f"output file ({', '.join(FORMATS)}); default: a name built from the selection")
    look.add_argument("--dpi", type=int, default=d.dpi, help=f"resolution of raster output (default {d.dpi})")
    look.add_argument("--figsize", type=_figsize, default=d.figsize,
                      help=f"WIDTH,HEIGHT in inches (default {PANEL_SIZE[0]:g},{PANEL_SIZE[1]:g} per panel)")
    look.add_argument("--font-size", type=float, default=d.font_size, help="base font size in points (default %(default)s)")
    look.add_argument("--linewidth", type=float, default=d.linewidth, help="band line width in points (default %(default)s)")
    look.add_argument("--marker-scale", type=float, default=d.marker_scale,
                      help="fat-band marker area in points squared per unit weight (default %(default)s)")
    look.add_argument("--alpha", type=float, default=d.alpha, help="fat-band marker opacity, 0 to 1 (default %(default)s)")
    look.add_argument("--cmap", default=d.cmap, help="Matplotlib colour map of the color mode (default %(default)s)")
    look.add_argument("--cmap-signed", default=d.cmap_signed,
                      help="diverging colour map for mx, my, mz; its middle is drawn grey (default %(default)s)")
    look.add_argument("--cmin", type=float, default=d.cmin, help="lower colour limit (default 0)")
    look.add_argument("--cmax", type=float, default=d.cmax, help="upper colour limit (default: largest weight in the window)")
    look.add_argument("--title", default=d.title, help="figure title (default: the selection)")
    look.add_argument("--no-title", action="store_true", help="draw no title")
    look.add_argument("--legend", action="store_true", help="draw a legend naming the selection (always on for overlays)")
    look.add_argument("--style", choices=tuple(STYLES), default=d.style,
                      help="paper (default) or talk, which scales fonts by 1.4, lines by 1.5 and marker area by 2")
    look.add_argument("--rcparams-file", default=d.rcparams_file, metavar="JSON",
                      help="JSON object of Matplotlib rcParams applied to the figure, for example {\"font.family\": \"serif\"}")
    look.add_argument("--despine", action="store_true", help="remove the top and right frame lines")
    look.add_argument("--rasterize", action="store_true",
                      help="embed markers and band lines of SVG and PDF output as an image at --dpi (default: all vector)")
    fast = parser.add_argument_group("large files")
    fast.add_argument("--kpoint-stride", type=int, default=1, metavar="N",
                      help="plot every Nth k-point of each segment, always keeping both segment ends (default 1)")
    fast.add_argument("--float32", action="store_true", help="store projections in single precision, half the memory")
    return parser


# PlotOptions from the parsed flags: every PlotOptions field has a flag of the same name, except those marked flag False.
def options_from_args(args: argparse.Namespace) -> PlotOptions:
    return PlotOptions(**{f.name: getattr(args, f.name) for f in dataclasses.fields(PlotOptions) if f.metadata.get("flag", True)})


# Output path: --output as given, or a default name in the current folder (_check_flags refuses the vasprun's folder).
def _output_path(args: argparse.Namespace, label: str, block: str = "") -> Path:
    if args.output:
        return Path(args.output).expanduser()
    return Path.cwd() / f"bands_{slug(label)}{'_' + block if block else ''}_{args.mode}{'_norm' if args.normalize else ''}.png"


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


# Refuse flag combinations that need no file to judge, before a large file is parsed; channel rules come from COMPATIBILITY.
def _check_flags(args: argparse.Namespace) -> None:
    vasprun = Path(args.vasprun).expanduser()
    if not args.output and Path.cwd().resolve() == vasprun.resolve().parent:
        raise UsageError(f"the default output would be written next to {vasprun.name}, in its data folder; "
                         "pass --output to choose the file, for example --output bands.png")
    if args.kpoint_stride < 1:
        raise UsageError(f"--kpoint-stride must be 1 or more, not {args.kpoint_stride}")
    if args.label_series and not args.overlay:
        raise UsageError("--label-series names the series of an overlay; add --overlay")
    if args.overlay and (args.atoms or args.element or args.orbitals):
        raise UsageError("with --overlay, give each series as --select ATOMS:ORBITALS; --atoms, --element and --orbitals "
                         "are not used")
    if args.overlay and not args.select:
        raise UsageError("--overlay needs --select ATOMS:ORBITALS terms, one per series, for example "
                         "--select Mn:d --select Te:p")
    if args.overlay and len(args.select) > MAX_SERIES:
        raise UsageError(f"--overlay draws at most {MAX_SERIES} series, but {len(args.select)} --select terms were given")
    if args.label_series and len(args.label_series) != len(args.select or []):
        raise UsageError(f"{len(args.label_series)} --label-series names for {len(args.select or [])} --select terms; "
                         "give one name per term")
    if args.no_labels and args.label_map:
        raise UsageError("--no-labels removes every label, so --label-map has nothing to rename; use one of them")
    for flag, given in (("--normalize", args.normalize), ("--overlay", args.overlay)):
        refusal = flag_refusal(args.component, flag) if given and args.component else None
        if refusal:
            raise UsageError(refusal)


# Selections to draw: one per --select term with --overlay, otherwise the single summed selection.
def _selections(args: argparse.Namespace, data: VasprunData) -> list[Selection]:
    if args.overlay:
        return build_series(data, args.select, args.label_series)
    return [build_selection(data, atoms=args.atoms, elements=args.element, orbitals=args.orbitals, select=args.select)]


# Panels to draw, already cut to the band window and the kept k-points, with energies shifted once to plotted units.
def _panels(data: VasprunData, choice: BlockChoice, selections: Sequence[Selection], overlay: bool, normalize: bool,
            kept: np.ndarray | None, shift: float, emin: float, emax: float) -> tuple[list["Panel"], slice]:
    """Bands are kept when they reach the window in some panel, judged on the strided view, so the weight sums never see the rest."""
    from .plotting import MARKER_COLOUR, Panel, Series
    # Series colours only for an overlay; a single selection keeps the plain marker colour.
    colours = SERIES_COLOURS if overlay else (MARKER_COLOUR,)
    rows = slice(None) if kept is None else kept
    views = [data.eigenvalues[choice.eigen(b)][rows] for b in choice.blocks]
    window = bands_reaching(views, emin + shift, emax + shift)
    panels = []
    for b, view in zip(choice.blocks, views):
        total = state_totals(data, b.projection, kept, window) if normalize else None
        series = [Series(select_weights(data, sel, b.projection, normalize=normalize, kpoints=kept, bands=window,
                                        total=total), sel.label, colours[i])
                  for i, sel in enumerate(selections)]
        panels.append(Panel(view[:, window] - shift, series, b.heading))
    return panels, window


# Run a plot: parse, choose the blocks, select, build the k-path, draw and save; print the output path and every note.
def run_plot(args: argparse.Namespace) -> int:
    vasprun = Path(args.vasprun).expanduser()
    if args.list_atoms or args.list_orbitals:
        _print_lists(args, vasprun)
        return 0
    _check_flags(args)
    from .kpath import build_kpath, parse_label_map, rename_labels, stride_indices, subsample, without_labels
    from .plotting import plot_panels, resolve_window
    options = options_from_args(args)
    options.validate(Path(args.output) if args.output else None)
    label_map = parse_label_map(args.label_map) if args.label_map else None
    data = parse_vasprun(vasprun, float32=args.float32)
    if args.efermi is not None:
        data = dataclasses.replace(data, efermi=args.efermi)
    choice = choose_blocks(data, spin=args.spin, component=args.component, normalize=args.normalize, overlay=args.overlay)
    selections = _selections(args, data)
    efermi, emin, emax = resolve_window(data, options.emin, options.emax, options.absolute)
    kpath = build_kpath(data, labels=args.high_symmetry)
    if label_map:
        kpath = rename_labels(kpath, label_map)
    if args.no_labels:
        kpath = without_labels(kpath)
    kept = None
    if args.kpoint_stride != 1:
        kept = stride_indices(data.nkpoints, kpath.segment_starts, args.kpoint_stride)
        kpath = subsample(kpath, kept)
    shift = 0.0 if options.absolute else efermi
    panels, window = _panels(data, choice, selections, args.overlay, options.normalize, kept, shift, emin, emax)
    label = " vs ".join(sel.label for sel in selections) if args.overlay else selections[0].label
    output = _output_path(args, label, choice.tag)
    title = None if args.no_title else (options.title or f"{label} ({choice.label}"
                                        f"{', normalised' if options.normalize else ''})")
    options = dataclasses.replace(options, emin=emin, emax=emax, title=title)
    notes = [n for n in (choice.note, choice.caption) if n]
    path = plot_panels(panels, kpath, output, options, zero=efermi - shift, notes=notes, signed=choice.signed,
                       quantity=choice.quantity, caption=choice.caption)
    print(f"wrote {path}")
    for sel in selections:
        nions, norbitals = len({ion for ion, _ in sel.pairs}), len({orb for _, orb in sel.pairs})
        print(f"selection: {sel.label}; {nions} ions, {norbitals} orbitals, block {choice.label}")
    print(f"Fermi energy {efermi:.6f} eV; k-path segments from {kpath.segment_source}; labels from {kpath.label_source}")
    if kept is not None:
        print(f"k-points: {kept.size} of {data.nkpoints} kept (--kpoint-stride {args.kpoint_stride}, segment ends kept)")
    print(f"bands: {window.stop - window.start} of {data.nbands} reach the window"
          + (f" (bands {window.start + 1} to {window.stop})" if window.stop > window.start else ""))
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
