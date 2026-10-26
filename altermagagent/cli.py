"""Command line entry point: altermag import | rescan | delete | summary | graph | suggest | literature | bands | structure | benchmark | parent | gui."""

import json

import argparse
import sys

from . import db
from .files import display_path


# Import a root folder as a project and print a one-line result.
def cmd_import(args):
    from .importer import import_project
    name, n = import_project(args.root, args.name)
    print(f"Imported {n} runs into project {name} (index {display_path(db.DB_PATH)}).")
    print(scan_line(db.get_project(name)))


# Re-import a project from its stored folder under the same name, rebuilding runs and summary with the current code.
def cmd_rescan(args):
    from .importer import rescan_project
    try:
        name, n = rescan_project(args.project)
    except (ValueError, FileNotFoundError) as exc:
        print(f"Rescan failed: {exc}; the index was not changed.", file=sys.stderr)
        return 1
    print(f"Rescanned {name}: {n} runs, summary rebuilt in {display_path(db.DB_PATH)}; {db.NO_DATA_CHANGE}.")
    return 0


# Remove a project and its rows from the index only; folders on disk are not touched.
def cmd_delete(args):
    n = db.delete_project(args.project)
    if n is None:
        print(f"Delete failed: no project named {args.project}; {db.NO_DATA_CHANGE}.", file=sys.stderr)
        return 1
    print(f"Removed {args.project} ({n} runs) from {display_path(db.DB_PATH)}; {db.NO_DATA_CHANGE}.")
    return 0


# Return the closing line with runs, files and seconds from the project details.
def scan_line(p):
    d = p["details"]
    return (f"Scanned {d.get('n_runs', 0)} runs and {d.get('n_files', 0)} files in "
            f"{d.get('scan_seconds', 0)} seconds; the folder was only read.")


# Print the eight summary points, the branches, the literature section (or why it is unavailable) and the runs in order.
def cmd_summary(args):
    from .literature import section, section_lines
    from .summary import as_text
    p = db.get_project(args.project)
    if p is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    runs = db.get_runs(args.project)
    brs = p["summary"].get("branches") or []
    lit = section_lines(section(p["summary"], runs, branch=brs[0] if brs else None))
    title, _, rest = as_text(args.project, p["summary"], runs, literature=lit, verbose=args.verbose).partition("\n")
    print(f"{title}\n{scan_line(p)}\n{rest}")
    return 0


# Print only the literature section of a project: published papers found by search, not evidence from its files.
def cmd_literature(args):
    from .literature import section, section_lines
    p = db.get_project(args.project)
    if p is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    brs = p["summary"].get("branches") or []
    print("\n".join(section_lines(section(p["summary"], db.get_runs(args.project), branch=brs[0] if brs else None))))
    return 0


# Print the branches of a project and their runs in historical order, as text.
def cmd_graph(args):
    from .branches import graph_text
    p = db.get_project(args.project)
    if p is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    print(graph_text(args.project, p["summary"].get("branches") or [], db.get_runs(args.project)))
    return 0


# Print the suggestions for the next DFT calculations of one branch, rebuilt from the index.
def cmd_suggest(args):
    from .branches import find_branch, suggestions, suggestions_text
    p = db.get_project(args.project)
    if p is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    brs = p["summary"].get("branches") or []
    b = find_branch(brs, args.branch)
    if b is None:
        print(f"No branch {args.branch} in project {args.project}. Branches: {', '.join(x['name'] for x in brs)}", file=sys.stderr)
        return 1
    runs = db.get_runs(args.project)
    sugg = suggestions(b, runs, project=args.project, root=p["root"], logs=p["details"].get("msg_logs"),
                       compute_band_match=True)
    print(suggestions_text(args.project, b, sugg))
    from .literature import section, section_lines
    print("\n" + "\n".join(section_lines(section(p["summary"], runs, branch=dict(b, suggestions=sugg)))))
    return 0


# Find one run of an indexed project by name, or None.
def find_run(project, run_name):
    return next((r for r in db.get_runs(project) if r["name"].rstrip("/") == run_name.rstrip("/")), None)


# Draw the band figure of one run and print the result dict and the output paths.
def cmd_bands(args):
    from .bands import plot_bands
    if db.get_project(args.project) is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    r = find_run(args.project, args.run_name)
    if r is None:
        names = ", ".join(x["name"] for x in db.get_runs(args.project))
        print(f"No run {args.run_name} in project {args.project}. Runs: {names}", file=sys.stderr)
        return 1
    try:
        res = plot_bands(r["path"], emin=args.emin, emax=args.emax, style=args.style, out_dir=args.out_dir,
                         project=args.project, run_name=r["name"])
    except (ValueError, FileNotFoundError) as exc:
        print(f"Band figure failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(res, indent=1, ensure_ascii=False))
    for kind, path in res["paths"].items():
        print(f"{kind.upper()}: {display_path(path)}")
    return 0


# Draw the structure figure with moment arrows of one run and print the moments and the output paths.
def cmd_structure(args):
    from .structure import plot_structure
    if db.get_project(args.project) is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    r = find_run(args.project, args.run_name)
    if r is None:
        names = ", ".join(x["name"] for x in db.get_runs(args.project))
        print(f"No run {args.run_name} in project {args.project}. Runs: {names}", file=sys.stderr)
        return 1
    try:
        res = plot_structure(r["path"], out_dir=args.out_dir, supercell=tuple(args.supercell),
                             project=args.project, run_name=r["name"])
    except (ValueError, FileNotFoundError) as exc:
        print(f"Structure figure failed: {exc}", file=sys.stderr)
        return 1
    print(f"{res['formula']}, {res['spacegroup']} ({res['spacegroup_number']}), from {res['structure_file']}")
    for m in res["moments"]:
        print(f"  {m['element']}{m['site']}: MAGMOM {m['magmom']}, vector {m['vector']}, size {m['size']:.3f} μB")
    print(f"  net: {res['net_moment']}, size {res['net_moment_size']:.3f} μB")
    for n in res["notes"]:
        print(f"Note: {n}")
    for kind, path in res["paths"].items():
        print(f"{kind.upper()}: {display_path(path)}")
    return 0


# Compare a project with a benchmark JSON file and print pass, fail or unknown per row; exit 1 when a row fails.
def cmd_benchmark(args):
    from .benchmark import as_text, run_benchmark
    p = db.get_project(args.project)
    if p is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    res = run_benchmark(args.benchmark, p["summary"], db.get_runs(args.project))
    print(f"Benchmark {args.benchmark} on project {args.project}")
    print(as_text(res))
    return 1 if any(r["result"] == "fail" for r in res) else 0


# Set, clear or mark unknown the parent of one run; the choice is stored in the index only.
def cmd_parent(args):
    if find_run(args.project, args.run_name) is None:
        print(f"No run {args.run_name} in project {args.project}.", file=sys.stderr)
        return 1
    val = None if args.parent == "auto" else "" if args.parent == "unknown" else args.parent
    if val and find_run(args.project, val) is None:
        print(f"No run {val} in project {args.project}.", file=sys.stderr)
        return 1
    db.set_parent_choice(args.project, args.run_name, val)
    print(f"Parent of {args.run_name} set to {args.parent} in {display_path(db.DB_PATH)}; {db.NO_DATA_CHANGE}.")
    return 0


# Start the browser interface.
def cmd_gui(args):
    from .app import run
    run(port=args.port)


# Compare a Wannier Hamiltonian with its ancestor band run and print the errors and the figure path.
def cmd_compare_wannier(args):
    from .wannier_compare import run_comparison
    if db.get_project(args.project) is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    try:
        res = run_comparison(args.project, hr_file=args.hr_file, run_name=args.run_name, reference_run=args.reference,
                             emin=args.emin, emax=args.emax, window=args.window, band_start=args.band_start,
                             structure_only=args.structure_only, use_ws_distance=False if args.no_ws_distance else None,
                             fermi_eV=args.fermi_eV, spin=args.spin, out_dir=args.out_dir)
    except (ValueError, OSError) as exc:
        print(f"Wannier comparison failed: {exc}", file=sys.stderr)
        return 1
    m = res["metrics"]
    print(f'{res["wannier_run"]} ({res["hr_kind"]} Hamiltonian {res["hamiltonian"]}) versus {res["reference_run"]} '
          f'({res["reference_source"]}), Fermi level {res["fermi_eV"]:.4f} eV from {res["fermi_source"]}')
    if m:
        print(f'MAE {m["mae_meV"]:.3f} meV, RMSE {m["rmse_meV"]:.3f} meV, 95th percentile {m["p95_meV"]:.3f} meV, '
              f'max {m["max_meV"]:.3f} meV, {m["samples_over_20_meV"]} of {m["samples"]} samples above 20 meV in the window '
              f'{res["metric_emin"]:g} to {res["metric_emax"]:g} eV ({m["mapping"]}; {m["mapping_source"]})')
    else:
        print("No errors computed: give --band-start N when the band correspondence is known, or --band-start auto.")
    if res["interpolation"].get("convention_note"):
        print(f'Interpolation: {res["interpolation"]["convention_note"]}')
    print(f'Evidence: {res["evidence_state"]}')
    for note in res["fallbacks"]:
        print(f"Note: {note}")
    print(f'Figure: {display_path(res["paths"]["png"])}')
    print(f'Details: {display_path(res["paths"]["json"])}')
    return 0


# Spin projection of a Wannier Hamiltonian along the band k-path; prints the figure, the spin order and the window means.
def cmd_spin_projection(args):
    from .spin_projection import spin_projection, summary_lines
    if db.get_project(args.project) is None:
        print(f"No project named {args.project}. Run: altermag import <root> --name {args.project}", file=sys.stderr)
        return 1
    try:
        res = spin_projection(args.project, hr_file=args.hr_file, run_name=args.run_name, reference_run=args.reference,
                              emin=args.emin, emax=args.emax, component=args.component, order=args.order,
                              out_dir=args.out_dir, fermi_eV=args.fermi_eV)
    except (ValueError, OSError) as exc:
        print(f"Spin projection failed: {exc}", file=sys.stderr)
        return 1
    print(f'Figure {display_path(res["paths"]["png"])}, data {display_path(res["paths"]["json"])}')
    print("\n".join(summary_lines(res)))
    return 0


# Parse arguments and dispatch to the sub-command.
def main(argv=None):
    ap = argparse.ArgumentParser(prog="altermag", description="Read-only storylines of VASP projects.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("import", help="scan a folder tree read-only and index it")
    p.add_argument("root")
    p.add_argument("--name", help="project name, default the folder name")
    p.set_defaults(func=cmd_import)
    p = sub.add_parser("rescan", help="re-import a project from its stored folder under the same name")
    p.add_argument("project")
    p.set_defaults(func=cmd_rescan)
    p = sub.add_parser("delete", help="remove a project from the index; folders on disk are not touched")
    p.add_argument("project")
    p.set_defaults(func=cmd_delete)
    p = sub.add_parser("summary", help="print the summary and runs of a project")
    p.add_argument("project")
    p.add_argument("--verbose", action="store_true", help="add the reason for each run's parent on an indented line")
    p.set_defaults(func=cmd_summary)
    p = sub.add_parser("graph", help="print the branches (structure, magnetic realization, SOC) and their runs")
    p.add_argument("project")
    p.set_defaults(func=cmd_graph)
    p = sub.add_parser("suggest", help="print the suggestions for the next DFT calculations of one branch")
    p.add_argument("project")
    p.add_argument("branch")
    p.set_defaults(func=cmd_suggest)
    p = sub.add_parser("literature", help="print published papers found by a Valency search (not evidence from the files)")
    p.add_argument("project")
    p.set_defaults(func=cmd_literature)
    p = sub.add_parser("bands", help="draw the band structure figure of a band run (PNG, SVG, HTML)")
    p.add_argument("project")
    p.add_argument("run_name")
    p.add_argument("--emin", type=float, default=-3.0, help="lowest energy relative to the Fermi level, eV")
    p.add_argument("--emax", type=float, default=3.0, help="highest energy relative to the Fermi level, eV")
    p.add_argument("--style", default="paper", choices=["paper", "screen"])
    p.add_argument("--out-dir", help="output folder, default ~/.altermagagent/figures/<project>")
    p.set_defaults(func=cmd_bands)
    p = sub.add_parser("structure", help="draw the crystal structure with MAGMOM arrows of a run (PNG, SVG, HTML)")
    p.add_argument("project")
    p.add_argument("run_name")
    p.add_argument("--supercell", type=int, nargs=3, default=[1, 1, 1], metavar=("NA", "NB", "NC"))
    p.add_argument("--out-dir", help="output folder, default ~/.altermagagent/figures/<project>")
    p.set_defaults(func=cmd_structure)
    p = sub.add_parser("benchmark", help="compare a project with a benchmark JSON file (pass, fail, unknown per row)")
    p.add_argument("project")
    p.add_argument("benchmark")
    p.set_defaults(func=cmd_benchmark)
    p = sub.add_parser("parent", help="set the parent of a run in the index: a run name, unknown, or auto")
    p.add_argument("project")
    p.add_argument("run_name")
    p.add_argument("parent")
    p.set_defaults(func=cmd_parent)
    p = sub.add_parser("compare-wannier", aliases=["wannier-compare"],
                       help="overlay a Wannier Hamiltonian on its ancestor VASP bands and print MAE, RMSE and max error")
    p.add_argument("project")
    p.add_argument("hr_file", nargs="?", help="path to wannier90_hr.dat or a variant such as wannier_sym/wannier90_hr.dat_*")
    p.add_argument("--run-name", "--run", dest="run_name", help="indexed Wannier run that owns wannier90_hr.dat")
    p.add_argument("--hr-file", dest="hr_file_opt", help="same as the positional hr_file")
    p.add_argument("--reference-run", "--reference", dest="reference", help="ancestor band run to compare with")
    p.add_argument("--emin", type=float, default=-3.0, help="lowest plotted energy relative to the Fermi level, eV")
    p.add_argument("--emax", type=float, default=3.0, help="highest plotted energy relative to the Fermi level, eV")
    p.add_argument("--window", type=float, default=1.0, help="errors are computed within -window to +window eV")
    p.add_argument("--band-start", default=None,
                   help="VASP band N matched to Wannier band 1 (errors only then), or auto: offset by smallest error, not verified")
    p.add_argument("--fermi-eV", dest="fermi_eV", type=float, help="known common Fermi energy, eV")
    p.add_argument("--spin", choices=["up", "down"], help="spin channel of a spin-polarised reference")
    p.add_argument("--out-dir", help="output folder outside all indexed roots, default ~/.altermagagent/figures/<project>")
    p.add_argument("--structure-only", action="store_true",
                   help="explicit preliminary mode: allow INCAR differences and missing charge files")
    p.add_argument("--no-ws-distance", action="store_true", help="direct-HR convention, without wsvec")
    p.set_defaults(func=cmd_compare_wannier)
    p = sub.add_parser("spin-projection", help="colour Wannier bands by <sigma_x>, <sigma_y>, <sigma_z> along the band k-path")
    p.add_argument("project")
    p.add_argument("hr_file", nargs="?", help="path to wannier90_hr.dat or a variant such as wannier_sym/wannier90_hr.dat_*")
    p.add_argument("--run-name", "--run", dest="run_name", help="indexed Wannier run that owns wannier90_hr.dat")
    p.add_argument("--reference-run", "--reference", dest="reference", help="ancestor band run whose k-path is used")
    p.add_argument("--emin", type=float, default=-1.0, help="lowest energy relative to the Fermi level, eV")
    p.add_argument("--emax", type=float, default=1.0, help="highest energy relative to the Fermi level, eV")
    p.add_argument("--component", choices=["x", "y", "z", "all"], default="z", help="colouring component; all draws three panels")
    p.add_argument("--order", choices=["udud", "uudd"],
                   help="Wannier spin order; default from the file header or name, else udud marked as assumed")
    p.add_argument("--fermi-eV", dest="fermi_eV", type=float, help="known common Fermi energy, eV")
    p.add_argument("--out-dir", help="output folder outside all indexed roots, default ~/.altermagagent/figures/<project>")
    p.set_defaults(func=cmd_spin_projection)
    p = sub.add_parser("gui", help="start the browser interface")
    p.add_argument("--port", type=int, default=8080)
    p.set_defaults(func=cmd_gui)
    args = ap.parse_args(argv)
    if getattr(args, "hr_file_opt", None):
        args.hr_file = args.hr_file or args.hr_file_opt
    if getattr(args, "band_start", "auto") not in ("auto", None):
        args.band_start = int(args.band_start)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
