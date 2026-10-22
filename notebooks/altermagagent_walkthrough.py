import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


# Title and purpose of the notebook.
@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # AltermagAgent walkthrough

        This notebook reproduces the AltermagAgent demo and its analysis in one file.
        It reads an indexed tree of VASP calculation folders and never writes inside it.
        It shows the summary points, the per-branch suggestions, the runs in order, the band and structure figures and the benchmark.
        No API key is needed, because the chat is not part of this notebook.
        """
    )
    return (mo,)


# Paths of the repository, the demo tree and the benchmark files.
@app.cell
def _():
    from pathlib import Path

    from altermagagent import db

    REPO = Path(__file__).resolve().parent.parent
    DEMO_ROOT = REPO / "demo" / "CrSe_demo"
    DEMO_BENCH = REPO / "demo" / "benchmark_crse.json"
    MNTE_BENCH = REPO / "docs" / "benchmark_mnte_template.json"
    return DEMO_BENCH, DEMO_ROOT, MNTE_BENCH, db


# Choose the project to show; the list comes from the index.
@app.cell
def _(db, index_note, mo):
    _names = [p["name"] for p in db.list_projects()]
    _default = "MnTe" if "MnTe" in _names else ("CrSe" if "CrSe" in _names else _names[0])
    project_pick = mo.ui.dropdown(options=_names, value=_default, label="Project")
    mo.vstack([mo.md("Pick a project from the index in `~/.altermagagent/index.sqlite`; every cell below follows it."),
               project_pick, mo.md(index_note)])
    return (project_pick,)


# Import the synthetic CrSe demo tree when the index is empty.
@app.cell
def _(DEMO_ROOT, db, mo):
    from altermagagent.importer import import_project

    if db.list_projects():
        index_note = "The index already holds projects, so no import was needed."
    else:
        if not DEMO_ROOT.is_dir():
            import subprocess
            import sys

            subprocess.run([sys.executable, str(DEMO_ROOT.parent / "make_demo.py")], check=True)
        _name, _n = import_project(DEMO_ROOT, name="CrSe")
        index_note = f"The index was empty, so `{DEMO_ROOT}` was imported read-only as project CrSe with {_n} runs."
    mo.md("This cell imports `demo/CrSe_demo` read-only when the index is empty.")
    return (index_note,)


# Load the summary and runs of the chosen project.
@app.cell
def _(db, project_pick):
    project = project_pick.value
    summ = db.get_project(project)["summary"]
    runs = db.get_runs(project)
    return project, runs, summ


# The nine summary points with the evidence tag of each point.
@app.cell
def _(mo, project, summ):
    _rows = [{"#": i, "point": p["label"], "value": p["value"], "evidence tag": p.get("state", ""),
              "evidence": p.get("state_note", "")} for i, p in enumerate(summ["points"], 1)]
    mo.vstack([mo.md(f"## Nine summary points of {project}\n\nEach point carries its evidence tag and the file it was read from."),
               mo.ui.table(_rows, selection=None, pagination=False)])
    return


# The per-branch suggestions for the next DFT calculations, six steps per branch.
@app.cell
def _(mo, summ):
    _parts = ["## Suggestions for the next DFT calculations, per branch",
              "A branch is one structure, one magnetic realization and one SOC setting."]
    for _b in summ.get("branches") or []:
        _s = _b.get("suggestions") or {}
        _parts.append(f"### {_b['name']} ({_b['realization']}, {'SOC' if _b['soc'] else 'non-SOC'}, {len(_b['runs'])} runs)")
        _parts.append("\n".join(f"- ({_t['step']}) {_t['title']}: {_t['state']}" for _t in _s.get("steps", [])))
        _parts.append(f"**{_s.get('next', '')}**")
    mo.md("\n\n".join(_parts))
    return


# The runs in historical sequence with phase, status and parent note.
@app.cell
def _(mo, runs):
    from altermagagent.summary import STATUS_TEXT, parent_text, phase_text

    _rows = [{"#": r["seq"], "run": r["name"], "phase": phase_text(r["details"]), "status": STATUS_TEXT[r["status"]],
              "parent": parent_text(r)} for r in runs]
    mo.vstack([mo.md("## Runs in historical sequence\n\nEach run shows its phase, its status and how it came from its parent."),
               mo.ui.table(_rows, selection=None, pagination=False)])
    return


# Controls for the band figure: the band run and the energy window.
@app.cell
def _(mo, runs):
    _bands = [r["name"] for r in runs if r["details"].get("phase") == "band"]
    _default = "bands_nosoc" if "bands_nosoc" in _bands else (_bands[0] if _bands else None)
    band_pick = mo.ui.dropdown(options=_bands, value=_default, label="Band run")
    emin = mo.ui.slider(start=-8.0, stop=-0.5, step=0.5, value=-3.0, label="emin (eV)", show_value=True)
    emax = mo.ui.slider(start=0.5, stop=8.0, step=0.5, value=3.0, label="emax (eV)", show_value=True)
    mo.vstack([mo.md("## Band figure\n\nMove the sliders to redraw the bands of the chosen run relative to the Fermi level."),
               mo.hstack([band_pick, emin, emax], justify="start")])
    return band_pick, emax, emin


# Draw the band figure with plot_bands and show the PNG.
@app.cell
def _(band_pick, emax, emin, mo, project, runs):
    from altermagagent.bands import plot_bands

    _run = next((r for r in runs if r["name"] == band_pick.value), None)
    if _run is None:
        _out = mo.md("This project has no band run.")
    else:
        try:
            _res = plot_bands(_run["path"], emin=emin.value, emax=emax.value, project=project, run_name=_run["name"])
            _split = _res["largest_splitting_eV"]
            _line = f"Largest spin splitting in the window: {_split:.2f} eV." if _split is not None else "No spin splitting measured."
            _out = mo.vstack([mo.image(src=_res["paths"]["png"]),
                              mo.md(f"{_line} Source `{_res['source']}`.\n\n" + "\n".join(f"- {n}" for n in _res["notes"]))])
        except Exception as _e:
            _out = mo.md(f"The band figure of {_run['name']} could not be drawn: {_e}")
    _out
    return


# Draw the structure figure of the first run with a structure file.
@app.cell
def _(mo, project, runs):
    from altermagagent.structure import plot_structure

    _run = next((r for r in runs if r["details"].get("structure_file")), None)
    if _run is None:
        _out = mo.md("No run of this project has a POSCAR or CONTCAR.")
    else:
        try:
            _res = plot_structure(_run["path"], project=project, run_name=_run["name"])
            _out = mo.vstack([mo.md(f"## Structure\n\nAtoms of {_run['name']} with the INCAR starting moments as arrows: "
                                    f"{_res['formula']}, {_res['spacegroup']} ({_res['spacegroup_number']})."),
                              mo.image(src=_res["paths"]["png"])])
        except Exception as _e:
            _out = mo.md(f"The structure figure of {_run['name']} could not be drawn: {_e}")
    _out
    return


# The benchmark table with pass, fail or unknown per row.
@app.cell
def _(DEMO_BENCH, DEMO_ROOT, MNTE_BENCH, db, mo, project, runs, summ):
    from altermagagent.benchmark import fmt, run_benchmark

    _root = db.get_project(project)["root"]
    _file = DEMO_BENCH if _root == str(DEMO_ROOT.resolve()) else (MNTE_BENCH if project == "MnTe" else None)
    if _file is None:
        _out = mo.md("## Benchmark\n\nNo benchmark file belongs to this project.")
    else:
        _res = run_benchmark(_file, summ, runs)
        _rows = [{"row": r["id"], "result": r["result"], "expected": fmt(r["expected"]), "actual": fmt(r["actual"]),
                  "units": r["units"], "reason": r["reason"], "source": r["source"]} for r in _res]
        _n = {k: sum(1 for r in _res if r["result"] == k) for k in ("pass", "fail", "unknown")}
        _out = mo.vstack([mo.md(f"## Benchmark\n\nThe project compared with `{_file.name}`: "
                                f"{_n['pass']} pass, {_n['fail']} fail, {_n['unknown']} unknown."),
                          mo.ui.table(_rows, selection=None, pagination=False)])
    _out
    return


# The recommended next calculation of each branch, one sentence each.
@app.cell
def _(mo, summ):
    _lines = [f"- {_b['name']}: {(_b.get('suggestions') or {}).get('next', '')}" for _b in summ.get("branches") or []]
    mo.md("## Next calculation per branch\n\n" + ("\n".join(_lines) or "No branches."))
    return


if __name__ == "__main__":
    app.run()
