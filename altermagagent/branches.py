"""Branches of a project: runs grouped by structure, magnetic realization and SOC, the git-like graph and the per-branch suggestions."""

import math
import re
from pathlib import Path

from .files import display_path

MOMENT_TOL_UB = 0.1          # moments below this count as zero, and patterns match within this, Bohr magneton
PARALLEL_TOL_DEG = 1.0       # largest angle from one common axis for collinear moments
PLANE_TOL_UB = 0.1           # largest out-of-plane component for coplanar moments, Bohr magneton
NET_TOL_UB = 0.01            # declared tolerance on the net moment per cell, Bohr magneton
CLASSES = ("collinear", "coplanar", "non-coplanar", "non-magnetic", "unknown")
STEP_TITLES = {
    "a": "Crystal structure and magnetic moments",
    "b": "Symmetries",
    "c": "Converged SCF",
    "d": "Band structure and band projections",
    "e": "Wannierisation",
    "f": "Transport",
}
DONE, PARTIAL, MISSING, NA = "done", "partial", "missing", "not applicable"
PHASE_COLOURS = {"relax": "#a78bfa", "scf": "#34d399", "band": "#4f8cff", "wannier": "#fbbf24",
                 "symmetrisation": "#f472b6", "transport": "#f87171", "projection": "#22d3ee", "other": "#a1a1aa"}


# Strict structure matcher used for branches and parent links: no volume rescale, 1 percent lattice tolerance.
def matcher():
    from pymatgen.analysis.structure_matcher import StructureMatcher
    return StructureMatcher(ltol=0.01, scale=False)


# Size of a 3-vector.
def norm(v):
    return math.sqrt(sum(x * x for x in v))


# Dot product of two 3-vectors.
def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


# Cross product of two 3-vectors.
def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


# True when a run is non-collinear (LSORBIT or LNONCOLLINEAR set in INCAR).
def noncollinear(d):
    raw = str(d.get("incar_raw", {}).get("LNONCOLLINEAR", "")).strip(". ").upper()
    return d.get("soc") or raw in ("T", "TRUE")


# Moment vectors of one run as (vectors, source text): INCAR MAGMOM first, converged OUTCAR moments otherwise.
def moment_vectors(run):
    d = run["details"]
    nat = len(d.get("species") or [])
    vals = d.get("magmom_list") or []
    if nat and vals:
        if len(vals) == 3 * nat:
            return [vals[3 * i:3 * i + 3] for i in range(nat)], "INCAR MAGMOM"
        if len(vals) == nat:
            return [[0.0, 0.0, v] for v in vals], "INCAR MAGMOM"
    mom = d.get("moments")
    if mom and mom.get("sites") and run.get("status") == "completed":
        vecs = [s.get("vector") or [0.0, 0.0, s["moment"]] for s in mom["sites"]]
        return vecs, "OUTCAR magnetization (converged run)"
    if nat and d.get("incar", {}).get("ISPIN") in (None, 1) and not noncollinear(d):
        return [[0.0, 0.0, 0.0] for _ in range(nat)], "ISPIN = 1, no MAGMOM"
    return None, "no MAGMOM in INCAR and no converged OUTCAR moments"


# Classify moment vectors as collinear, coplanar, non-coplanar or non-magnetic.
def classify(vecs):
    if vecs is None:
        return "unknown"
    mag = [v for v in vecs if norm(v) >= MOMENT_TOL_UB]
    if not mag:
        return "non-magnetic"
    ref = mag[0]
    off = [v for v in mag if abs(dot(ref, v)) / (norm(ref) * norm(v)) < math.cos(math.radians(PARALLEL_TOL_DEG))]
    if not off:
        return "collinear"
    n = cross(ref, off[0])
    n = [x / norm(n) for x in n]
    if all(abs(dot(n, v)) <= PLANE_TOL_UB for v in mag):
        return "coplanar"
    return "non-coplanar"


# True when two moment patterns are equal up to one global rotation, within MOMENT_TOL_UB per site.
def same_pattern(a, b):
    if a is None or b is None:
        return a is None and b is None
    if len(a) != len(b):
        return False
    if any(abs(norm(x) - norm(y)) > MOMENT_TOL_UB for x, y in zip(a, b)):
        return False
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            tol = MOMENT_TOL_UB * (norm(a[i]) + norm(a[j]) + MOMENT_TOL_UB)
            if abs(dot(a[i], a[j]) - dot(b[i], b[j])) > tol:
                return False
    return True


# SOC word used in branch names.
def soc_word(soc):
    return "SOC" if soc else "noSOC"


# Short pattern text such as "Cr ±3.00 μB (6 sites)" from moment vectors and species.
def pattern_text(vecs, species):
    if vecs is None:
        return "moments unknown"
    by_el = {}
    for v, el in zip(vecs, species or ["?"] * len(vecs)):
        if norm(v) >= MOMENT_TOL_UB:
            by_el.setdefault(el, []).append(norm(v))
    if not by_el:
        return "all moments zero"
    return "; ".join(f"{el} {max(s):.2f} μB ({len(s)} sites)" for el, s in by_el.items())


# Group runs into realizations (structure, class, pattern) and branches (realization, SOC); returns the branch list.
def build_branches(runs):
    from pymatgen.core import Structure
    m = matcher()
    reals, branches, of_run = [], [], {}
    for r in runs:
        d = r["details"]
        if not d.get("structure"):
            continue
        s = Structure.from_dict(d["structure"])
        vecs, src = moment_vectors(r)
        cls = classify(vecs)
        real = next((x for x in reals if x["class"] == cls and same_pattern(x["vectors"], vecs)
                     and m.fit(x["structure"], s)), None)
        if real is None:
            real = {"index": len(reals) + 1, "class": cls, "vectors": vecs, "structure": s, "source": src,
                    "pattern": pattern_text(vecs, d.get("species")), "runs": []}
            reals.append(real)
        real["runs"].append(r["name"])
        soc = bool(d.get("soc"))
        br = next((b for b in branches if b["realization"] is real and b["soc"] == soc), None)
        if br is None:
            br = {"realization": real, "soc": soc, "runs": []}
            branches.append(br)
        br["runs"].append(r["name"])
        of_run[r["name"]] = br
    by_name = {r["name"]: r for r in runs}
    for r in runs:
        if r["name"] in of_run:
            continue
        home = of_run.get(r.get("parent") or "")
        if home is None:
            folder = str(Path(r["name"]).parent)
            home = next((of_run[e["name"]] for e in reversed(runs[:runs.index(r)])
                         if e["name"] in of_run and str(Path(e["name"]).parent) == folder), None)
        if home is None:
            home = next((b for b in branches if b.get("other")), None)
            if home is None:
                home = {"realization": None, "soc": bool(r["details"].get("soc")), "runs": [], "other": True}
                branches.append(home)
        home["runs"].append(r["name"])
        of_run[r["name"]] = home
    name_branches(branches)
    out = []
    for b in branches:
        real = b["realization"]
        names = sorted(b["runs"], key=lambda n: by_name[n]["seq"])
        out.append({"name": b["name"], "soc": b["soc"], "runs": names,
                    "realization": real["class"] if real else "unknown",
                    "realization_index": real["index"] if real else None,
                    "pattern": real["pattern"] if real else "no structure in these runs",
                    "moment_source": real["source"] if real else "",
                    "vectors": real["vectors"] if real else None})
    return out


# Give every branch a name: the shared first folder plus noSOC or SOC when obvious, otherwise realization-N.
def name_branches(branches):
    used = set()
    for b in branches:
        if b.get("other"):
            b["name"] = "other runs"
            continue
        firsts = {n.split("/")[0] for n in b["runs"]}
        deep = all("/" in n for n in b["runs"])
        rivals = [x for x in branches if x is not b and not x.get("other") and x["soc"] == b["soc"]
                  and {n.split("/")[0] for n in x["runs"]} & firsts]
        if len(firsts) == 1 and deep and not rivals:
            name = f"{firsts.pop()}-{soc_word(b['soc'])}"
        else:
            name = f"realization-{b['realization']['index']}-{soc_word(b['soc'])}"
        while name in used:
            name += "'"
        used.add(name)
        b["name"] = name


# Branch dict of one run name, or None.
def branch_of(branches, run_name):
    return next((b for b in branches if run_name in b["runs"]), None)


# Find one branch by name, accepting a case-insensitive match.
def find_branch(branches, name):
    return next((b for b in branches if b["name"] == name), None) or \
        next((b for b in branches if b["name"].lower() == (name or "").lower()), None)


# ---------- the git-like graph ----------

# Common folder prefix of the run names of one branch, such as "coplanar/1-noSOC/", or "".
def common_folder(names):
    parts = [n.split("/")[:-1] for n in names]
    if not parts:
        return ""
    pre = []
    for comp in zip(*parts):
        if len(set(comp)) != 1:
            break
        pre.append(comp[0])
    return "/".join(pre) + "/" if pre else ""


# Inline SVG of the branches: one lane per branch, nodes in time order coloured by phase, edges for parent links.
def graph_svg(branches, runs, selected=None, event=None):
    import html
    esc = html.escape
    lane_w = max(120, 8 * max((len(b["name"]) for b in branches), default=10) + 16)
    row_h, top, left = 30, 64, 8
    seq = {r["name"]: i for i, r in enumerate(runs)}
    lane = {n: i for i, b in enumerate(branches) for n in b["runs"]}
    short = {}
    for b in branches:
        pre = common_folder(b["runs"])
        for n in b["runs"]:
            short[n] = n[len(pre):] if pre else n
    click = f' onclick="emitEvent(\'{event}\', this.dataset.branch)"' if event else ""
    xs = lambda n: left + lane[n] * lane_w + lane_w // 2
    ys = lambda n: top + seq[n] * row_h
    tx = left + len(branches) * lane_w + 10
    width = tx + max(260, 8 * max((len(v) for v in short.values()), default=10) + 110)
    height = top + len(runs) * row_h + 10
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" preserveAspectRatio="xMinYMin meet" '
             f'style="width:100%;min-width:{min(width, 640)}px;height:auto;font-family:system-ui,sans-serif" '
             'role="img" aria-label="branch graph: one lane per branch, one node per run">']
    for i, b in enumerate(branches):
        x = left + i * lane_w + lane_w // 2
        sel = b["name"] == selected
        parts.append(f'<g class="glane" data-branch="{esc(b["name"])}" style="cursor:pointer"{click}>'
                     f'<rect x="{x - lane_w // 2 + 3}" y="4" width="{lane_w - 6}" height="{height - 8}" rx="8" '
                     f'fill="{"#14213d" if sel else "#0f1115"}" stroke="{"#4f8cff" if sel else "#262a33"}"/>'
                     f'<text x="{x}" y="24" text-anchor="middle" fill="#e6e8ec" font-size="13" font-weight="700">'
                     f'{esc(b["name"])}</text>'
                     f'<text x="{x}" y="42" text-anchor="middle" fill="#9aa3ad" font-size="11">'
                     f'{esc(b["realization"])}, {len(b["runs"])} run{"s" if len(b["runs"]) != 1 else ""}</text></g>')
        if b["runs"]:
            y0, y1 = ys(b["runs"][0]), ys(b["runs"][-1])
            parts.append(f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1}" stroke="#4b5260" stroke-width="2" '
                         'stroke-dasharray="3 5" pointer-events="none"/>')
    for r in runs:
        p = r.get("parent")
        if p and p in lane and r["name"] in lane:
            x0, y0, x1, y1 = xs(p), ys(p), xs(r["name"]), ys(r["name"])
            if x0 == x1:
                d = f"M{x0},{y0} L{x1},{y1}"
            else:
                d = f"M{x0},{y0} C{x0},{(y0 + y1) / 2} {x1},{(y0 + y1) / 2} {x1},{y1}"
            parts.append(f'<path d="{d}" fill="none" stroke="#c7ccd4" stroke-width="2.5" pointer-events="none"/>')
    for r in runs:
        n = r["name"]
        if n not in lane:
            continue
        ph = r["details"].get("phase") or "other"
        col = PHASE_COLOURS.get(ph, PHASE_COLOURS["other"])
        br = branches[lane[n]]["name"]
        ring = {"completed": "#e6e8ec", "error": "#f87171", "running": "#fbbf24"}.get(r["status"], "#52525b")
        label = f'{r["seq"]}. {short[n]}'
        parts.append(f'<g class="gnode" data-branch="{esc(br)}" data-run="{esc(n)}" style="cursor:pointer"{click}>'
                     f'<title>{esc(n)}: {esc(ph)}, {esc(r["status"])}, branch {esc(br)}</title>'
                     f'<line x1="{xs(n) + 10}" y1="{ys(n)}" x2="{tx - 4}" y2="{ys(n)}" stroke="#262a33" stroke-width="1"/>'
                     f'<circle cx="{xs(n)}" cy="{ys(n)}" r="9" fill="{col}" stroke="{ring}" stroke-width="2"/>'
                     f'<text x="{tx}" y="{ys(n) + 4}" fill="#e6e8ec" font-size="13">{esc(label)}</text>'
                     f'<text x="{tx + 16 + 7.5 * len(label)}" y="{ys(n) + 4}" fill="{col}" font-size="12">{esc(ph)}</text></g>')
    parts.append("</svg>")
    return "".join(parts)


# Legend markup for the phase colours of the graph.
def graph_legend():
    return " ".join(f'<span style="display:inline-flex;align-items:center;gap:4px;margin-right:10px">'
                    f'<span style="width:12px;height:12px;border-radius:50%;background:{c};display:inline-block"></span>{p}</span>'
                    for p, c in PHASE_COLOURS.items())


# Text form of the graph: each branch with its realization and its runs in historical order.
def graph_text(project, branches, runs):
    by = {r["name"]: r for r in runs}
    lines = [f"Branches of {project}: {len(branches)} branches, {len(runs)} runs"]
    for b in branches:
        lines.append("")
        lines.append(f'* {b["name"]}  ({b["realization"]}, {"SOC" if b["soc"] else "non-SOC"}, {len(b["runs"])} runs; '
                     f'{b["pattern"]})')
        for n in b["runs"]:
            r = by[n]
            par = r.get("parent")
            pb = branch_of(branches, par) if par else None
            link = f"parent {par}" + (f" [branch {pb['name']}]" if pb and pb is not b else "") if par else "no parent link"
            lines.append(f'  | {r["seq"]:>2}. {n}  {r["details"].get("phase") or "other"}  {r["status"]}  {link}')
    return "\n".join(lines)


# ---------- suggestions for the next DFT calculations ----------

# True when a path (relative to the root) lies in or below a folder (relative to the root, "." for the root).
def under(path, folder):
    return folder in (".", "") or path == folder or path.startswith(folder.rstrip("/") + "/")


# MSG evidence of one branch: findsym or MSG logs in the folders of its runs (or above them), as relative paths.
def msg_logs(branch, logs):
    out = []
    for lg in logs or []:
        folder = str(Path(lg).parent)
        if any(under(n, folder) for n in branch["runs"]):
            out.append(lg)
    return out


# Magnetic space group label and BNS number from one findsym.log through the parsers module, or (None, None).
def findsym_msg(root, rel):
    try:
        from .parsers_extra import parse_findsym_log
        res = parse_findsym_log(Path(root) / rel)
    except Exception:
        return None, None
    if not res or res.get("error") or not res.get("bns_label"):
        return None, None
    num = res.get("bns_number")
    return (f'{res["bns_label"]} (BNS {num})' if num else res["bns_label"]), (str(num) if num else None)


# Define MSG for one branch: spglib on its first converged run, the findsym.log label when one exists, and the findsym.in.
def define_msg(branch, runs, project=None, root=None, logs=None):
    by = {r["name"]: r for r in runs}
    br = [by[n] for n in branch["runs"] if n in by]
    conv = [r for r in br if r["status"] == "completed" and (r["details"].get("scf_converged") or {}).get("converged")]
    run = conv[0] if conv else next((r for r in br if r["details"].get("structure_file")), None)
    out = {"run": run["name"] if run else None, "computed": None, "bns": None, "error": None, "source": None,
           "findsym": None, "findsym_bns": None, "findsym_log": None, "agreement": None, "findsym_in": None}
    if run is not None:
        try:
            from .symmetry import PENDING, magnetic_space_group
            m = magnetic_space_group(run["path"])
            out.update(computed=f'{m.get("bns_symbol") or "symbol not in table"} (BNS {m["bns_number"]})', bns=str(m["bns_number"]),
                       source=f'spglib on {run["name"]}/{m["structure_file"]} with moments from ' +
                       (f'{run["name"]}/OUTCAR' if "OUTCAR" in m["moment_source"].split("(")[0] else f'{run["name"]}/INCAR MAGMOM'),
                       state=PENDING, n_operations=m["n_operations"], n_unitary=m["n_unitary"], n_antiunitary=m["n_antiunitary"])
        except Exception as exc:
            out["error"] = f"{type(exc).__name__}: {exc}"
    found = msg_logs(branch, logs)
    if found and root:
        out["findsym_log"] = found[0]
        out["findsym"], out["findsym_bns"] = findsym_msg(root, found[0])
        if out["findsym_bns"] and out["bns"]:
            out["agreement"] = "agree" if out["findsym_bns"] == out["bns"] else "differ"
    path = write_findsym(project, branch, runs)
    out["findsym_in"] = str(path) if path else None
    return out


# Evidence lines of step (b) from a define_msg result.
def msg_evidence(msg):
    pend = "parsed, pending validation"
    comp = (f'magnetic space group computed: {msg["computed"]} [{pend}] ({msg["source"]}; {msg["n_operations"]} operations, '
            f'{msg["n_unitary"]} unitary, {msg["n_antiunitary"]} antiunitary)' if msg["computed"] else
            f'magnetic space group computed: not computed ({msg["error"] or "no run with a structure"})')
    lines = []
    if msg["findsym"]:
        lines.append(f'magnetic space group: {msg["findsym"]} from findsym.log ({msg["findsym_log"]})')
        lines.append(comp + (f'; {msg["agreement"]} with findsym.log' if msg["agreement"] else ""))
    elif msg["findsym_log"]:
        lines.append(f'magnetic space group: findsym.log present ({msg["findsym_log"]}), no magnetic label read')
        lines.append(comp)
    else:
        lines.append(comp)
    lines.append(f'findsym input written to {display_path(msg["findsym_in"])}' if msg["findsym_in"] else "findsym input could not be written")
    return lines


# Path of the findsym input written for one branch under the index folder.
def findsym_path(project, branch_name):
    from . import db
    safe = re.sub(r"[^A-Za-z0-9._'-]+", "_", branch_name)
    return db.HOME_DIR / "msg" / project / safe / "findsym.in"


# FINDSYM input text (structure plus moment vectors in lattice-vector components) for the first run of a branch.
def findsym_input(run, title):
    from pymatgen.core import Structure
    import numpy as np
    from .structure import saxis_to_cartesian
    d = run["details"]
    s = Structure.from_dict(d["structure"])
    vecs, _ = moment_vectors(run)
    if vecs is None:
        vecs = [[0.0, 0.0, 0.0]] * len(s)
    saxis = [float(x) for x in str(d.get("incar_raw", {}).get("SAXIS", "0 0 1")).split()[:3]] or [0, 0, 1]
    if len(saxis) != 3:
        saxis = [0.0, 0.0, 1.0]
    cart = [list(saxis_to_cartesian(v, saxis)) for v in vecs]
    inv = np.linalg.inv(s.lattice.matrix)
    frac = [list(np.array(v) @ inv) for v in cart]
    big = max((abs(x) for v in frac for x in v), default=0.0) or 1.0
    a, b, c = s.lattice.abc
    al, be, ga = s.lattice.angles
    lines = [title, "0.0001", "2", f"{a:.5f} {b:.5f} {c:.5f} {al:.4f} {be:.4f} {ga:.4f}", "2", "P", str(len(s)),
             " ".join(sp.symbol for sp in s.species), "magnetic"]
    for site, f in zip(s, frac):
        x, y, z = site.frac_coords
        lines.append(f"{x:.6f} {y:.6f} {z:.6f}   " + " ".join(f"{v / big:.5f}" if abs(v / big) > 1e-6 else "0" for v in f))
    return "\n".join(lines) + "\n"


# Write the FINDSYM input of one branch under ~/.altermagagent/msg/<project>/<branch>/ and return its path, or None.
def write_findsym(project, branch, runs):
    by = {r["name"]: r for r in runs}
    first = next((by[n] for n in branch["runs"] if by[n]["details"].get("structure")), None)
    if first is None or not project:
        return None
    path = findsym_path(project, branch["name"])
    try:
        text = findsym_input(first, f"{first['details'].get('formula', 'structure')} {branch['name']} from {first['name']}")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    except Exception:
        return None
    return path


# One suggestion step dict.
def step(key, state, evidence, action=None):
    return {"step": key, "title": STEP_TITLES[key], "state": state, "evidence": evidence, "action": action}


# INCAR text of a boolean.
def tf(v):
    return ".TRUE." if v else ".FALSE."


# The six steps of the physicist's loop for one branch, each with a state and its evidence, plus the next calculation.
def suggestions(branch, runs, project=None, root=None, logs=None, compute_band_match=False):
    by = {r["name"]: r for r in runs}
    br = [by[n] for n in branch["runs"] if n in by]
    soc = branch["soc"]
    steps = []
    # (a) structure and moments
    st_runs = [r for r in br if r["details"].get("structure_file")]
    first = st_runs[0] if st_runs else None
    if first:
        ev = [f'{first["name"]}/{first["details"]["structure_file"]}: {first["details"].get("formula")}, '
              f'{first["details"].get("spacegroup")} ({first["details"].get("spacegroup_number")})']
        moms = [r for r in br if r["details"].get("moments") and r["status"] == "completed"]
        moms = [r for r in moms if r["details"].get("phase") == "scf"] or moms
        known = branch["realization"] not in ("unknown",)
        ev.append(f'moments: {branch["realization"]}, {branch["pattern"]} (from {branch["moment_source"]})')
        if moms:
            src = moms[-1]
            m = src["details"]["moments"]
            net = norm(m["net_vector"]) if m.get("net_vector") else abs(m["net"])
            comp = "compensated within" if net < NET_TOL_UB else "not compensated within"
            ev.append(f'{src["name"]}/OUTCAR: net moment {net:.2f} μB per cell, {comp} the tolerance {NET_TOL_UB} μB')
        steps.append(step("a", DONE if known else PARTIAL, ev,
                          None if known else "set MAGMOM in INCAR for this branch"))
    else:
        steps.append(step("a", MISSING, ["no POSCAR or CONTCAR in this branch"], "write POSCAR and MAGMOM"))
    # (b) symmetries
    sg = f'{first["details"].get("spacegroup")} ({first["details"].get("spacegroup_number")})' if first else None
    sev = [f"crystal space group {sg}: spglib, symprec 0.01 ({first['name']}/{first['details']['structure_file']})"] if sg else \
        ["crystal space group: no readable structure"]
    if branch["realization"] == "non-magnetic":
        steps.append(step("b", DONE if sg else MISSING, sev + ["magnetic space group: not applicable, no moments"]))
    else:
        msg = define_msg(branch, runs, project=project, root=root, logs=logs)
        path = msg["findsym_in"]
        label = msg["findsym"] or msg["computed"]
        if label:
            from .summary import PENDING_LABEL, find_validation, msg_label, validation_tag
            val = find_validation(project, "magnetic space group", branch["name"], "branch")
            src = f'from findsym.log ({msg["findsym_log"]})' if msg["findsym"] else f'computed ({msg["source"]})'
            agree = f'; computed {"agrees" if msg["agreement"] == "agree" else "differs"}' if msg["agreement"] else ""
            sev.append(f"magnetic space group: {msg_label(label)} {src}{agree}; {validation_tag(val, label, PENDING_LABEL)}")
            steps.append(step("b", DONE, sev))
        elif msg["findsym_log"]:
            sev += msg_evidence(msg)
            steps.append(step("b", DONE if sg else PARTIAL, sev))
        else:
            sev += msg_evidence(msg)
            act = (f"define MSG: validate the computed label by running findsym on {path}" if msg["computed"] and path else
                   f"define MSG: run findsym on {path} or the symmetry branch" if path else "define MSG")
            steps.append(step("b", PARTIAL if sg else MISSING, sev, act))
        steps[-1]["msg"] = msg
    # (c) converged SCF
    scfs = [r for r in br if r["details"].get("phase") == "scf"]
    conv = [r for r in scfs if r["status"] == "completed" and (r["details"].get("scf_converged") or {}).get("converged")]
    if conv:
        c = conv[-1]
        steps.append(step("c", DONE, [f'{c["name"]}: completion marker in {c["details"].get("completion_file")}; '
                                       f'{c["details"]["scf_converged"]["text"]}']))
    elif scfs:
        ev = []
        for r in scfs:
            cv = r["details"].get("scf_converged") or {}
            ev.append(f'{r["name"]}: {r["status"]}; {cv.get("text", "convergence not read")}')
        steps.append(step("c", PARTIAL, ev))
    else:
        steps.append(step("c", MISSING, ["no SCF run in this branch"]))
    # (d) bands and projections
    bands = [r for r in br if r["details"].get("phase") == "band"]
    with_vr = [r for r in bands if r["details"].get("has_vasprun")]
    if with_vr:
        b = next((r for r in reversed(with_vr) if r["status"] == "completed"), with_vr[-1])
        names = {f["name"] for f in b["details"]["files"]}
        pb = sorted(n for n in names if n.startswith("PBAND_"))
        ev = [f'{b["name"]}/vasprun.xml: band run, {b["status"]}']
        if pb:
            ev.append(f'projections: {len(pb)} PBAND files ({b["name"]}/{pb[0]}, ...)')
        elif b["details"]["incar"].get("LORBIT"):
            ev.append(f'projections: LORBIT = {b["details"]["incar"]["LORBIT"]}, projected weights in {b["name"]}/vasprun.xml')
        else:
            ev.append("projections: none (no PBAND files, LORBIT not set)")
        ok = b["status"] == "completed" and (pb or b["details"]["incar"].get("LORBIT"))
        steps.append(step("d", DONE if ok else PARTIAL, ev))
    elif bands:
        steps.append(step("d", MISSING, [f'{r["name"]}: band run without vasprun.xml ({r["status"]})' for r in bands]))
    else:
        steps.append(step("d", MISSING, ["no band run (ICHARG = 11 with a line-mode k-path) in this branch"]))
    # (e) Wannierisation
    hr = [r for r in br if r["details"].get("has_hr")]
    wan = [r for r in br if r["details"].get("phase") == "wannier" or r["details"].get("has_win")]
    if hr:
        w = hr[-1]
        sym = [e for e in w["details"].get("steps", {}).get("symmetrisation", [])]
        sub = [p for p in w["details"].get("sub_files", [])]
        names = {f["name"] for f in w["details"]["files"]}
        ev = [f'{w["name"]}/wannier90_hr.dat present']
        ev.append(f'symmetry-reduced: yes, wannhr_symm outputs in {w["name"]}/{sym[-1]["folder"]}' if sym else
                  "symmetry-reduced: no wannier_sym outputs")
        bm = None
        if project:
            from .wannier_compare import band_match, band_match_text
            try:
                bm = band_match(project, w["name"], compute=compute_band_match)
            except Exception as exc:
                bm = {"status": "error", "reason": f"{type(exc).__name__}: {exc}"}
            ev.append(band_match_text(bm))
        else:
            ev.append(f'band match check: comparison available ({w["name"]}/wannier90_band.dat)' if "wannier90_band.dat" in names
                      else "band match check: not checked (no wannier90_band.dat)")
        from .spin_projection import spin_entry
        spin = spin_entry(w["name"], bm)
        ev.append(spin["line"])
        winfo = win_info(w)
        if winfo:
            ev.append(f"wannier90.win: {winfo}")
        steps.append(step("e", DONE if sym else PARTIAL, ev))
        steps[-1]["band_match"] = {"run": w["name"], "result": bm, "pending": bool(project) and bm is None}
        steps[-1]["spin"] = spin
    elif wan:
        steps.append(step("e", PARTIAL, [f'{r["name"]}: Wannier90 input, no wannier90_hr.dat' for r in wan]))
    else:
        steps.append(step("e", MISSING, ["no wannier90_hr.dat and no Wannier90 input in this branch"]))
    # (f) transport
    tr = [(r, f) for r in br for f in r["details"].get("transport_files", [])]
    if tr:
        r, f = tr[0]
        figs = [p for x in br for p in x["details"].get("sub_files", []) + [g["name"] for g in x["details"]["files"]]
                if p.lower().endswith((".png", ".pdf", ".svg")) and ("ahc" in p.lower() or "shc" in p.lower() or "berry" in p.lower()
                                                                     or "cond" in p.lower())]
        kind = "AHC" if "ahc" in f.lower() else "conductivity"
        ev = [f'{kind}: {r["name"]}/{f}' + (f" and {len(tr) - 1} more output files" if len(tr) > 1 else ""),
              f"plot: {r['name']}/{figs[0]}" if figs else "plot: not available"]
        steps.append(step("f", DONE, ev))
    else:
        steps.append(step("f", MISSING, ["no AHC or other conductivity output in this branch"]))
    nxt = next_calculation(branch, steps, br, soc)
    return {"branch": branch["name"], "steps": steps, "next": nxt}


# Key information of wannier90.win through the parsers module (energy windows, projections, num_wann), or None.
def win_info(run):
    try:
        from .parsers_extra import parse_wannier_win
    except ImportError:
        return None
    path = Path(run["path"]) / "wannier90.win"
    if not run["details"].get("has_win") or not path.exists():
        return None
    try:
        res = parse_wannier_win(path)
    except Exception:
        return None
    return None if not res or res.get("error") else res.get("summary")


# One sentence naming the next calculation: the first missing step (partial counts for c, d, e), then a partial one.
def next_calculation(branch, steps, br, soc):
    pick = next((s for s in steps if s["state"] == MISSING or (s["state"] == PARTIAL and s["step"] in "cde")), None) or \
        next((s for s in steps if s["state"] == PARTIAL), None)
    if pick is None:
        return "Every step of the loop has outputs in this branch; the next step is the physicist's check, not a new run."
    k = pick["step"]
    ref = next((r for r in br if r["details"].get("structure_file")), None)
    refn = ref["name"] if ref else "the first run"
    scfs = [r for r in br if r["details"].get("phase") == "scf"]
    conv = [r for r in scfs if (r["details"].get("scf_converged") or {}).get("converged") and r["status"] == "completed"]
    ncl = soc or branch["realization"] in ("coplanar", "non-coplanar") or any(noncollinear(r["details"]) for r in br)
    nc = "ISPIN = 2, LNONCOLLINEAR = .TRUE., " if soc else ("LNONCOLLINEAR = .TRUE., " if ncl else "ISPIN = 2, ")
    if k == "a":
        return f"Next: write POSCAR and INCAR MAGMOM for {branch['name']} so the structure and moments are known."
    if k == "b":
        return f"Next: {pick['action'] or 'define MSG'}; this needs no INCAR change."
    if k == "c":
        if scfs:
            last = scfs[-1]
            nelm = last["details"].get("incar_raw", {}).get("NELM", "the default 60")
            return (f"Next: restart the SCF {last['name']} from its WAVECAR with ISTART = 1, ICHARG = 1, a larger NELM (now {nelm}) "
                    f"and EDIFF unchanged, until OSZICAR's last step is below EDIFF.")
        return (f"Next: run a SCF for {branch['name']} with the POSCAR and MAGMOM of {refn}: ICHARG = 2, ISTART = 0, NSW = 0, "
                f"{nc}LSORBIT = {tf(soc)}, EDIFF = 1E-6, LCHARG = .TRUE., so that a band run can read its CHGCAR.")
    if k == "d":
        src = conv[-1]["name"] if conv else (scfs[-1]["name"] if scfs else None)
        frm = f"from the CHGCAR of {src}" if src else "from a converged SCF of this branch"
        bands = [r for r in br if r["details"].get("phase") == "band"]
        kp = f"the line-mode KPOINTS of {bands[-1]['name']}" if bands else "a labelled line-mode KPOINTS"
        return (f"Next: run the band structure {frm} with ICHARG = 11, {nc}LSORBIT = {tf(soc)}, LORBIT = 11 for projections, "
                f"the same POSCAR and MAGMOM, along {kp}, writing vasprun.xml.")
    if k == "e":
        src = conv[-1]["name"] if conv else refn
        if pick["state"] == PARTIAL and any(r["details"].get("has_hr") for r in br):
            w = [r for r in br if r["details"].get("has_hr")][-1]["name"]
            return f"Next: symmetrise {w}/wannier90_hr.dat with wannhr_symm (wannier_sym folder); this needs no INCAR change."
        return (f"Next: run the Wannier step from {src} with LWANNIER90 = .TRUE., ISTART = 1, ICHARG = 11, NUM_WANN set to the "
                f"number of projected orbitals and LSORBIT = {tf(soc)}, with projections and energy windows in wannier90.win.")
    hr = [r for r in br if r["details"].get("has_hr")]
    w = hr[-1]["name"] if hr else "the Wannier run"
    return (f"Next: compute the anomalous Hall conductivity from {w}/wannier90_hr.dat with postw90 (berry = true, berry_task = ahc "
            f"in wannier90.win); this needs no INCAR change.")


# Suggestions as plain text lines for the command line.
def suggestions_text(project, branch, sugg):
    lines = [f'Suggestions for the next DFT calculations: {project}, branch {branch["name"]} '
             f'({branch["realization"]}, {"SOC" if branch["soc"] else "non-SOC"}, {len(branch["runs"])} runs)']
    for s in sugg["steps"]:
        lines.append(f'  ({s["step"]}) {s["title"]}: {s["state"]}')
        lines += [f"      - {e}" for e in s["evidence"]]
        if s.get("action"):
            lines.append(f"      action: {s['action']}")
    lines.append(f"  {sugg['next']}")
    return "\n".join(lines)
