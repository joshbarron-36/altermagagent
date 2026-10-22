"""Eight summary points, the branches with their suggestions and a gaps list built from the indexed runs of one project."""

import json
import re
import os
import warnings
from datetime import datetime
from pathlib import Path

STATUS_TEXT = {"completed": "completed", "error": "finished with error", "running": "possibly running", "missing": "missing"}
STATUS_DOT = {"completed": "ok", "error": "err", "running": "run", "missing": "miss"}
FERMI_WINDOW_EV = 1.5
FOUND, PARSED, VALID, INSUFF = "output found", "parsed", "validated", "insufficient evidence"
STATE_CLASS = {FOUND: "found", PARSED: "parsed", VALID: "valid", INSUFF: "insuff"}
NET_TOL_UB = 0.01           # declared tolerance on the net moment per cell, Bohr magneton
SPLIT_NONE_EV = 0.01        # up-down differences below this are reported as "none detected"
SPLIT_WORDS = "largest up-down difference among the sampled k-points and bands in the window"
PROCEDURE_LABEL = "local procedure check (Diana's notes), pending her confirmation"
LITERATURE_NOTE = "Literature classification is not read from these files."
PENDING_LABEL = "parsed, pending validation"
VALIDATIONS_DEFAULT = Path(__file__).resolve().parents[1] / "docs" / "validations.json"


# Physicist validations: ALTERMAG_VALIDATIONS (else docs/validations.json) plus ALTERMAG_HOME/validations.json if present.
def load_validations():
    home = Path(os.environ.get("ALTERMAG_HOME", Path.home() / ".altermagagent"))
    out = []
    for p in (Path(os.environ.get("ALTERMAG_VALIDATIONS", VALIDATIONS_DEFAULT)), home / "validations.json"):
        if p.is_file():
            try:
                out += json.loads(p.read_text())
            except (OSError, ValueError):
                warnings.warn(f"validations file not readable: {p}")
    return out


# The validation entry for one project, item and branch or run name, or None.
def find_validation(project, item, name, key):
    for v in load_validations():
        names = v.get(key) if isinstance(v.get(key), list) else [v.get(key)]
        if v.get("project") == project and v.get("item") == item and name in names:
            return v
    return None


# Comparable form of a label: no spaces, no "BNS", no brackets.
def norm_label(t):
    return "".join(str(t).replace("BNS", "").split()).replace("(", "").replace(")", "")


# Evidence tag after comparing a computed value with a validation: validated, differs, or the fallback tag.
def validation_tag(v, computed, fallback):
    if not v:
        return fallback
    who = v.get("who", "physicist").split()[0]
    try:
        when = datetime.strptime(v.get("date", ""), "%Y-%m-%d").strftime("%-d %b")
    except ValueError:
        when = v.get("date", "")
    if norm_label(computed) == norm_label(v["value"]):
        return f"validated by physicist ({who}, {when})"
    return f'{fallback}; differs from validated value {v["value"]} ({who}, {when})'


# Format a POSIX time as "YYYY-MM-DD HH:MM".
def fmt_time(t):
    return datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M") if t else "unknown"


# Format a moment with an explicit sign, using a real minus sign.
def signed(v):
    return f"+{v:.2f}" if v >= 0 else f"−{abs(v):.2f}"


# Build one summary point dict with its evidence state (output found, parsed, validated, insufficient evidence).
# bold lists the substrings of value, in order, that the interface shows in bold.
def point(key, label, value, status, evidence, state=INSUFF, state_note="", data=None, bold=None):
    return {"key": key, "label": label, "value": value, "status": status, "evidence": evidence,
            "state": state, "state_note": state_note or state, "data": data, "bold": bold or []}


# Split a stored procedure check into (level, text, state), accepting the older two-field form.
def check_parts(c):
    return (c[0], c[1], c[2]) if len(c) > 2 else (c[0], c[1], PARSED)


# Describe the site moments of one run as "Cr1 +3.12, Cr2 -3.12".
def site_text(mom, min_abs=0.05):
    parts = []
    for s in mom["sites"]:
        if abs(s["moment"]) >= min_abs:
            parts.append(f'{s.get("element", "?")}{s["site"]} {signed(s["moment"])}')
    return ", ".join(parts) or "all sites below 0.05 μB"


# Per-element moment text such as "Cr ±3.12 μB" when the sublattices mirror each other.
def element_text(mom, min_abs=0.5):
    by_el = {}
    for s in mom["sites"]:
        if abs(s["moment"]) >= min_abs:
            by_el.setdefault(s.get("element", "?"), []).append(s["moment"])
    parts = []
    for el, vals in by_el.items():
        mags = {round(abs(v), 2) for v in vals}
        if len(mags) == 1 and any(v > 0 for v in vals) and any(v < 0 for v in vals):
            parts.append(f"{el} ±{mags.pop():.2f} μB")
        else:
            parts.append(f"{el} " + ", ".join(signed(v) for v in vals) + " μB")
    return "; ".join(parts)


# Size of a 3-vector.
def norm(v):
    return sum(x * x for x in v) ** 0.5


# Largest site moment per magnetic element from one OUTCAR moment block, using the vector size when x, y, z are written.
def element_sizes(mom, min_abs=0.5):
    out = {}
    for st in mom["sites"]:
        size = norm(st["vector"]) if st.get("vector") else abs(st["moment"])
        if size >= min_abs:
            el = st.get("element", "?")
            out[el] = max(out.get(el, 0.0), size)
    return out


# One realization entry for the Magnetic point: class, run count, net and element moments with their source.
def realization_entry(cls, names, runs):
    rs = [r for r in runs if r["name"] in names]
    done = [r for r in rs if r["status"] == "completed" and r["details"].get("moments")]
    pref = [r for r in done if r["details"].get("phase") == "scf"] or done
    ent = {"class": cls, "n_runs": len(rs), "runs": [r["name"] for r in rs]}
    if not pref:
        from .branches import moment_vectors, pattern_text
        src = next((r for r in rs if r["details"].get("magmom_list")), rs[0] if rs else None)
        pat = pattern_text(moment_vectors(src)[0], src["details"].get("species")) if src else "moments unknown"
        ent.update(text=f"{cls} ({len(rs)} runs): INCAR MAGMOM only, {pat}, no OUTCAR moments", net=None, run=None,
                   evidence=f'{src["name"]}/INCAR MAGMOM: {pat}' if src else "no run")
        return ent
    src = pref[-1]
    mom = src["details"]["moments"]
    net = norm(mom["net_vector"]) if mom.get("net_vector") else abs(mom["net"])
    els = element_sizes(mom)
    el_txt = ", ".join(f"{el} {v:.2f} μB" for el, v in els.items()) or "no site above 0.5 μB"
    ent.update(text=f"{cls} ({len(rs)} runs): net {net:.2f} μB per cell, {el_txt} ({src['name']}/OUTCAR)", net=net,
               run=src["name"], elements=el_txt, compensated=net < NET_TOL_UB,
               evidence=f'{src["name"]}/OUTCAR magnetization ({", ".join(mom.get("components", [mom["component"]]))}): '
                        f'{site_text(mom)}, net {net:.3f} μB')
    return ent


# Magnetic point: every magnetic realization found, with its run count, net moment and element moments.
def magnetic_point(runs, branches=()):
    reals = {}
    for b in branches:
        if b.get("realization_index") is None:
            continue
        e = reals.setdefault(b["realization_index"], {"class": b["realization"], "runs": []})
        e["runs"] += b["runs"]
    if not reals:
        return point("magnetic", "Magnetic", "unknown, no structure with moments in any run", "miss", [], INSUFF,
                     "insufficient evidence: no POSCAR with MAGMOM or OUTCAR moments")
    ents = [realization_entry(e["class"], set(e["runs"]), runs) for _, e in sorted(reals.items())]
    val = "; ".join(e["text"] for e in ents)
    ev = [e["evidence"] for e in ents]
    parsed = [e for e in ents if e.get("net") is not None]
    tol = f"declared tolerance {NET_TOL_UB} μB per cell (net = size of the summed site moments)"
    if parsed and len(parsed) == len(ents) and all(e["compensated"] for e in parsed):
        st, note = VALID, f"parsed from OUTCAR; validated: every realization has a net moment below the {tol}"
    elif parsed:
        st, note = PARSED, ("parsed from OUTCAR and INCAR MAGMOM; " +
                            ", ".join(f'{e["class"]}: net {e["net"]:.3f} μB' for e in parsed) + f"; {tol}")
    else:
        st, note = PARSED, "parsed from INCAR MAGMOM only; no OUTCAR magnetization block in these runs"
    first = parsed[0] if parsed else ents[0]
    data = {"run": first.get("run"), "net": first.get("net"), "compensated": first.get("compensated"),
            "elements": first.get("elements"), "realizations": ents}
    return point("magnetic", "Magnetic", val, "ok", ev, st, note, data, bold=[e["class"] for e in ents])


# Largest |E_up - E_down| near E_F from one spin-polarised vasprun.xml, or None.
def splitting_from_vasprun(path):
    from pymatgen.electronic_structure.core import Spin
    from pymatgen.io.vasp import Vasprun
    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vr = Vasprun(str(path), parse_projected_eigen=False, parse_potcar_file=False)
    if not vr.is_spin or Spin.down not in vr.eigenvalues:
        return None
    up, dn = vr.eigenvalues[Spin.up][..., 0], vr.eigenvalues[Spin.down][..., 0]
    ef = vr.efermi if vr.efermi is not None else 0.0
    near = np.minimum(np.abs(up - ef), np.abs(dn - ef)) <= FERMI_WINDOW_EV
    if not near.any():
        return None
    diff = np.where(near, np.abs(up - dn), -1.0)
    k, b = np.unravel_index(int(diff.argmax()), diff.shape)
    kp = vr.actual_kpoints[k]
    return {"split": float(diff[k, b]), "kpoint": [round(float(x), 3) for x in kp], "k_index": int(k) + 1,
            "band": int(b) + 1, "e_rel": float((up[k, b] + dn[k, b]) / 2 - ef), "efermi": ef,
            "nk": int(up.shape[0]), "nb": int(up.shape[1])}


# Collinear non-SOC spin-polarised runs (ISPIN = 2, LSORBIT off, no LNONCOLLINEAR).
def collinear_nosoc_runs(runs):
    out = []
    for r in runs:
        d = r["details"]
        nc = str(d.get("incar_raw", {}).get("LNONCOLLINEAR", "")).strip(". ").upper() in ("T", "TRUE")
        if not d.get("soc") and not nc and d["incar"].get("ISPIN") == 2:
            out.append(r)
    return out


# Spin splitting scan from the non-SOC bands or SCF vasprun.xml; used by the benchmark, not shown in the summary.
def splitting_point(runs):
    lab = "Spin splitting scan"
    coll = collinear_nosoc_runs(runs)
    cands = [r for r in coll if r["details"].get("has_vasprun")]
    cands.sort(key=lambda r: (r["details"].get("kind") == "bands", r["seq"]), reverse=True)
    if not coll:
        why = "no collinear non-SOC run (ISPIN = 2, LSORBIT off) exists in the folder"
        return point("splitting", lab, f"insufficient evidence: {why}", "miss", [], INSUFF, f"insufficient evidence: {why}",
                     {"reason": why, "collinear_runs": []})
    if not cands:
        names = ", ".join(r["name"] for r in coll)
        why = f"the collinear non-SOC runs ({names}) have no vasprun.xml with eigenvalues"
        return point("splitting", lab, f"insufficient evidence: {why}", "miss", [], INSUFF, f"insufficient evidence: {why}",
                     {"reason": why, "collinear_runs": [r["name"] for r in coll]})
    tried = []
    for r in cands:
        try:
            res = splitting_from_vasprun(Path(r["path"]) / "vasprun.xml")
        except Exception as exc:
            tried.append(f'{r["name"]}/vasprun.xml: could not be parsed ({type(exc).__name__})')
            continue
        if not res:
            tried.append(f'{r["name"]}/vasprun.xml: no spin-polarised eigenvalues within ±{FERMI_WINDOW_EV} eV of E_F')
            continue
        line = r["details"].get("kpoints_line_mode")
        where = "sampled along the k-path" if line else "sampled on the k-mesh"
        k = ", ".join(f"{x:.3f}" for x in res["kpoint"])
        window = f"±{FERMI_WINDOW_EV} eV of E_F"
        res.update(run=r["name"], line_mode=line, sampled=where, detected=res["split"] >= SPLIT_NONE_EV)
        if res["detected"]:
            val = (f'{SPLIT_WORDS} (window {window}): {res["split"]:.2f} eV at k-point {res["k_index"]} ({k}), band {res["band"]}, '
                   f'{r["name"]}; {res["nk"]} k-points {where}')
        else:
            val = (f'none detected: {SPLIT_WORDS} (window {window}) is {res["split"]:.3f} eV, below {SPLIT_NONE_EV} eV; '
                   f'{res["nk"]} k-points {where}, {r["name"]}')
        ev = [f'{r["name"]}/vasprun.xml: {SPLIT_WORDS}: {res["split"]:.3f} eV at k-point {res["k_index"]} ({k}), '
              f'band {res["band"]}, {res["e_rel"]:+.2f} eV from E_F',
              f'{r["name"]}/vasprun.xml: E_F {res["efermi"]:.3f} eV, {res["nk"]} k-points {where}, {res["nb"]} bands per spin, '
              f'only states within {window} considered'] + tried
        note = (f"parsed from {r['name']}/vasprun.xml; not validated: only the {res['nk']} k-points {where} are read, "
                "no convergence check over the full Brillouin zone")
        return point("splitting", lab, val, "ok", ev, PARSED, note, res)
    why = "; ".join(tried)
    return point("splitting", lab, "insufficient evidence: no non-SOC eigenvalues could be read", "miss", tried, FOUND,
                 f"output found but not parsed: {why}", {"reason": why, "collinear_runs": [r["name"] for r in coll]})


# SOC point: "non-SOC: N runs, SOC: M runs".
def soc_point(runs):
    soc = [r["name"] for r in runs if r["details"].get("soc")]
    nosoc = [r["name"] for r in runs if not r["details"].get("soc")]
    ev = [f"{n}/INCAR: LSORBIT = .TRUE." for n in soc] + [f"{n}/INCAR: no LSORBIT" for n in nosoc]
    return point("soc", "SOC", f"non-SOC: {len(nosoc)} runs, SOC: {len(soc)} runs", "ok" if soc else "miss", ev,
                 PARSED if runs else INSUFF, "parsed from INCAR LSORBIT; no check applies")


# Modification time of one named file in a run, or None.
def file_time(run, name):
    for f in run["details"]["files"]:
        if f["name"] == name:
            return f["mtime"]
    return None


# Step entries of one kind (symmetrisation, transport, projection, other) across runs, as (run, entry) pairs.
def step_entries(runs, step):
    return [(r, e) for r in runs for e in r["details"].get("steps", {}).get(step, [])]


# Run-relative folder text of one step entry, such as "9.6-wannier/wannier_sym".
def step_folder(r, e):
    return r["name"].rstrip("/") + ("" if e["folder"] == "." else "/" + e["folder"])


# Text naming the wannhr_symm logs found, or None when there are none.
def symm_text(runs):
    sy = step_entries(runs, "symmetrisation")
    if not sy:
        return None, []
    r, e = sy[-1]
    logs = [f for f in e["files"] if f == "protmat.log" or (f.startswith("magmom_") and f.endswith(".log"))]
    ops = f', mysymmop.dat lists {e["symmops"]} operations' if e.get("symmops") else ", mysymmop.dat" if "mysymmop.dat" in e["files"] else ""
    ev = [f'{step_folder(x, y)}: {", ".join(y["files"])}' for x, y in sy]
    return f"symmetrised Hamiltonian present, wannhr_symm logs found ({step_folder(r, e)}: {', '.join(logs)}{ops})", ev


# Evidence note for the wannhr_symm logs: output found, parsed when mysymmop.dat lists operations, never validated.
def symm_note(runs):
    sy = step_entries(runs, "symmetrisation")
    if not sy:
        return None
    n = sy[-1][1].get("symmops")
    return (f"symmetry logs: output found, parsed: {n} operations listed; not validated against the magnetic state" if n
            else "symmetry logs: output found; not parsed (no operations counted); not validated against the magnetic state")


# Hamiltonian point: "generated: <branch names>", "prepared" or "none", one line; the files stay in the evidence.
def hamiltonian_point(runs, branches=()):
    from .branches import branch_of
    hr = [r for r in runs if r["details"].get("has_hr")]
    if hr:
        names = []
        for r in hr:
            b = branch_of(branches, r["name"])
            if b and b["name"] not in names:
                names.append(b["name"])
        ev = [f'{x["name"]}/wannier90_hr.dat, {fmt_time(file_time(x, "wannier90_hr.dat"))}' for x in hr]
        st, sev = symm_text(runs)
        note = "output found (wannier90_hr.dat); not parsed, its matrix elements are not read"
        sn = symm_note(runs)
        return point("hamiltonian", "Hamiltonian", "generated: " + (", ".join(names) or "no branch"), "ok",
                     ev + (sev or []), FOUND, note + (f"; {sn}" if sn else ""))
    prep = [r for r in runs if r["details"].get("has_win") or r["details"]["incar"].get("LWANNIER90") is True]
    if prep:
        return point("hamiltonian", "Hamiltonian", "prepared", "run",
                     [f'{x["name"]}: wannier90.win or LWANNIER90 present, no wannier90_hr.dat' for x in prep],
                     INSUFF, "insufficient evidence: Wannier90 input found, no wannier90_hr.dat")
    return point("hamiltonian", "Hamiltonian", "none", "miss", [], INSUFF,
                 "insufficient evidence: no Wannier90 files in any run")


# Group runs into distinct structures with StructureMatcher.
def structure_groups(runs):
    from pymatgen.core import Structure
    from .branches import matcher as strict_matcher
    matcher, groups = strict_matcher(), []
    for r in runs:
        sd = r["details"].get("structure")
        if not sd:
            continue
        s = Structure.from_dict(sd)
        for g in groups:
            if matcher.fit(g["structure"], s):
                g["runs"].append(r)
                break
        else:
            groups.append({"structure": s, "runs": [r]})
    return groups


# Structures point: distinct structures with formula and space group.
def structures_point(groups):
    if not groups:
        return point("structures", "Structures investigated", "none readable", "miss", [], INSUFF,
                     "insufficient evidence: no readable POSCAR or CONTCAR")
    descs, ev, bold = [], [], []
    for g in groups:
        d = g["runs"][0]["details"]
        sg = f'{d.get("spacegroup")} ({d.get("spacegroup_number")})' if d.get("spacegroup") else "space group unknown"
        bold.append(f'{d["formula"]} {sg}')
        descs.append(f'{d["formula"]} {sg}, {d.get("crystal_system") or "unknown"}')
        ev.append(f'{descs[-1]}, a = {d.get("a")} Å, c = {d.get("c")} Å: {", ".join(r["name"] + "/" + (r["details"].get("structure_file") or "?") for r in g["runs"])}')
    return point("structures", "Structures investigated", f"{len(groups)}, " + "; ".join(descs), "ok", ev, PARSED,
                 "parsed from POSCAR or CONTCAR with pymatgen and spglib (symprec 0.01); not validated against a reference",
                 bold=bold)


# Magnetic space group label without the BNS number, for example "P-62'm' (BNS 189.225)" -> "P-62'm'".
def msg_label(text):
    return re.sub(r"\s*\(BNS[^)]*\)", "", str(text or "")).strip()


# Magnetic point group of a magnetic space group symbol: lattice letter dropped, screws to rotations, glides to mirrors.
def magnetic_point_group(label):
    lab = msg_label(label)
    if not lab or not lab[0] in "PABCFIR":
        return None
    body = re.sub(r"_\d", "", lab[1:])
    return re.sub(r"[abcnde]", "m", body) or None


# MSG label, its source, agreement and validation for one branch from its step (b) result, or None without a label.
def branch_msg(b, sug, project):
    sb = next((x for x in sug["steps"] if x["step"] == "b"), None) if sug else None
    msg = (sb or {}).get("msg") or {}
    label = msg.get("findsym") or msg.get("computed")
    if not label:
        return None
    val = find_validation(project, "magnetic space group", b["name"], "branch")
    tag = validation_tag(val, label, PENDING_LABEL)
    if msg.get("findsym"):
        how = "from findsym.log" + (f', computed {"agrees" if msg.get("agreement") == "agree" else "differs"}'
                                    if msg.get("agreement") else "")
    else:
        how = "computed" + (", validated by " + val.get("who", "physicist").split()[0] if tag.startswith("validated") else "")
    return {"label": msg_label(label), "how": how, "tag": tag, "findsym": bool(msg.get("findsym")), "msg": msg}


# Symmetries point: the magnetic space group once per magnetic realization (it does not depend on SOC), with its point group.
def symmetry_point(groups, runs=(), branches=(), suggs=None, project=None):
    ev, per, seen, bold = [], [], {}, []
    for b in branches:
        idx = b.get("realization_index")
        if idx is None or b["realization"] == "non-magnetic":
            continue
        info = branch_msg(b, (suggs or {}).get(b["name"]), project)
        best = seen.get(idx)
        if best is None or (info and (best["info"] is None or (info["findsym"] and not best["info"]["findsym"]))):
            seen[idx] = {"class": b["realization"], "info": info, "branch": b["name"]}
        if info:
            ev.append(f'{b["name"]}: magnetic space group {info["label"]} ({info["how"]}; {info["tag"]})')
    classes = [e["class"] for e in seen.values()]
    for idx, e in seen.items():
        name = e["class"] if classes.count(e["class"]) == 1 else f'{e["class"]} {idx}'
        info = e["info"]
        bold.append(name)
        if info is None:
            per.append(f"{name}: define MSG")
            continue
        mpg = magnetic_point_group(info["label"])
        bold += [info["label"]] + ([mpg] if mpg else [])
        per.append(f'{name}: {info["label"]} ({info["how"]})' + (f", magnetic point group {mpg}" if mpg else ""))
    for r in runs:
        w = (r["details"].get("symmetry") or {}).get("wannhr") or {}
        if w.get("state") == PENDING_LABEL:
            ev.append(f'wannhr_symm ({r["name"]}/{w["folder"]}/mysymmop.dat): {w["n_operations"]} operations, '
                      f'{w["n_magnetic_unitary"]} kept without and {w["n_magnetic_antiunitary"]} kept with time reversal; '
                      f'{w.get("comparison", "not compared")} [{PENDING_LABEL}]')
    if not per:
        return point("symmetries", "Symmetries", "no magnetic realization", "miss", ev, INSUFF,
                     "insufficient evidence: no branch with magnetic moments")
    sn = symm_note(runs)
    tags = {e["info"]["tag"] for e in seen.values() if e["info"]}
    state = VALID if tags and all(t.startswith("validated") for t in tags) else PARSED
    note = ("magnetic space groups from findsym.log, else computed with spglib from the branch's first converged run; "
            "one per magnetic realization, since it does not depend on SOC; point group from the MSG symbol "
            "(screws to rotations, glides to mirrors); the crystal space group is in the Structures point") + (f"; {sn}" if sn else "")
    return point("symmetries", "Symmetries", "; ".join(per), "ok", ev, state, note, bold=bold)


# Transport point: postw90 or WannierBerri outputs present or not.
def transport_point(runs):
    tr = [(r["name"], f) for r in runs for f in r["details"].get("transport_files", [])]
    tasks = [f'{e["task"]} in {step_folder(r, e)}' for r, e in step_entries(runs, "transport") if e.get("task")]
    if tr:
        parsed = [(r["name"], f, t) for r in runs for f, t in r["details"].get("transport_parsed", {}).items()]
        if parsed:
            st = PARSED
            what = ", ".join(f"{n.rstrip('/')}/{f}: {t['rows']} rows of {t['cols']} numbers" for n, f, t in parsed)
            note = f"output found and parsed ({what}); not validated (no convergence or unit check available)"
        else:
            st, note = FOUND, "output found; not parsed (no numeric table read); not validated"
        kinds = sorted({"AHC" if "ahc" in f.lower() else "SHC" if "shc" in f.lower() else "conductivity" for _, f in tr})
        val = f"calculated: {', '.join(kinds)} ({len(tr)} output files)"
        ev = [f"{n}/{f}" for n, f in tr] + ([f"setting task: {', '.join(tasks)}"] if tasks else [])
        return point("transport", "Transport", val, "ok", ev, st, note)
    prep = step_entries(runs, "transport")
    if prep:
        return point("transport", "Transport", "prepared, no output files yet", "run",
                     [f'{step_folder(r, e)}: {", ".join(e["files"])}' for r, e in prep], INSUFF,
                     "insufficient evidence: transport inputs found, no output file")
    return point("transport", "Transport", "not calculated", "miss", ["no postw90 or WannierBerri output files in any run"],
                 INSUFF, "insufficient evidence: no transport output file in any run")


# True when a run is an scf whose electronic convergence was read and not reached.
def unconverged(d):
    return d.get("phase") == "scf" and (d.get("scf_converged") or {}).get("converged") is False


# One "Convergence: ..." line for a run row, with the source file and the evidence tag (physicist validation if any).
def convergence_line(r):
    d = r["details"]
    if d.get("phase") != "scf":
        return "Convergence: not applicable"
    cv = d.get("scf_converged") or {}
    ok, src = cv.get("converged"), cv.get("file")
    if ok is None:
        return f"Convergence: not read ({cv.get('text', 'no OSZICAR or OUTCAR')}) [{INSUFF}]"
    if cv.get("steps") is not None:
        nums = f"|dE| {cv['de']:.1E} {'below' if ok else 'above'} EDIFF {cv['ediff']:.0E} ({src})"
        body = f"reached in {cv['steps']} electronic steps, {nums}" if ok else \
            f"NOT reached, stopped at {cv['steps']} steps, {nums}"
    else:
        body = f"{'reached' if ok else 'NOT reached'} ({cv.get('text', src)})"
    val = find_validation(r.get("project"), "scf convergence", r.get("name"), "run")
    return f"Convergence: {body} [{validation_tag(val, 'reached' if ok else 'NOT reached', PARSED)}]"


# Status point: counts per status.
def status_point(runs):
    counts = {}
    for r in runs:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    order = ["completed", "error", "running", "missing"]
    val = ", ".join(f"{counts[s]} {STATUS_TEXT[s]}" for s in order if counts.get(s)) or "no runs"
    st = "run" if counts.get("running") else "err" if counts.get("error") else "miss" if counts.get("missing") else "ok"
    ev = [f'{r["name"]}: {STATUS_TEXT[r["status"]]}, {r["details"].get("status_note", "")}' for r in runs]
    unc = [r["name"] for r in runs if unconverged(r["details"])]
    if unc:
        val += f", {len(unc)} unconverged scf"
        ev.append(f'{len(unc)} unconverged scf: {", ".join(unc)} (OSZICAR last |dE| not below EDIFF)')
    return point("status", "Status", val, st, ev, PARSED if runs else INSUFF,
                 "parsed from completion and error markers in OUTCAR, vasprun.xml and logs")


# Last update point: newest file time across all runs.
def last_update_point(runs):
    if not runs:
        return point("last_update", "Last update", "unknown", None, [], INSUFF, "insufficient evidence: no runs")
    r = max(runs, key=lambda x: x["last_time"])
    newest = max(r["details"]["files"], key=lambda f: f["mtime"])
    return point("last_update", "Last update", f'{fmt_time(r["last_time"])} ({r["name"]})', None,
                 [f'{r["name"]}/{newest["name"]}, {fmt_time(newest["mtime"])}'], FOUND, "output found: file modification times")


# Gaps from simple rules over the runs.
def gaps(runs):
    out = []
    hr = [r for r in runs if r["details"].get("has_hr")]
    tr = any(r["details"].get("transport_files") for r in runs)
    if hr and not tr:
        out.append({"gap": "Transport not calculated although a Wannier Hamiltonian exists",
                    "evidence": [f'{r["name"]}/wannier90_hr.dat' for r in hr]})
    soc_bands = [r for r in runs if r["details"].get("soc") and r["details"].get("kind") == "bands"]
    if not any(r["status"] == "completed" for r in soc_bands):
        if soc_bands:
            ev = [f'{r["name"]}: {STATUS_TEXT[r["status"]]}' for r in soc_bands]
            out.append({"gap": "SOC bands missing: no SOC bands run has completed", "evidence": ev})
        else:
            out.append({"gap": "SOC bands missing: no SOC bands run (LSORBIT with ICHARG=11) found", "evidence": []})
    if not hr:
        out.append({"gap": "Hamiltonian not generated: no wannier90_hr.dat in any run", "evidence": []})
    elif not step_entries(runs, "symmetrisation"):
        out.append({"gap": "Hamiltonian not symmetrised: no wannhr_symm logs (protmat.log, mysymmop.dat, magmom_*.log) found",
                    "evidence": [f'{r["name"]}/wannier90_hr.dat' for r in hr]})
    return out


# Build the full summary dict (eight points, branches with suggestions, gaps) from index rows.
def build(runs, project=None, root=None, msg_logs=None):
    from .branches import build_branches, suggestions
    groups = structure_groups(runs)
    branches = build_branches(runs)
    suggs = {b["name"]: suggestions(b, runs, project=project, root=root, logs=msg_logs) for b in branches}
    points = [structures_point(groups), magnetic_point(runs, branches), symmetry_point(groups, runs, branches, suggs, project),
              soc_point(runs), hamiltonian_point(runs, branches), transport_point(runs), status_point(runs),
              last_update_point(runs)]
    for b in branches:
        b["suggestions"] = suggs[b["name"]]
    return {"points": points, "branches": branches, "gaps": gaps(runs), "msg_logs": msg_logs or [],
            "built_at": datetime.now().isoformat(timespec="seconds")}


# Parent text of a run: the parent with its note, "candidates: a, b", "parent unknown" or "none (first run)".
def parent_text(r):
    if r["parent"]:
        return f'{r["parent"]} ({r["parent_note"]})'
    note = r.get("parent_note") or ""
    if note.startswith("first run"):
        return f"none ({note})"
    if note.startswith(("candidates", "parent unknown")):
        return note
    return f"parent unknown ({note})" if note else "parent unknown"


# Short parent label for a card: the parent name, "unknown", "candidates a, b" or "none (first run)".
def parent_short(r):
    if r["parent"]:
        return r["parent"]
    note = r.get("parent_note") or ""
    if note.startswith("first run"):
        return "none (first run)"
    if note.startswith("candidates"):
        return note.split(" (")[0].split(";")[0].replace("candidates:", "candidates").strip()
    return "unknown"


# Reason sentence behind the parent choice, or "" when none is recorded.
def parent_reason(r):
    note = (r.get("parent_note") or "").strip()
    if note.startswith("parent unknown"):
        note = note[len("parent unknown"):].lstrip(" :(").rstrip(")")
    return note


# Render the summary and run list of one project as plain text lines; verbose adds the reason for each parent.
def as_text(project, summ, runs, literature=None, verbose=False):
    lines = [f"Storyline of {project}", "", "Summary"]
    for i, p in enumerate(summ["points"], 1):
        lines.append(f"  {i}. {p['label']} [{p.get('state', INSUFF)}]: {p['value']}")
        lines.append(f"       evidence: {p.get('state_note', '')}")
    if summ.get("branches"):
        lines += ["", "Branches (suggestions: altermag suggest <project> <branch>)"]
        for b in summ["branches"]:
            lines.append(f'  - {b["name"]}: {b["realization"]}, {"SOC" if b["soc"] else "non-SOC"}, {len(b["runs"])} runs')
            if b.get("suggestions"):
                lines.append(f'      {b["suggestions"]["next"]}')
    if literature:
        lines += [""] + literature
    lines += ["", "Gaps"] + ([f"  - {g['gap']}" for g in summ["gaps"]] or ["  none"])
    lines += ["", "Runs in historical sequence"]
    for r in runs:
        d = r["details"]
        lines.append(f'  {r["seq"]}. {r["name"]}  phase: {phase_text(d)}  {fmt_time(r["first_time"])[:10]}  '
                     f'{STATUS_TEXT[r["status"]]}  {params_line(d)}  parent: {parent_short(r)}')
        if verbose and parent_reason(r):
            lines.append(f"       why this parent: {parent_reason(r)}")
        lines.append(f"       {convergence_line(r)}")
        for c in d.get("checks", []):
            lvl, t, st = check_parts(c)
            if lvl == "warn":
                lines.append(f"       {PROCEDURE_LABEL} [{st}]: {t}")
    return "\n".join(lines)


# Phase of a run plus the extra steps found in it, such as "wannier + symmetrisation, transport".
def phase_text(d):
    ph = d.get("phase") or d.get("kind") or "other"
    extra = [k for k in ("symmetrisation", "projection", "transport") if k in d.get("steps", {}) and k != ph]
    return ph + (" + " + ", ".join(extra) if extra else "")


# One-line parameter digest for a run card.
def params_line(d):
    inc = d["incar"]
    keys = {"relax": ["ISIF", "ENCUT", "NSW"], "scf": ["LSORBIT", "SAXIS", "ENCUT"] if d.get("soc") else ["ISPIN", "MAGMOM", "ENCUT"],
            "bands": ["ICHARG", "ISPIN", "LSORBIT"], "wannier": ["LWANNIER90", "LSORBIT", "ISTART"]}.get(d.get("kind"), [])
    parts = []
    for k in keys:
        if k in inc:
            v = inc[k]
            v = ".TRUE." if v is True else ".FALSE." if v is False else v
            parts.append(f"{k}={v}" if k != "MAGMOM" else f"MAGMOM {v}")
    if d.get("error"):
        parts.append(f'error: {d["error"]["line"][:40]}')
    return " ".join(str(p) for p in parts)
