"""Read-only walk of a VASP project tree into the SQLite index."""

import fnmatch
import os
import re
import time
import warnings
from datetime import datetime
from pathlib import Path

from . import db

TAGS = ["ISPIN", "MAGMOM", "LSORBIT", "ICHARG", "ISTART", "ISIF", "ENCUT", "LWANNIER90", "NSW", "IBRION", "SAXIS", "NBANDS",
        "NUM_WANN", "EDIFF", "NELMDL", "ISYM", "LREAL", "NWRITE", "LORBIT", "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ"]
ERROR_MARKERS = ["ZBRENT", "ERROR", "Error EDDDAV", "VERY BAD NEWS", "BAD TERMINATION", "fatal error",
                 "Sub-Space-Matrix is not hermitian", "Segmentation fault", "DUE TO TIME LIMIT", "CANCELLED"]
ERROR_GLOBS = ["*.log", "*.out", "*.err", "OUTCAR"]
FIGURE_GLOBS = ["*.png", "*.pdf", "*.svg"]
TRANSPORT_GLOBS = ["*.wpout", "*-ahc*", "*_kubo*", "*boltzwann*", "*_elcond*", "*_seebeck*", "*-berry*",
                   "*wberri*", "*WannierBerri*", "*_tdf.dat", "*-morb*", "ahc*.dat", "shc*.dat", "*result/*"]
# Phase markers from Diana's procedure (docs/detection_rules_from_procedure.md).
BAND_MARKERS = {"KLABELS", "BAND.dat", "REFORMATTED_BAND.dat", "REFORMATED_BAND.dat", "KLINES.dat", "BAND_GAP"}
SYMM_GLOBS = ["protmat.log", "mysymmop.dat", "magmom_*.log", "wannier90_hr.dat_*nsymm*", "symmhr_addrptblock*.py",
              "new_rpt_block.dat"]
SYMM_INPUTS = {"poscar.in", "wann.in", "locaxis.in", "magmom.in"}
PROJ_GLOBS = ["PBAND_*", "DOS_atom*", "SELECTED_ATOMS_LIST", "parse_vaspxml.ini"]
TRANSPORT_INPUTS = {"setting", "1-submit.sh", "wpos"}
OTHER_MARKERS = {"wt.in": "WannierTools input (surface states)", "trace.txt": "vasp2trace output (band representations)",
                 "band.conf": "phonopy band input (phonon modes)", "band.yaml": "phonopy phonon bands",
                 "phonopy_disp.yaml": "phonopy displacements (phonon modes)", "findsym.in": "symmetr input (linear response symmetry)"}
TOOL_DEPTH = 3
NEVER_OPEN = {"CHGCAR", "WAVECAR", "PROCAR", "CHG", "WAVEDER", "AECCAR0", "AECCAR1", "AECCAR2", "LOCPOT", "ELFCAR"}
RUNNING_WINDOW_S = 3600
PARENT_LTOL = 0.01  # parent links: no volume rescaling, lattice lengths within 1 %, so a 2 % strain is another structure
SKIP_DIRS = {".git", ".venv", "__pycache__", "node_modules"}


# Parse an INCAR file into an upper-case tag dictionary of raw strings.
def parse_incar(path):
    tags = {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = re.split(r"[#!]", line, maxsplit=1)[0]
            for part in line.split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    if k.strip():
                        tags[k.strip().upper()] = v.strip()
    return tags


# Turn an INCAR value into int, float, bool or the raw string.
def tag_value(raw):
    if raw is None:
        return None
    s = raw.strip()
    low = s.lower().strip(".")
    if low in ("true", "t"):
        return True
    if low in ("false", "f"):
        return False
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


# Expand a MAGMOM string such as "2*3 -3" into a list of floats.
def expand_magmom(raw):
    out = []
    for tok in (raw or "").split():
        try:
            if "*" in tok:
                n, v = tok.split("*", 1)
                out += [float(v)] * int(float(n))
            else:
                out.append(float(tok))
        except ValueError:
            return []
    return out


# Yield lines of a text file without loading it whole.
def iter_lines(path):
    with open(path, "r", errors="replace") as fh:
        yield from fh


# Return the last bytes of a file decoded as text.
def tail_text(path, nbytes=4096):
    with open(path, "rb") as fh:
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - nbytes))
        return fh.read().decode("utf-8", errors="replace")


# List the files of a folder with size and modification time, without opening them.
def list_files(folder):
    out = []
    with os.scandir(folder) as it:
        for e in it:
            if e.is_file(follow_symlinks=False):
                st = e.stat(follow_symlinks=False)
                out.append({"name": e.name, "size": st.st_size, "mtime": st.st_mtime})
    return sorted(out, key=lambda f: f["name"])


# Read the structure from CONTCAR (when non-empty) or POSCAR and describe its symmetry.
def read_structure(folder, names):
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    for fname in ("CONTCAR", "POSCAR"):
        p = folder / fname
        if fname not in names or p.stat().st_size == 0:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = Structure.from_file(str(p))
        except Exception:
            continue
        info = {"structure_file": fname, "formula": s.composition.reduced_formula, "full_formula": s.formula,
                "species": [sp.symbol for sp in s.species], "structure": s.as_dict(),
                "a": round(s.lattice.a, 4), "c": round(s.lattice.c, 4)}
        try:
            sga = SpacegroupAnalyzer(s, symprec=0.01)
            info.update(spacegroup=sga.get_space_group_symbol(), spacegroup_number=sga.get_space_group_number(),
                        crystal_system=sga.get_crystal_system(), lattice_type=sga.get_lattice_type())
        except Exception:
            info.update(spacegroup=None, spacegroup_number=None, crystal_system=None, lattice_type=None)
        return info
    return {"structure_file": None}


# Find the first error marker line in log, out, err files or OUTCAR.
def find_error(folder, names):
    for name in names:
        if name in NEVER_OPEN or not any(fnmatch.fnmatch(name, g) for g in ERROR_GLOBS):
            continue
        for line in iter_lines(folder / name):
            for m in ERROR_MARKERS:
                if m in line:
                    return {"file": name, "line": line.strip()[:200]}
    return None


# Decide whether the run finished: OUTCAR timing, the vasprun.xml closing tag, or the Wannier90 "All done" line.
def find_completion(folder, names):
    if "OUTCAR" in names:
        for line in iter_lines(folder / "OUTCAR"):
            if "General timing and accounting" in line:
                return "OUTCAR"
    if "vasprun.xml" in names and tail_text(folder / "vasprun.xml").rstrip().endswith("</modeling>"):
        return "vasprun.xml"
    if "wannier90.wout" in names and "All done: wannier90 exiting" in tail_text(folder / "wannier90.wout"):
        return "wannier90.wout"
    return None


# Return the last total energy from OSZICAR F= or OUTCAR TOTEN, with its source file.
def read_energy(folder, names):
    if "OSZICAR" in names:
        last = None
        for line in iter_lines(folder / "OSZICAR"):
            m = re.search(r"F=\s*([-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?)", line)
            if m:
                last = float(m.group(1).replace("D", "E").replace("d", "e"))
        if last is not None:
            return last, "OSZICAR"
    if "OUTCAR" in names:
        last = None
        for line in iter_lines(folder / "OUTCAR"):
            m = re.search(r"free\s+energy\s+TOTEN\s*=\s*([-+]?\d+\.\d+)", line)
            if m:
                last = float(m.group(1))
        if last is not None:
            return last, "OUTCAR"
    return None, None


# Electronic convergence of the last ionic step: OSZICAR last |dE| below EDIFF, else OUTCAR "EDIFF is reached".
def read_convergence(folder, names, ediff):
    ediff = ediff if isinstance(ediff, (int, float)) else 1e-4
    if "OSZICAR" in names:
        last = None
        for line in iter_lines(folder / "OSZICAR"):
            m = re.match(r"\s*([A-Z]{2,4}):\s+(\d+)\s+(\S+)\s+(\S+)", line)
            if m:
                last = (m.group(1), int(m.group(2)), m.group(4))
        if last:
            try:
                de = abs(float(last[2].replace("D", "E")))
            except ValueError:
                de = None
            if de is not None:
                ok = de < ediff
                return {"converged": ok, "file": "OSZICAR", "steps": last[1], "de": de, "ediff": ediff,
                        "text": f"OSZICAR last electronic step {last[0]} {last[1]}: |dE| = {de:.1E} "
                                f"{'below' if ok else 'not below'} EDIFF = {ediff:.0E}"}
    if "OUTCAR" in names:
        for line in iter_lines(folder / "OUTCAR"):
            if "EDIFF is reached" in line or "reached required accuracy" in line:
                found = line.strip().strip("-").strip()
                return {"converged": True, "file": "OUTCAR", "text": f'OUTCAR: "{found}"'}
        return {"converged": False, "file": "OUTCAR", "text": "OUTCAR has no \"EDIFF is reached\" line"}
    return {"converged": None, "file": None, "text": "no OSZICAR or OUTCAR, convergence not read"}


# Parse the last magnetization (x/y/z) blocks of OUTCAR into per-site moments and the total.
def read_moments(folder, names):
    if "OUTCAR" not in names:
        return None
    blocks, comp, rows = {}, None, None
    for line in iter_lines(folder / "OUTCAR"):
        m = re.match(r"\s*magnetization \(([xyz])\)", line)
        if m:
            comp, rows = m.group(1), []
            continue
        if comp is None:
            continue
        s = line.split()
        if s and s[0].isdigit() and len(s) >= 2:
            try:
                rows.append((int(s[0]), float(s[-1])))
            except ValueError:
                pass
        elif s and s[0] == "tot":
            try:
                total = float(s[-1])
            except ValueError:
                total = sum(v for _, v in rows)
            blocks[comp] = {"sites": rows, "total": total}
            comp = None
    if not blocks:
        return None
    pick = "x"
    for c in ("z", "y", "x"):
        if c in blocks and max((abs(v) for _, v in blocks[c]["sites"]), default=0) > 0.05:
            pick = c
            break
    if pick not in blocks:
        pick = next(iter(blocks))
    b = blocks[pick]
    out = {"component": pick, "sites": [{"site": i, "moment": v} for i, v in b["sites"]], "net": b["total"],
           "components": sorted(blocks)}
    if len(blocks) > 1:
        comps = [c for c in ("x", "y", "z") if c in blocks]
        per = {c: dict(blocks[c]["sites"]) for c in comps}
        for s in out["sites"]:
            s["vector"] = [per[c].get(s["site"], 0.0) if c in per else 0.0 for c in ("x", "y", "z")]
        out["net_vector"] = [blocks[c]["total"] if c in blocks else 0.0 for c in ("x", "y", "z")]
    return out


# Return True when one file name matches any of the glob patterns.
def matches(name, globs):
    return any(fnmatch.fnmatch(name, g) for g in globs)


# Return the name of the main job script in a run folder: run.wexac.sh first, else any *.sh with #BSUB or #SBATCH lines.
def find_job_script(folder, names):
    if "run.wexac.sh" in names:
        return "run.wexac.sh"
    for name in sorted(names):
        if not name.endswith(".sh") or name == "submit_wannier.sh":
            continue
        try:
            path = folder / name
            if path.stat().st_size > 200_000:
                continue
            if any(re.match(r"\s*#(BSUB|SBATCH)\b", line) for line in iter_lines(path)):
                return name
        except OSError:
            continue
    return None


# Read the slot count ("#BSUB -n N" or "#SBATCH -n/--ntasks N") and the VASP binary (std, ncl, gam) from a job script, ignoring commented lines.
def read_job_script(path):
    out = {"slots": None, "vasp": None}
    try:
        if path.stat().st_size > 200_000:
            return out
        for line in iter_lines(path):
            m = re.match(r"\s*#BSUB\s+-n\s+(\d+)", line) or re.match(r"\s*#SBATCH\s+(?:-n\s*|--ntasks[=\s]\s*)(\d+)", line)
            if m:
                out["slots"] = int(m.group(1))
                continue
            if line.lstrip().startswith("#"):
                continue
            m = re.search(r"vasp[._-]?(\d+\.\d+\.\d+)?[^\s/]*/bin/vasp_(std|ncl|gam)|\bvasp_(std|ncl|gam)\b", line)
            if m:
                out["vasp"] = m.group(2) or m.group(3)
                ver = re.search(r"vasp[._-]?(\d+\.\d+\.\d+)", line)
                out["vasp_version"] = ver.group(1) if ver else None
    except OSError:
        pass
    return out


# Parse the key = value lines of a small wannier90.win into a lower-case dictionary.
def read_win(path):
    out = {}
    if path.stat().st_size > 1_000_000:
        return out
    for line in iter_lines(path):
        line = re.split(r"[#!]", line, maxsplit=1)[0]
        m = re.match(r"\s*([A-Za-z_]+)\s*[=:]\s*(\S+)", line)
        if m:
            out[m.group(1).lower()] = m.group(2)
    return out


# Return True when KPOINTS is a line-mode file (third line starts with L).
def kpoints_line_mode(folder, names):
    if "KPOINTS" not in names:
        return False
    try:
        lines = [ln for _, ln in zip(range(3), iter_lines(folder / "KPOINTS"))]
    except OSError:
        return False
    return len(lines) == 3 and lines[2].strip()[:1].upper() == "L"


# List files in sub-folders of a run that are not runs themselves, as relative paths, without opening them.
def sub_files(folder):
    out = []
    base = len(folder.parts)
    for dirpath, dirnames, filenames in os.walk(folder, followlinks=False):
        here = Path(dirpath)
        depth = len(here.parts) - base
        if depth and is_run_folder(set(filenames)):
            dirnames[:] = []
            continue
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")) if depth < TOOL_DEPTH else []
        if depth:
            rel = here.relative_to(folder).as_posix()
            out += [f"{rel}/{n}" for n in sorted(filenames)]
    return out


# Return True when a folder holds a VASP run (INCAR plus POSCAR or CONTCAR).
def is_run_folder(names):
    return "INCAR" in names and ("POSCAR" in names or "CONTCAR" in names)


# Find the tool steps (symmetrisation, projection, transport, other) present in a list of relative file paths.
def tool_steps(paths):
    steps = {}
    for p in paths:
        name = p.rsplit("/", 1)[-1]
        sub = p.rsplit("/", 1)[0] if "/" in p else "."
        if matches(name, SYMM_GLOBS) or name in SYMM_INPUTS:
            steps.setdefault("symmetrisation", {}).setdefault(sub, []).append(name)
        elif matches(name, PROJ_GLOBS):
            steps.setdefault("projection", {}).setdefault(sub, []).append(name)
        elif matches(name, TRANSPORT_GLOBS) or matches(p, TRANSPORT_GLOBS) or name in TRANSPORT_INPUTS:
            steps.setdefault("transport", {}).setdefault(sub, []).append(name)
        elif name in OTHER_MARKERS:
            steps.setdefault("other", {}).setdefault(sub, []).append(name)
    out = {}
    for step, by_dir in steps.items():
        for sub, names in by_dir.items():
            if step == "symmetrisation" and not any(matches(n, SYMM_GLOBS) for n in names):
                continue
            if step == "transport" and not (any(matches(n, TRANSPORT_GLOBS) for n in names) or TRANSPORT_INPUTS <= set(names) | {"wpos"}):
                continue
            out.setdefault(step, []).append({"folder": sub, "files": sorted(names)})
    return out


# Read the task line (ahc, shc, fullahc, fullshc) from an SHC/AHC setting file.
def read_setting_task(path):
    try:
        lines = [ln.strip() for ln in iter_lines(path)][:200]
    except OSError:
        return None
    for i, ln in enumerate(lines):
        if ln.lower().startswith("#task") and i + 1 < len(lines):
            return lines[i + 1]
    return None


# Count the symmetry operations listed in a wannhr_symm mysymmop.dat.
def count_symmops(path):
    n = 0
    for line in iter_lines(path):
        if re.match(r"\s*-{5,}\s+\d+\s+-{5,}", line):
            n += 1
    return n


# Count the numeric rows and columns of a small transport table such as ahc.dat; None when nothing numeric reads.
def read_numeric_table(path, max_bytes=5_000_000):
    try:
        if path.stat().st_size > max_bytes:
            return None
        rows, cols = 0, 0
        for line in iter_lines(path):
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            rows, cols = rows + 1, max(cols, len(vals))
        return {"rows": rows, "cols": cols} if rows else None
    except OSError:
        return None


# Evidence state of one procedure check text: output found, parsed or insufficient evidence.
def check_state(text):
    if " missing" in text or "not found" in text or "not set" in text:
        return "insufficient evidence"
    if " present" in text:
        return "output found"
    return "parsed"


# Primary phase of a run (relax, scf, band, wannier, other) from INCAR tags and files, with the reason.
def run_phase(incar, names, line_mode):
    if incar.get("LWANNIER90") is True or "wannier90.win" in names or "NUM_WANN" in incar:
        return "wannier", "LWANNIER90 = .TRUE. or wannier90.win or NUM_WANN"
    hits = sorted(names & BAND_MARKERS)
    if incar.get("ICHARG") == 11:
        if hits or line_mode:
            return "band", "ICHARG = 11 with " + (", ".join(hits) if hits else "line-mode KPOINTS")
        return "other", "ICHARG = 11 on a regular mesh (non-self-consistent, for example the Fermi surface step)"
    for m, what in OTHER_MARKERS.items():
        if m in names and m != "findsym.in":
            return "other", what
    if any(n.startswith("POSCAR-0") for n in names):
        return "other", "POSCAR-001 and following: phonopy displacement supercells"
    nsw = incar.get("NSW") or 0
    ibrion = incar.get("IBRION", 0 if nsw else -1)
    if isinstance(nsw, int) and nsw > 0 and ibrion not in (-1,):
        return "relax", f"NSW = {nsw}, IBRION = {ibrion}"
    return "scf", f'ICHARG = {incar.get("ICHARG", "default")}, NSW = {nsw}'


# Sanity checks from Diana's procedure for one run; each is (ok or warn or info, text).
def procedure_checks(folder, names, incar, phase, steps):
    out = []
    ediff = incar.get("EDIFF")
    if isinstance(ediff, (int, float)) and phase in ("scf", "band", "wannier"):
        out.append(("ok" if ediff <= 1e-6 else "warn", f"INCAR EDIFF = {ediff:.0E}; procedure: 1E-6 is enough for good convergence"))
    jname = find_job_script(folder, names)
    job = read_job_script(folder / jname) if jname else {}
    if job.get("vasp"):
        need = "ncl" if incar.get("LSORBIT") is True else None
        if need:
            out.append(("ok" if job["vasp"] == need else "warn", f"job script {jname} uses vasp_{job['vasp']}; SOC (LSORBIT) needs vasp_ncl"))
        else:
            out.append(("info", f"job script {jname} uses vasp_{job['vasp']}; procedure: ncl for SOC or noncollinear, otherwise std"))
    if phase == "band":
        for m, why in (("KLABELS", "k-point labels"), ("BAND.dat", "vaspkit option 211 output")):
            out.append(("ok" if m in names else "info", f"{m} {'present' if m in names else 'not found'} ({why})"))
    if phase == "wannier":
        for ext in ("amn", "mmn", "eig"):
            f = f"wannier90.{ext}"
            out.append(("ok" if f in names else "warn", f"{f} {'present' if f in names else 'missing'}"))
        nb, nw = incar.get("NBANDS"), incar.get("NUM_WANN")
        out.append(("ok" if nw else "warn", f"INCAR NUM_WANN = {nw}" if nw else "INCAR NUM_WANN not set; procedure asks to set it"))
        if "wannier90.win" in names:
            win = read_win(folder / "wannier90.win")
            if nb is not None and "num_bands" in win:
                same = str(nb) == win["num_bands"]
                out.append(("ok" if same else "warn", f"wannier90.win num_bands = {win['num_bands']}, INCAR NBANDS = {nb}"))
            if nw is not None and "num_wann" in win:
                same = str(nw) == win["num_wann"]
                out.append(("ok" if same else "warn", f"wannier90.win num_wann = {win['num_wann']}, INCAR NUM_WANN = {nw}"))
            wh = win.get("write_hr", "").lower().strip(".")
            out.append(("ok" if wh in ("true", "t") else "warn", f"wannier90.win write_hr = {win.get('write_hr', 'not set')}"))
        if job.get("slots") and nb:
            out.append(("ok" if job["slots"] == nb else "warn",
                        f"job script {jname} slots = {job['slots']}, INCAR NBANDS = {nb}; procedure: number of nodes equal to NBANDS"))
        if "submit_wannier.sh" in names and nw:
            sw = read_job_script(folder / "submit_wannier.sh")
            if sw.get("slots"):
                out.append(("ok" if sw["slots"] == nw else "warn",
                            f"job script submit_wannier.sh slots = {sw['slots']}, NUM_WANN = {nw}; procedure: nodes equal to NUM_WANN"))
        if steps.get("symmetrisation") and job.get("vasp_version"):
            ok = job["vasp_version"].startswith("5.4.4")
            out.append(("ok" if ok else "info", f"VASP {job['vasp_version']} in job script {jname}; the group's recommendation is 5.4.4 before wannhr_symm"))
    return [(lvl, text, check_state(text)) for lvl, text in out]


# Classify a run as relax, scf, bands or wannier from its INCAR and files.
def run_kind(incar, names):
    if incar.get("LWANNIER90") is True or "wannier90.win" in names:
        return "wannier"
    if incar.get("ICHARG") == 11:
        return "bands"
    nsw = incar.get("NSW") or 0
    ibrion = incar.get("IBRION", 0 if nsw else -1)
    if isinstance(nsw, int) and nsw > 0 and ibrion not in (-1,):
        return "relax"
    return "scf"


# Collect everything known about one run folder; never writes anything.
def scan_run(root, folder, now):
    files = list_files(folder)
    names = {f["name"] for f in files}
    raw = parse_incar(folder / "INCAR") if "INCAR" in names else {}
    incar = {k: tag_value(raw.get(k)) for k in TAGS if k in raw}
    if "MAGMOM" in raw:
        incar["MAGMOM"] = raw["MAGMOM"]
    times = [f["mtime"] for f in files]
    first, last = min(times), max(times)
    det = {"incar": incar, "incar_raw": raw, "files": files, "kind": run_kind(incar, names)}
    det.update(read_structure(folder, names))
    det["magmom_list"] = expand_magmom(raw.get("MAGMOM"))
    err = find_error(folder, names)
    done = find_completion(folder, names)
    if err:
        status, note = "error", f'{err["file"]}: "{err["line"]}"'
    elif done:
        status, note = "completed", f"completion marker in {done}"
    elif now - last < RUNNING_WINDOW_S:
        status, note = "running", f"files changed {datetime.fromtimestamp(last):%Y-%m-%d %H:%M}, no completion marker yet"
    else:
        status, note = "missing", "no completion marker, no error marker, no recent change"
    det.update(status_note=note, completion_file=done, error=err)
    det["energy"], det["energy_source"] = read_energy(folder, names)
    det["scf_converged"] = read_convergence(folder, names, incar.get("EDIFF"))
    mom = read_moments(folder, names)
    if mom:
        sp = det.get("species") or []
        for s in mom["sites"]:
            s["element"] = sp[s["site"] - 1] if 0 < s["site"] <= len(sp) else "?"
    det["moments"] = mom
    det["has_hr"] = "wannier90_hr.dat" in names
    det["has_win"] = "wannier90.win" in names
    det["has_wout"] = "wannier90.wout" in names
    det["has_vasprun"] = "vasprun.xml" in names
    det["figures"] = sorted(n for n in names if any(fnmatch.fnmatch(n, g) for g in FIGURE_GLOBS))
    det["transport_files"] = sorted(n for n in names if any(fnmatch.fnmatch(n, g) for g in TRANSPORT_GLOBS))
    det["soc"] = incar.get("LSORBIT") is True
    subs = sub_files(folder)
    steps = tool_steps(sorted(names) + subs)
    det["sub_files"] = subs
    det["transport_files"] += [p for p in subs if matches(p.rsplit("/", 1)[-1], TRANSPORT_GLOBS) or matches(p, TRANSPORT_GLOBS)]
    if "INCAR" in names:
        det["phase"], det["phase_reason"] = run_phase(incar, names, kpoints_line_mode(folder, names))
    else:
        main = next((p for p in ("symmetrisation", "transport", "projection", "other") if p in steps), "other")
        det["phase"], det["phase_reason"] = main, "no INCAR; tool outputs: " + ", ".join(steps.get(main, [{}])[0].get("files", []))[:200]
        if main in ("symmetrisation", "transport") and not det["error"]:
            status, det["status_note"] = "completed", f"{main} outputs present, no error marker"
    for st in steps.values():
        for e in st:
            base = folder if e["folder"] == "." else folder / e["folder"]
            if "mysymmop.dat" in e["files"]:
                e["symmops"] = count_symmops(base / "mysymmop.dat")
            if "setting" in e["files"]:
                e["task"] = read_setting_task(base / "setting")
    det["steps"] = steps
    try:
        from .symmetry import analyse
        det["symmetry"] = analyse(folder, "OUTCAR" in names, any("mysymmop.dat" in e["files"] for e in steps.get("symmetrisation", [])))
    except Exception as exc:
        det["symmetry"] = {"error": f"{type(exc).__name__}: {exc}"}
    det["kpoints_line_mode"] = kpoints_line_mode(folder, names)
    det["transport_parsed"] = {}
    for tf in det["transport_files"]:
        if tf.endswith(".dat"):
            tab = read_numeric_table(folder / tf)
            if tab:
                det["transport_parsed"][tf] = tab
    det["checks"] = procedure_checks(folder, names, incar, det["phase"], steps)
    rel = folder.relative_to(root).as_posix()
    return {"name": rel if rel != "." else folder.name, "path": str(folder), "status": status,
            "first_time": first, "last_time": last, "details": det}


# Return True when two runs hold the same structure according to StructureMatcher.
def same_structure(a, b, cache, matcher):
    key = (a["name"], b["name"])
    if key not in cache:
        sa, sb = a["details"].get("structure"), b["details"].get("structure")
        if not sa or not sb:
            cache[key] = False
        else:
            from pymatgen.core import Structure
            cache[key] = bool(matcher.fit(Structure.from_dict(sa), Structure.from_dict(sb)))
    return cache[key]


# Hubbard U settings of one run as a comparable tuple; runs without LDAU = .TRUE. give None.
def hubbard_key(incar):
    if incar.get("LDAU") is not True:
        return None
    nums = lambda k: tuple(round(float(x), 4) for x in str(incar.get(k, "")).split() if re.fullmatch(r"[-+]?\d*\.?\d+", x))
    return (incar.get("LDAUTYPE", 2), nums("LDAUL"), nums("LDAUU"), nums("LDAUJ"))


# Link every run to its parent when the files decide it; otherwise "parent unknown" or "candidates: a, b", with the reason.
def assign_parents(runs):
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from .branches import moment_vectors, same_pattern, classify
    matcher, cache = StructureMatcher(scale=False, ltol=PARENT_LTOL, angle_tol=1), {}
    for i, r in enumerate(runs):
        d = r["details"]
        r["parent"], d["parent_candidates"] = None, []
        if i == 0:
            r["parent_note"] = "first run in the folder"
            continue
        same = [e for e in runs[:i] if same_structure(r, e, cache, matcher)]
        if not same:
            r["parent_note"] = "parent unknown: no earlier run with the same structure"
            continue
        inc = d["incar"]
        not_nscf = [e for e in same if e["details"]["incar"].get("ICHARG") != 11]
        fnames = lambda e: {f["name"] for f in e["details"]["files"]}
        istart = inc.get("ISTART")
        if inc.get("ICHARG") == 11:
            tag, need = "ICHARG=11 reads CHGCAR", "CHGCAR"
        elif isinstance(istart, int) and istart > 0:
            tag, need = f"ISTART={istart} reads WAVECAR", "WAVECAR"
        else:
            tag, need = None, None
        socs = ""
        if need:
            # A restart file is only usable from a run with the same SOC setting and the same Hubbard U settings.
            hk = hubbard_key(inc)
            compat = [e for e in not_nscf if e["details"].get("soc") == d.get("soc") and hubbard_key(e["details"]["incar"]) == hk]
            pool = [e for e in compat if need in fnames(e)]
            socs = " and the same SOC setting" + (" and same U" if hk else "")
        else:
            mine = moment_vectors(r)[0]
            pool = [e for e in not_nscf if classify(moment_vectors(e)[0]) == classify(mine) and same_pattern(moment_vectors(e)[0], mine)]
            if not pool:
                r["parent_note"] = "parent unknown: no restart tag and no earlier same-structure run with the same magnetic realization"
                continue
            if len(pool) > 1:
                soc = [e for e in pool if e["details"].get("soc") == d.get("soc")]
                if soc:
                    pool, socs = soc, " and the same SOC setting"
        notag = "no restart tag (ICHARG=11 or ISTART>0)"
        if len(pool) == 1:
            r["parent"] = pool[0]["name"]
            r["parent_note"] = (f"{tag}; only earlier same-structure run with {need}{socs}; decided by the restart tag" if need else
                                f"{notag}; only earlier same-structure run{socs}; by structure only, unconfirmed")
        elif pool:
            d["parent_candidates"] = [e["name"] for e in pool]
            r["parent_note"] = (f"candidates: {', '.join(d['parent_candidates'])}; {tag or notag}, same structure{socs}; "
                                "nothing in the files decides between them")
        elif need and compat:
            # The restart file is gone now: link the only compatible earlier run (same SOC and U), preferring the same realization.
            mine = moment_vectors(r)[0]
            real = [e for e in compat if classify(moment_vectors(e)[0]) == classify(mine)
                    and same_pattern(moment_vectors(e)[0], mine)]
            pick = real or compat
            what = f"same structure{socs}{' and magnetic realization' if real else ''}"
            if len(pick) == 1:
                r["parent"] = pick[0]["name"]
                r["parent_note"] = f"{tag}; only compatible earlier run ({what}); restart file not present now"
            else:
                d["parent_candidates"] = [e["name"] for e in pick]
                r["parent_note"] = (f"candidates: {', '.join(d['parent_candidates'])}; {tag}, {what}; restart file not present now; "
                                    "nothing in the files decides between them")
        elif need:
            r["parent_note"] = f"parent unknown: {tag}, but no earlier same-structure run has the same SOC setting{' and same U' if hk else ''}"
        else:
            r["parent_note"] = "parent unknown: earlier same-structure runs are all non-self-consistent (ICHARG=11)"


# Walk a root folder read-only and return every VASP run found.
def find_runs(root, stats=None):
    root = Path(root).resolve()
    now = time.time()
    runs = []
    run_dirs = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS and not d.startswith("."))
        names, here = set(filenames), Path(dirpath)
        if stats is not None:
            stats["files"] = stats.get("files", 0) + len(filenames)
        inside = any(here.is_relative_to(d) for d in run_dirs)
        tool = not inside and tool_steps(sorted(names)) and any(k != "other" for k in tool_steps(sorted(names)))
        if is_run_folder(names) or tool:
            run_dirs.append(here)
            runs.append(scan_run(root, here, now))
        if stats is not None:
            for n in sorted(names):
                if n == "findsym.log" or fnmatch.fnmatch(n.lower(), "*msg*.log"):
                    rel = (here / n).relative_to(root).as_posix()
                    stats.setdefault("msg_logs", []).append(rel)
    runs.sort(key=lambda r: (r["first_time"], r["last_time"], r["name"]))
    for i, r in enumerate(runs):
        r["seq"] = i + 1
    assign_parents(runs)
    return runs


# Import one root as a project: scan, build the summary in memory, then store runs and summary in one write.
def import_project(root, name=None):
    from . import summary as summary_mod
    root_p = Path(root).expanduser().resolve()
    if not root_p.is_dir():
        raise FileNotFoundError(f"not a folder: {root_p}")
    if db.HOME_DIR.resolve().is_relative_to(root_p):
        raise ValueError("the index folder lies inside the scanned root; refusing to scan")
    name = name or root_p.name
    t0 = time.perf_counter()
    stats = {}
    runs = find_runs(root_p, stats)
    stamp = datetime.now().isoformat(timespec="seconds")
    details = {"n_runs": len(runs), "n_files": stats.get("files", 0), "scan_seconds": 0.0,
               "msg_logs": stats.get("msg_logs", [])}
    summ = summary_mod.build(db.with_parent_choices(name, runs), project=name, root=str(root_p),
                             msg_logs=details["msg_logs"])
    details["scan_seconds"] = round(time.perf_counter() - t0, 2)
    db.save_project(name, str(root_p), stamp, runs, summary=summ, details=details)
    return name, len(runs)


# Re-import an indexed project from its stored root under the same name; a failed rescan leaves the index as it was.
def rescan_project(name):
    p = db.get_project(name)
    if p is None:
        raise ValueError(f"no project named {name}")
    return import_project(p["root"], name)
