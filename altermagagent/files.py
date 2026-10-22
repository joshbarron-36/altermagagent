"""Read-only file helpers: excerpts, folder listings and the pill links of each run."""

import os
from datetime import datetime
from pathlib import Path

NEVER_OPEN = {"CHGCAR", "WAVECAR", "PROCAR", "CHG", "WAVEDER", "AECCAR0", "AECCAR1", "AECCAR2", "LOCPOT", "ELFCAR"}
HEAD_LINES, TAIL_LINES = 8, 8
CAP_BYTES = 64 * 1024
LINE_CAP = 300


# Path as shown on screen: project name plus relative path inside an imported root, else the home folder as "~".
def display_path(path):
    if not path:
        return path
    text = str(path)
    tail = "/" if text.endswith("/") and len(text) > 1 else ""
    full = os.path.abspath(os.path.expanduser(text))
    best = None
    try:
        from . import db
        for proj in db.list_projects():
            root = os.path.abspath(os.path.expanduser(proj["root"]))
            if (full == root or full.startswith(root + os.sep)) and (best is None or len(root) > len(best[1])):
                best = (proj["name"], root)
    except Exception:
        best = None
    if best:
        rel = os.path.relpath(full, best[1])
        return best[0] + ("" if rel == "." else "/" + rel) + tail
    home = os.path.expanduser("~")
    if full == home or full.startswith(home + os.sep):
        return "~" + full[len(home):] + tail
    return text


# Human readable size such as 1.2 KB or 412 MB.
def human_size(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n} B"


# Modification time of a file as "YYYY-MM-DD HH:MM".
def fmt_mtime(t):
    return datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M")


# Return True when the file may be opened for an excerpt.
def may_open(path):
    return Path(path).name not in NEVER_OPEN


# Return metadata and the first and last lines of a text file; large binary outputs are never opened.
def excerpt(path, head=HEAD_LINES, tail=TAIL_LINES):
    p = Path(path)
    st = p.stat()
    info = {"name": p.name, "path": str(p), "size": st.st_size, "size_text": human_size(st.st_size),
            "mtime": st.st_mtime, "mtime_text": fmt_mtime(st.st_mtime), "head": [], "tail": [], "opened": False}
    if not may_open(p):
        info["note"] = f"Large file, not opened, reference only, {human_size(st.st_size)}"
        return info
    with open(p, "rb") as fh:
        chunk = fh.read(CAP_BYTES)
        if b"\x00" in chunk[:4096]:
            info["note"] = "Binary file, not shown"
            return info
        head_lines = chunk.decode("utf-8", errors="replace").splitlines()
        info["opened"] = True
        if st.st_size <= CAP_BYTES and len(head_lines) <= head + tail:
            info["head"] = [l[:LINE_CAP] for l in head_lines]
            return info
        info["head"] = [l[:LINE_CAP] for l in head_lines[:head]]
        fh.seek(max(0, st.st_size - CAP_BYTES))
        tail_lines = fh.read().decode("utf-8", errors="replace").splitlines()
        info["tail"] = [l[:LINE_CAP] for l in tail_lines[-tail:]]
    return info


# Excerpt as one display string with a separator between head and tail.
def excerpt_text(info):
    if info.get("note") and not info["opened"]:
        return info["note"]
    if not info["tail"]:
        return "\n".join(info["head"])
    return "\n".join(info["head"] + ["", "    ... (middle of the file not shown) ...", ""] + info["tail"])


# List a folder's files with sizes and dates, without opening any file.
def folder_listing(path):
    out = []
    with os.scandir(path) as it:
        for e in it:
            st = e.stat(follow_symlinks=False)
            out.append({"name": e.name + ("/" if e.is_dir() else ""), "size": st.st_size, "size_text": human_size(st.st_size),
                        "mtime": st.st_mtime, "mtime_text": fmt_mtime(st.st_mtime)})
    return sorted(out, key=lambda f: f["name"])


# Pick the main log file of a run: vasp.log first, then any *.log, *.out or *.err.
def log_file(names):
    if "vasp.log" in names:
        return "vasp.log"
    for n in sorted(names):
        if n.endswith((".log", ".out", ".err")):
            return n
    return None


# Pill link set of one run, by status first (error, running) and then by run kind.
def pills(run):
    d = run["details"]
    files = {f["name"]: f for f in d.get("files", [])}
    names = set(files)
    log = log_file(names)
    kind, status = d.get("kind"), run["status"]
    if status == "error":
        spec = [("Log", log, "file"), ("Output", "OUTCAR", "file")]
    elif status == "running":
        spec = [("Log tail", log, "tail"), ("Energy steps tail", "OSZICAR", "tail")]
    elif kind == "relax":
        spec = [("Final structure", "CONTCAR", "file"), ("Output", "OUTCAR", "file"), ("Energy steps", "OSZICAR", "file")]
    elif kind == "scf":
        spec = [("Structure file", "POSCAR", "file"), ("Output", "OUTCAR", "file"), ("Energy steps", "OSZICAR", "file"),
                ("Charge density", "CHGCAR", "ref")]
    elif kind == "bands":
        spec = [("Plot", "vasprun.xml", "plot"), ("Bands data", "vasprun.xml", "file"), ("Output", "OUTCAR", "file")]
    elif kind == "wannier":
        spec = [("Hamiltonian", "wannier90_hr.dat", "file"), ("Compare", "wannier90_hr.dat", "compare"),
                ("Input", "wannier90.win", "file"), ("Output", "wannier90.wout", "file")]
    else:
        spec = [("Output", "OUTCAR", "file")]
    out = []
    for label, fname, action in spec:
        if not fname or fname not in names:
            continue
        mono = fname
        if action == "ref":
            mono = f"{fname}, {human_size(files[fname]['size'])}"
        elif action == "plot":
            mono = "band figure"
        elif action == "compare":
            mono = "VASP versus Wannier bands"
        elif action == "tail":
            mono = f"{fname} (tail)"
        out.append({"label": label, "file": fname, "mono": mono, "action": action})
    sfile = next((f for f in ("POSCAR", "CONTCAR") if f in names), None)
    if sfile:
        out.append({"label": "Structure", "file": sfile, "mono": "moments figure", "action": "structure"})
    short = run["name"].rstrip("/") + "/"
    out.append({"label": "Folder", "file": "", "mono": short, "action": "folder"})
    return out
