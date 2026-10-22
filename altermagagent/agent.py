"""Project chat: answers only from the indexed summary and runs, through the Anthropic SDK or a local Ollama server, with a band figure tool and a literature search tool."""

import json
import os
import urllib.error
import urllib.request
from datetime import datetime

from . import db
from .summary import convergence_line

MODEL = "claude-sonnet-5-5"
MAX_TOKENS = 4000
NO_KEY_MESSAGE = ("No API key set. Export ANTHROPIC_API_KEY to enable the agent, "
                  "or run Ollama locally (see README).")
OLLAMA_DEFAULT_HOST = "http://localhost:11434"
OLLAMA_DEFAULT_MODEL = "llama3.1:8b"
OLLAMA_MAX_TOOL_ROUNDS = 3
OLLAMA_NUM_CTX = 8192
OLLAMA_TIMEOUT = 600
NO_TOOLS_NOTE = ("Note: this local model cannot call tools, so it cannot draw figures; "
                 "use the Redraw button on a band run to draw the band figure.")
LOG_PATH = db.HOME_DIR / "agent_calls.log"
MAX_TOOL_ROUNDS = 6
EXCERPT_MAX_BYTES = 2 * 1024 * 1024
EXCERPT_MAX_LINES = 60
TOOLS = [{
    "name": "plot_bands",
    "description": ("Draw the band structure of an existing band run of this project from its vasprun.xml and KPOINTS, "
                    "along the high-symmetry k-path when KPOINTS is a labelled line-mode file. "
                    "Writes PNG, SVG and interactive HTML under ~/.altermagagent/figures, never into the run folder. "
                    "Returns the file paths, whether labels were found, spin polarisation, SOC and the largest up-down difference "
                    "among the sampled k-points and bands in the window (k-point index is 1-based)."),
    "input_schema": {"type": "object", "properties": {
        "run_name": {"type": "string", "description": "Name of the band run; omit to use the latest band run of the project."},
        "emin": {"type": "number", "description": "Lowest energy relative to the Fermi level, eV. Default -3."},
        "emax": {"type": "number", "description": "Highest energy relative to the Fermi level, eV. Default 3."},
        "style": {"type": "string", "enum": ["paper", "screen"],
                  "description": "paper = white, serif, thin lines, 300 dpi; screen = larger sans-serif figure."}},
        "required": []},
}, {
    "name": "search_literature",
    "description": ("Search published scientific papers through the Valency literature service. "
                    "Returns up to `limit` papers with title, first author, year, source, url and citation count, "
                    "or an empty list with a note 'literature unavailable: <reason>'. "
                    "These are papers found by search, not evidence from the project folder."),
    "input_schema": {"type": "object", "properties": {
        "query": {"type": "string", "description": "Search text, for example 'MnTe altermagnet spin splitting'."},
        "limit": {"type": "integer", "description": "Number of papers, 1 to 5. Default 3."}},
        "required": ["query"]},
}, {
    "name": "read_file_excerpt",
    "description": ("Read-only: return the first and last lines of a small text input or output file (for example INCAR, "
                    "KPOINTS, POSCAR, OUTCAR) inside an indexed run folder of this project. At most 60 lines in total. "
                    "Refuses files over 2 MB, binary files, CHGCAR, WAVECAR, PROCAR and paths outside the project folders."),
    "input_schema": {"type": "object", "properties": {
        "run_name": {"type": "string", "description": "Name of the run as in RUNS JSON, for example '2-SOC/2-bands'."},
        "file_name": {"type": "string", "description": "File name inside the run folder, for example 'KPOINTS' or 'INCAR'."},
        "max_lines": {"type": "integer", "description": "Total lines returned, 2 to 60. Default 60."}},
        "required": ["run_name", "file_name"]},
}]


WORKFLOW_KNOWLEDGE = (
    "Workflow knowledge (from Diana's procedure, docs/detection_rules_from_procedure.md). "
    "Each run has a phase: relax, scf, band, wannier, symmetrisation, projection, transport or other. "
    "Her order is SCF, then BAND, then WANNIER90, then symmetrisation with wannhr_symm, then transport (SHC or AHC). "
    "SCF: ICHARG=2, ISTART=0, NSW=0, LCHARG=.TRUE.; writes CHGCAR, CHG and WAVECAR; EDIFF=1E-6 is enough. "
    "If dE at the start did not reach 1E-5, NELMDL should be enlarged; ENCUT is the POTCAR maximum times 1.3. "
    "SOC or noncollinear runs need vasp_ncl, others vasp_std; the binary is read from the job script (run.wexac.sh or another *.sh with #BSUB or #SBATCH lines). "
    "BAND: ICHARG=11 with KLABELS or a line-mode KPOINTS; vaspkit option 211 writes BAND.dat, REFORMATTED_BAND.dat and KLABELS; "
    "band_plot.m plots them on the PC. "
    "WANNIER90: LWANNIER90=.TRUE., ISYM=-1, NUM_WANN in INCAR, wannier90.win with num_bands equal to NBANDS and write_hr=.true. "
    "Wannier90 must write wannier90.amn, .mmn and .eig must exist; job script slots should equal NBANDS and submit_wannier.sh slots NUM_WANN. "
    "Symmetrisation: wannhr_symm from WannierTools writes protmat.log, mysymmop.dat, magmom_primed.log, magmom_unprimed.log "
    "and wannier90_hr.dat_*nsymm* files. Its inputs are poscar.in, wann.in, locaxis.in and magmom.in; the basis must be in uudd order. "
    "Projection: PBAND_*, DOS_atom*, SELECTED_ATOMS_LIST (vaspkit 213 or 214, parse_vaspxml_v7.py, projband.py). "
    "Transport: the SHC/AHC code uses setting, wpos, 1-submit.sh to 3-submit.sh and writes ahc.dat or shc data and a result folder. "
    "Other: ICHARG=11 on a regular mesh is a Fermi surface step (pyprocar); wt.in is WannierTools surface states; "
    "trace.txt is vasp2trace; band.conf and POSCAR-001 are phonopy. "
    "Each run lists its phase, phase_reason, extra steps found in sub-folders, and local procedure checks marked ok, warn or info "
    "with an evidence state; warn means it differs from Diana's notes, pending her confirmation, not a physics error. "
    "Use these to explain what a run is and which step of her order is missing next.")


# Compact run list for the prompt, without the bulky structure dicts.
def run_digest(runs):
    out = []
    for r in runs:
        d = r["details"]
        out.append({"seq": r["seq"], "name": r["name"], "path": r["path"], "status": r["status"],
                    "status_note": d.get("status_note"), "kind": d.get("kind"), "phase": d.get("phase"),
                    "phase_reason": d.get("phase_reason"), "steps": d.get("steps"),
                    "local_procedure_checks_pending_diana": d.get("checks"),
                    "first_time": datetime.fromtimestamp(r["first_time"]).isoformat(timespec="minutes"),
                    "last_time": datetime.fromtimestamp(r["last_time"]).isoformat(timespec="minutes"),
                    "parent": r["parent"], "parent_note": r["parent_note"],
                    "parent_candidates": d.get("parent_candidates"), "incar": d.get("incar"),
                    "formula": d.get("formula"), "spacegroup": d.get("spacegroup"), "energy_eV": d.get("energy"),
                    "energy_source": d.get("energy_source"), "moments": d.get("moments"),
                    "convergence": convergence_line(r), "scf_converged": d.get("scf_converged"),
                    "files": [f["name"] for f in d.get("files", [])]})
    return out


# System prompt holding the project's summary JSON and run list.
def system_prompt(project, summary, runs):
    return (
        f"You are the AltermagAgent assistant for the VASP project '{project}'. "
        "Answer only from the project data below, which was read from the calculation folders. "
        "Cite the run names and file names that support each statement, for example 'scf_nosoc/OUTCAR'. "
        "If the data does not contain the answer, say so plainly and name the file or run that would be needed. "
        "Do not invent values. Keep answers short. "
        "Never state that a material is an altermagnet; state what the files support, using the branches in the summary "
        "and their suggestions (six steps marked done, partial, missing or not applicable, and the next calculation). "
        "A branch is one structure, one magnetic realization and one SOC setting. "
        "Literature classification is not read from these files. "
        "Give each fact its evidence state: output found, parsed, validated (with the stated tolerance) or insufficient evidence. "
        "For spin splitting say 'largest up-down difference among the sampled k-points and bands in the window', "
        "give the number of sampled k-points, and say 'none detected' rather than 'no splitting'. "
        "Procedure checks are local procedure checks from Diana's notes, pending her confirmation, not physics errors. "
        "A parent may be unknown or have several candidates; give the reason from parent_note. "
        "When asked for a band figure, call the plot_bands tool and report the PNG path and the notes it returns. "
        "When asked what is published or known about a material, you may call search_literature. "
        "Label every literature answer 'from the literature', cite each paper with its title and year, and keep it in a "
        "separate paragraph; never merge it with file evidence, and never let it change the suggestions from the files. "
        "If the tool returns 'literature unavailable', say so. "
        "You may call read_file_excerpt to read small input files such as INCAR, KPOINTS or POSCAR of an indexed run. "
        "For a question asking to prepare the next step or next calculation, the answer is a preview only: "
        "start with a line 'Preview only: nothing is written or submitted.', and never claim that a file was written or a job submitted. "
        "For such a question, pick the reference run (the similar finished run of the sibling branch) and the target run, "
        "read the KPOINTS and the INCAR of both with read_file_excerpt, and state whether the two KPOINTS paths match "
        "and which INCAR tags differ, citing the run and file names. "
        "Then give the order of steps and which files to copy from where.\n\n"
        f"{WORKFLOW_KNOWLEDGE}\n\n"
        f"SUMMARY JSON:\n{json.dumps(summary, ensure_ascii=False)}\n\n"
        f"RUNS JSON:\n{json.dumps(run_digest(runs), ensure_ascii=False)}")


# Append one call record with provider, model and token counts to the agent call log.
def log_call(project, input_tokens, output_tokens, note="", provider="anthropic", model=MODEL):
    db.HOME_DIR.mkdir(parents=True, exist_ok=True)
    rec = {"time": datetime.now().isoformat(timespec="seconds"), "project": project, "provider": provider,
           "model": model, "input_tokens": input_tokens, "output_tokens": output_tokens, "note": note}
    with open(LOG_PATH, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


# First and last lines of a small text file inside an indexed run folder of the project, as a JSON string.
def read_file_excerpt(project, run_name, file_name, max_lines=EXCERPT_MAX_LINES):
    from . import files
    runs = db.get_runs(project)
    want = (run_name or "").strip().rstrip("/")
    r = next((x for x in runs if x["name"].rstrip("/") == want), None)
    ends = [x for x in runs if want and x["name"].rstrip("/").endswith("/" + want)]
    if r is None and len(ends) == 1:
        r = ends[0]
    if r is None:
        return json.dumps({"error": f"no run named {want}; runs: {[x['name'] for x in runs]}"})
    roots = [os.path.realpath(os.path.expanduser(q["root"])) for q in db.list_projects()
             if q.get("name") == project and q.get("root")]
    run_dir = os.path.realpath(os.path.expanduser(r["path"]))
    target = os.path.realpath(os.path.join(run_dir, str(file_name or "")))
    inside = lambda path, root: path == root or path.startswith(root.rstrip(os.sep) + os.sep)
    if not inside(target, run_dir) or not any(inside(target, root) for root in roots):
        return json.dumps({"error": "refused: path is outside the indexed run folder"})
    if not files.may_open(target):
        return json.dumps({"error": f"refused: {os.path.basename(target)} is a large binary output, never opened"})
    if not os.path.isfile(target):
        return json.dumps({"error": f"no file {file_name} in run {r['name']}",
                           "files": [f["name"] for f in r["details"].get("files", [])]})
    size = os.path.getsize(target)
    if size > EXCERPT_MAX_BYTES:
        return json.dumps({"error": f"refused: file is {files.human_size(size)}, over the 2 MB limit"})
    try:
        n = max(2, min(EXCERPT_MAX_LINES, int(max_lines or EXCERPT_MAX_LINES)))
    except (TypeError, ValueError):
        n = EXCERPT_MAX_LINES
    info = files.excerpt(target, head=n - n // 2, tail=n // 2)
    if not info["opened"]:
        return json.dumps({"error": f"refused: {info.get('note', 'not opened')}"})
    return json.dumps({"run": r["name"], "file": files.display_path(target), "size": info["size_text"],
                       "modified": info["mtime_text"], "text": files.excerpt_text(info)}, ensure_ascii=False)


# Run one tool call from the model and return its result as a JSON string.
def run_tool(project, name, args):
    if name == "read_file_excerpt":
        return read_file_excerpt(project, args.get("run_name"), args.get("file_name"),
                                 args.get("max_lines") or EXCERPT_MAX_LINES)
    if name == "search_literature":
        from .literature import search_literature
        try:
            limit = max(1, min(5, int(args.get("limit") or 3)))
        except (TypeError, ValueError):
            limit = 3
        papers = search_literature(str(args.get("query") or ""), limit=limit)
        return json.dumps({"label": "from the literature, not from these files", "papers": list(papers),
                           "note": papers.note}, ensure_ascii=False)
    if name != "plot_bands":
        return json.dumps({"error": f"unknown tool {name}"})
    from .bands import plot_bands
    runs = db.get_runs(project)
    want = (args.get("run_name") or "").rstrip("/")
    if want:
        r = next((x for x in runs if x["name"].rstrip("/") == want), None)
    else:
        r = next((x for x in reversed(runs) if x["details"].get("kind") == "bands" and x["status"] == "completed"), None)
    if r is None:
        band_runs = [x["name"] for x in runs if x["details"].get("kind") == "bands"]
        return json.dumps({"error": f"no run named {want or '(band run)'}; band runs: {band_runs}"})
    try:
        res = plot_bands(r["path"], emin=args.get("emin", -3.0), emax=args.get("emax", 3.0),
                         style=args.get("style", "paper"), project=project, run_name=r["name"])
    except Exception as exc:
        return json.dumps({"error": f"{type(exc).__name__}: {exc}", "run": r["name"]})
    res.pop("splitting_per_band", None)
    res["run"] = r["name"]
    return json.dumps(res, ensure_ascii=False)


# Answer one question about a project; history is a list of {"role", "content"} dicts; the model may call plot_bands, search_literature or read_file_excerpt.
def chat(project, question, history=None):
    if not os.environ.get("ANTHROPIC_API_KEY"):
        if ollama_available():
            return chat_ollama(project, question, history)
        log_call(project, 0, 0, "no API key and no Ollama server, not sent", provider="none", model="")
        return NO_KEY_MESSAGE
    import anthropic
    p = db.get_project(project)
    if p is None:
        return f"Project {project} is not in the index."
    runs = db.get_runs(project)
    msgs = [m for m in (history or []) if m.get("role") in ("user", "assistant")]
    msgs.append({"role": "user", "content": question})
    system = system_prompt(project, p["summary"], runs)
    # A user-level key (sk-ant-usr) must name a workspace in a header; a workspace key needs nothing.
    headers = {}
    if os.environ.get("ANTHROPIC_WORKSPACE_ID"):
        headers["anthropic-workspace-id"] = os.environ["ANTHROPIC_WORKSPACE_ID"]
    client = anthropic.Anthropic(default_headers=headers) if headers else anthropic.Anthropic()
    for _ in range(MAX_TOOL_ROUNDS + 1):
        try:
            resp = client.messages.create(model=MODEL, max_tokens=MAX_TOKENS, system=system, messages=msgs, tools=TOOLS)
        except Exception as exc:
            log_call(project, 0, 0, f"error: {type(exc).__name__}")
            return f"The agent call failed: {type(exc).__name__}: {exc}"
        calls = [b for b in resp.content if getattr(b, "type", "") == "tool_use"]
        log_call(project, resp.usage.input_tokens, resp.usage.output_tokens, ",".join(b.name for b in calls))
        if resp.stop_reason != "tool_use" or not calls:
            return "".join(b.text for b in resp.content if getattr(b, "type", "") == "text").strip()
        msgs.append({"role": "assistant", "content": resp.content})
        msgs.append({"role": "user", "content": [{"type": "tool_result", "tool_use_id": b.id,
                                                  "content": run_tool(project, b.name, b.input)} for b in calls]})
    return "The agent stopped after too many tool calls."


# Base address of the Ollama server, from OLLAMA_HOST or the default.
def ollama_host():
    host = os.environ.get("OLLAMA_HOST", "").strip() or OLLAMA_DEFAULT_HOST
    if not host.startswith(("http://", "https://")):
        host = "http://" + host
    return host.rstrip("/")


# Model name for Ollama, from ALTERMAG_OLLAMA_MODEL or the default.
def ollama_model():
    return os.environ.get("ALTERMAG_OLLAMA_MODEL", "").strip() or OLLAMA_DEFAULT_MODEL


# True when an Ollama server answers GET /api/tags within 2 seconds.
def ollama_available():
    try:
        with urllib.request.urlopen(ollama_host() + "/api/tags", timeout=2) as resp:
            return resp.status == 200
    except Exception:
        return False


# The tool definitions in Ollama's tools format.
def ollama_tools():
    return [{"type": "function", "function": {"name": t["name"], "description": t["description"],
                                              "parameters": t["input_schema"]}} for t in TOOLS]


# POST one request to Ollama's /api/chat and return the decoded JSON reply; raises RuntimeError with the server message.
def ollama_post(payload):
    req = urllib.request.Request(ollama_host() + "/api/chat", data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"}, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=OLLAMA_TIMEOUT) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        body = exc.read().decode(errors="replace")
        try:
            body = json.loads(body).get("error", body)
        except ValueError:
            pass
        raise RuntimeError(f"HTTP {exc.code}: {body}") from None


# Answer one question through a local Ollama server, with up to three rounds of tool calls (plot_bands, search_literature).
def chat_ollama(project, question, history=None):
    p = db.get_project(project)
    if p is None:
        return f"Project {project} is not in the index."
    model = ollama_model()
    runs = db.get_runs(project)
    msgs = [{"role": "system", "content": system_prompt(project, p["summary"], runs)}]
    msgs += [{"role": m["role"], "content": m["content"]} for m in (history or [])
             if m.get("role") in ("user", "assistant") and isinstance(m.get("content"), str)]
    msgs.append({"role": "user", "content": question})
    tools = ollama_tools()
    for _ in range(OLLAMA_MAX_TOOL_ROUNDS + 1):
        payload = {"model": model, "messages": msgs, "stream": False, "options": {"num_ctx": OLLAMA_NUM_CTX}}
        if tools:
            payload["tools"] = tools
        try:
            resp = ollama_post(payload)
        except RuntimeError as exc:
            if tools and "does not support tools" in str(exc):
                log_call(project, 0, 0, "model does not support tools, plain chat", provider="ollama", model=model)
                tools = None
                continue
            log_call(project, 0, 0, f"error: {exc}", provider="ollama", model=model)
            return f"The local Ollama call failed: {exc}"
        except Exception as exc:
            log_call(project, 0, 0, f"error: {type(exc).__name__}", provider="ollama", model=model)
            return f"The local Ollama call failed: {type(exc).__name__}: {exc}"
        msg = resp.get("message") or {}
        calls = msg.get("tool_calls") or []
        names = [c.get("function", {}).get("name", "") for c in calls]
        log_call(project, resp.get("prompt_eval_count"), resp.get("eval_count"), ",".join(names),
                 provider="ollama", model=model)
        if not calls:
            text = (msg.get("content") or "").strip()
            return text if tools else f"{text}\n\n{NO_TOOLS_NOTE}".strip()
        msgs.append({"role": "assistant", "content": msg.get("content") or "", "tool_calls": calls})
        for c, name in zip(calls, names):
            args = c.get("function", {}).get("arguments") or {}
            if isinstance(args, str):
                try:
                    args = json.loads(args)
                except ValueError:
                    args = {}
            for key in ("emin", "emax"):
                try:
                    if key in args:
                        args[key] = float(args[key])
                except (TypeError, ValueError):
                    args.pop(key)
            msgs.append({"role": "tool", "tool_name": name, "content": run_tool(project, name, args)})
    return "The agent stopped after too many tool calls."
