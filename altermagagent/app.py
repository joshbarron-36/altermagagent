"""NiceGUI three-panel interface: projects, storyline, agent chat; reads only from the index."""

import base64
import html
from datetime import datetime
from pathlib import Path

from nicegui import run as ng_run
from nicegui import app, ui

from . import agent, bands, branches as brmod, db, files, literature, structure
from .summary import (INSUFF, PARSED, SPLIT_NONE_EV, STATE_CLASS, STATUS_DOT, STATUS_TEXT, check_parts,
                      convergence_line, parent_reason, parent_short, parent_text, phase_text, unconverged)

BADGE = {"completed": "s-ok", "error": "s-err", "running": "s-run", "missing": "s-miss"}

CSS = """
<style>
:root { --bg:#000; --panel:#0f1115; --line:#262a33; --text:#e6e8ec; --muted:#9aa3ad; --soft:#181b21;
  --accent:#4f8cff; --accent-soft:#14213d; --ok:#34d399; --ok-bg:#0f2a1f; --err:#f87171; --err-bg:#3a1414;
  --run:#fbbf24; --run-bg:#3a2a08; --miss:#a1a1aa; --miss-bg:#26272b; }
html, body { margin:0; padding:0; background:var(--bg) !important; color:var(--text);
  font-family: system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif; font-size:16px; line-height:1.45; }
.nicegui-content { padding:0 !important; gap:0 !important; }
.q-page { background:var(--bg); }
button { cursor:pointer; font:inherit; color:inherit; }
.appbar { display:flex; align-items:center; gap:12px; padding:10px 16px; background:var(--panel);
  border-bottom:1px solid var(--line); width:100%; }
.appbar h1 { font-size:1.15rem; margin:0; font-weight:700; line-height:1.4; }
.appbar .dot { width:12px; height:12px; border-radius:50%; background:var(--accent); flex:none; }
.appbar .sub { color:var(--muted); font-size:0.95rem; }
.layout { display:grid; grid-template-columns:20% 55% 25%; width:100%; height:calc(100vh - 54px - 34px); min-height:480px; }
.panel { background:var(--panel); border-right:1px solid var(--line); overflow-y:auto; min-width:0; display:flex; flex-direction:column; }
.panel:last-child { border-right:0; }
.panel-head { padding:12px 16px; border-bottom:1px solid var(--line); font-weight:700; position:sticky; top:0;
  background:var(--panel); z-index:2; }
.panel-body { padding:12px 16px; flex:1; min-width:0; }
.proj-btn { width:100%; min-height:44px; display:flex; align-items:center; justify-content:space-between; gap:8px;
  padding:8px 12px; border:1px solid transparent; border-radius:8px; background:none; text-align:left; overflow-wrap:anywhere; }
.proj-btn:hover { background:var(--soft); }
.proj-btn.current { background:var(--accent-soft); border-color:var(--accent); font-weight:600; }
.proj-btn .meta { font-size:0.85rem; color:var(--muted); font-weight:400; flex:none; }
.proj-foot { padding:12px 16px; border-top:1px solid var(--line); }
.btn { min-height:44px; padding:8px 14px; border-radius:8px; border:1px solid var(--line); background:var(--panel); font-weight:600; }
.btn:hover { background:var(--soft); }
.btn.primary { background:var(--accent); border-color:var(--accent); color:#fff; }
.btn.primary:hover { background:#3b76e6; }
.link { background:none; border:0; padding:0 2px; color:var(--accent); text-decoration:underline; }
.badge { display:inline-block; font-size:0.85rem; font-weight:600; padding:2px 10px; border-radius:999px; white-space:nowrap; }
.s-ok { color:var(--ok); background:var(--ok-bg); } .s-err { color:var(--err); background:var(--err-bg); }
.s-run { color:var(--run); background:var(--run-bg); } .s-miss { color:var(--miss); background:var(--miss-bg); }
h2.sec { font-size:1.05rem; margin:8px 0; font-weight:700; line-height:1.4; }
.summary { margin:0 0 20px; padding:0; border:1px solid var(--line); border-radius:10px; overflow:hidden; width:100%; }
.sumrow { width:100%; min-height:40px; display:grid; grid-template-columns:minmax(110px, 30%) 1fr 20px; gap:10px;
  align-items:center; padding:7px 12px; background:none; border:0; text-align:left; }
.sumrow + .sumrow { border-top:1px solid var(--line); }
.sumrow:hover { background:var(--soft); }
.sumrow .lab { font-weight:600; } .sumrow .val { overflow-wrap:anywhere; min-width:0; }
.chev { color:var(--muted); font-size:1.2rem; text-align:right; }
.etag { display:inline-block; font-size:0.78rem; font-weight:600; padding:1px 8px; margin-right:6px; border-radius:999px;
  border:1px solid var(--line); white-space:nowrap; vertical-align:middle; }
.e-found { color:#93c5fd; border-color:#1e3a5f; } .e-parsed { color:var(--run); border-color:#5a4410; }
.e-valid { color:var(--ok); border-color:#145c3f; } .e-insuff { color:var(--miss); border-style:dashed; }
.diag { margin:0 0 20px; padding:10px 12px; border:1px solid var(--line); border-radius:10px; background:var(--soft); }
.diag h3 { font-size:1rem; margin:8px 0 2px; font-weight:700; } .diag h3:first-child { margin-top:0; }
.diag li { margin:2px 0 2px 18px; overflow-wrap:anywhere; }
.graph-wrap { overflow-x:auto; border:1px solid var(--line); border-radius:10px; background:var(--soft); padding:6px; max-width:100%; }
.graph-wrap svg { display:block; }
.legend { color:var(--muted); font-size:0.85rem; margin:6px 0 12px; }
.brhead { margin:14px 0 10px; padding:10px 12px; border:1px solid var(--accent); border-radius:10px; background:var(--accent-soft); }
.brhead h3 { margin:0 0 4px; font-size:1.1rem; } .brhead p { margin:2px 0; overflow-wrap:anywhere; }
.view3d { width:100%; height:440px; border:1px solid var(--line); border-radius:10px; background:#000; display:block; }
.sugg { margin:8px 0 20px; padding:10px 12px; border:1px solid var(--line); border-radius:10px; background:var(--soft); }
.sugg h4 { margin:10px 0 2px; font-size:1rem; } .sugg h4:first-child { margin-top:0; }
.sugg li { margin:2px 0 2px 18px; overflow-wrap:anywhere; }
.sugg .next { margin-top:10px; padding:8px 10px; border-left:3px solid var(--accent); background:var(--panel); overflow-wrap:anywhere; }
.st { display:inline-block; font-size:0.78rem; font-weight:700; padding:1px 8px; margin-left:6px; border-radius:999px; border:1px solid var(--line); }
.st-done { color:var(--ok); border-color:#145c3f; } .st-partial { color:var(--run); border-color:#5a4410; }
.st-missing { color:var(--err); border-color:#5c1d1d; } .st-na { color:var(--miss); border-style:dashed; }
details.fallback summary { cursor:pointer; color:var(--muted); min-height:32px; }
.lit a { color:var(--accent); text-decoration:underline; }
.lit h3 { font-size:1rem; line-height:1.4; margin:10px 0 4px; }
.psel { display:inline-flex; align-items:center; gap:4px; } .psel select { min-height:32px; max-width:100%;
  background:var(--panel); color:var(--text); border:1px solid var(--line); border-radius:6px; font:inherit; font-size:0.9rem; }
.caption { font-weight:600; margin:6px 0; }
.gaps { margin:8px 0 20px; padding:10px 12px; border:1px dashed var(--line); border-radius:10px; background:var(--soft); }
.gaps li { margin:2px 0 2px 18px; }
.timeline { margin:0; padding:0 0 0 22px; position:relative; width:100%; }
.timeline::before { content:""; position:absolute; left:8px; top:6px; bottom:6px; width:2px; background:var(--line); }
.tl-item { position:relative; margin-bottom:10px; }
.card, .files { line-height:1.35; }
.tl-item::before { content:""; position:absolute; left:-19px; top:18px; width:12px; height:12px; border-radius:50%;
  background:var(--panel); border:2px solid var(--muted); }
.tl-item.ok::before { border-color:var(--ok); background:var(--ok-bg); } .tl-item.err::before { border-color:var(--err); background:var(--err-bg); }
.tl-item.run::before { border-color:var(--run); background:var(--run-bg); }
.card { border:1px solid var(--line); border-radius:10px; padding:7px 12px; background:var(--panel); cursor:pointer;
  min-width:0; overflow-wrap:anywhere; }
.card:hover { border-color:var(--accent); box-shadow:0 0 0 2px var(--accent-soft); }
.card .top { display:flex; flex-wrap:wrap; align-items:center; gap:6px 10px; }
.mono, .card .name, .card .params { font-family: ui-monospace, Menlo, Consolas, monospace; }
.card .name { font-weight:700; } .card .date { color:var(--muted); }
.card .params { font-size:0.95rem; margin:2px 0; }
.card .parent { color:var(--muted); display:flex; align-items:center; flex-wrap:wrap; gap:4px; }
.card .fig { font-size:0.85rem; color:var(--accent); font-weight:600; }
.placeholder { border:1px dashed var(--line); border-radius:10px; padding:20px; color:var(--muted); background:var(--soft); }
.chat-head .ctx { font-weight:400; color:var(--muted); display:block; font-size:0.95rem; }
.msgs { display:flex; flex-direction:column; gap:10px; }
.msg { padding:10px 12px; border-radius:10px; max-width:95%; overflow-wrap:anywhere; white-space:pre-wrap; }
.msg .who { font-size:0.85rem; font-weight:700; color:var(--muted); display:block; margin-bottom:2px; }
.msg.mdmsg { white-space:normal; } .msg.mdmsg table { border-collapse:collapse; font-size:0.85rem; } .msg.mdmsg td, .msg.mdmsg th { border:1px solid var(--line); padding:3px 6px; vertical-align:top; overflow-wrap:normal; }
.msg.user { align-self:flex-end; background:var(--accent-soft); } .msg.agent { align-self:flex-start; background:var(--soft); }
.chat-form { display:flex; gap:8px; padding:12px 16px; border-top:1px solid var(--line); background:var(--panel); position:sticky; bottom:0; }
.files { display:flex; flex-wrap:wrap; gap:5px 6px; margin:5px 0 0; min-width:0; }
.fbtn > span { display:contents; }
.fbtn { min-height:28px; display:inline-flex; flex-wrap:wrap; align-items:center; gap:2px 6px; padding:2px 10px; border-radius:999px;
  border:1px solid var(--line); background:var(--soft); text-align:left; max-width:100%; overflow-wrap:anywhere; }
.fbtn:hover { border-color:var(--accent); }
.fbtn .fl { font-weight:600; font-size:0.88rem; }
.fbtn .fn { font-family: ui-monospace, Menlo, Consolas, monospace; font-size:0.8rem; color:var(--muted); min-width:0; }
.fbtn.ref { border-style:dashed; }
.dlg { background:var(--panel) !important; color:var(--text); width:100%; max-width:760px !important; min-width:0;
  box-shadow:0 10px 40px rgba(0,0,0,0.7); border-radius:12px; padding:0 !important; gap:0 !important; }
.dlg-head { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:10px 12px 10px 16px;
  border-bottom:1px solid var(--line); width:100%; }
.dlg-head h3 { margin:0; font-size:1.1rem; font-weight:700; line-height:1.4; overflow-wrap:anywhere; min-width:0; }
.close { min-width:44px; min-height:44px; border:1px solid var(--line); border-radius:8px; background:var(--panel); font-size:1.2rem; flex:none; }
.dlg-body { padding:12px 16px 16px; min-width:0; width:100%; }
.dlg-body h4 { margin:16px 0 6px; font-size:1rem; font-weight:700; line-height:1.4; }
.dlg-body p { margin:6px 0; }
.note { color:var(--muted); }
.tbl-wrap { overflow-x:auto; border:1px solid var(--line); border-radius:8px; max-width:100%; }
.tbl-wrap table { border-collapse:collapse; width:100%; min-width:420px; }
.tbl-wrap th, .tbl-wrap td { text-align:left; padding:8px 10px; border-bottom:1px solid var(--line); white-space:nowrap; }
.tbl-wrap tr:last-child td { border-bottom:0; }
.tbl-wrap th { background:var(--soft); font-weight:600; }
ul.lines { margin:4px 0; padding-left:20px; list-style:disc; }
pre.excerpt { margin:8px 0; padding:10px; background:#181b21; border:1px solid var(--line); border-radius:8px;
  overflow:auto; font-size:0.9rem; max-width:100%; max-height:18em; white-space:pre; }
.figure { border:1px solid var(--line); border-radius:10px; padding:10px; background:#fff; }
.figure img { display:block; width:100%; height:auto; }
.path { overflow-wrap:anywhere; }
.footer { text-align:center; color:var(--muted); font-size:0.85rem; padding:8px; height:34px; width:100%; }
@media (max-width: 799px) {
  .layout { display:block; height:auto; } .panel { border-right:0; min-height:auto; } .appbar .sub { display:none; }
  .sumrow { grid-template-columns:1fr 20px; } .sumrow .val { grid-column:1; } .sumrow .chev { grid-column:2; grid-row:1 / span 2; }
}
</style>
"""

esc = html.escape


# Short date of a POSIX time.
def day(t):
    return datetime.fromtimestamp(t).strftime("%Y-%m-%d") if t else "unknown"


# Badge markup for a run status.
def badge_html(status):
    return f'<span class="badge {BADGE[status]}">{esc(STATUS_TEXT[status])}</span>'


# Status badge of a run plus a red "unconverged" badge when its scf did not reach EDIFF.
def run_badges(r):
    extra = ' <span class="badge s-err">unconverged</span>' if unconverged(r["details"]) else ""
    return badge_html(r["status"]) + extra


# Short convergence text of an scf run for its card ("converged in 52 steps", "NOT converged, 200 steps"), else None.
def short_convergence(r):
    d = r["details"]
    if d.get("phase") != "scf":
        return None
    cv = d.get("scf_converged") or {}
    ok, n = cv.get("converged"), cv.get("steps")
    if ok is None:
        return "convergence not read"
    if n is None:
        return "converged" if ok else "NOT converged"
    return f"converged in {n} steps" if ok else f"NOT converged, {n} steps"


# Plain HTML element with trusted markup that we escaped ourselves.
def raw(markup, tag="div"):
    return ui.html(markup, sanitize=False, tag=tag)


# Button element with escaped text, a CSS class and a click handler.
def button(label, cls, handler, stop=False):
    return raw(esc(label), tag="button").classes(cls).on("click.stop" if stop else "click", handler)


# Escape a summary value and wrap its bold substrings (in order of appearance) in <strong>.
def bold_value(value, terms):
    out, pos = [], 0
    for t in terms or []:
        i = value.find(t, pos) if t else -1
        if i < 0:
            continue
        out += [esc(value[pos:i]), f"<strong>{esc(t)}</strong>"]
        pos = i + len(t)
    return "".join(out) + esc(value[pos:])


# Small evidence-state tag (output found, parsed, validated, insufficient evidence).
def etag(state):
    state = state or INSUFF
    return f'<span class="etag e-{STATE_CLASS.get(state, "insuff")}">{esc(state)}</span>'


# Two-column key/value table markup inside a horizontal scroll wrapper.
def kv_table(rows, h1, h2):
    body = "".join(f'<tr><td class="mono">{esc(str(a))}</td><td>{esc(str(b))}</td></tr>' for a, b in rows)
    return f'<div class="tbl-wrap"><table><thead><tr><th>{h1}</th><th>{h2}</th></tr></thead><tbody>{body}</tbody></table></div>'


class Page:
    """State and rendering of one browser tab."""

    # Set up per-tab state and draw the page.
    def __init__(self, project=None):
        projs = db.list_projects()
        names = [p["name"] for p in projs]
        self.project = project if project in names else (names[0] if names else None)
        self.chats = {}
        self.context = {}
        self.branch = {}
        self.build()

    # Current project's runs from the index.
    def runs(self):
        return db.get_runs(self.project) if self.project else []

    # Find one run of the current project by name.
    def find_run(self, name):
        return next((r for r in self.runs() if r["name"] == name), None)

    # Draw the static frame and the three refreshable panels.
    def build(self):
        ui.add_head_html(CSS)
        ui.on("altermag_branch", lambda e: self.select_branch(e.args))
        ui.query("body").style("background:#000")
        raw('<span class="dot"></span><h1>AltermagAgent</h1><span class="sub">Storylines of magnetic DFT projects</span>').classes("appbar")
        with ui.element("main").classes("layout"):
            with ui.element("section").classes("panel"):
                raw("Projects").classes("panel-head")
                with ui.element("div").classes("panel-body"):
                    self.projects_panel()
                with ui.element("div").classes("proj-foot"):
                    button("Import folder", "btn", self.open_import)
            with ui.element("section").classes("panel"):
                self.story_panel()
            with ui.element("section").classes("panel") as self.chat_section:
                self.chat_panel()
        raw(f"Read-only index at {esc(files.display_path(db.DB_PATH))}").classes("footer")

    # Left panel: one button per indexed project.
    @ui.refreshable_method
    def projects_panel(self):
        projs = db.list_projects()
        if not projs:
            raw('<p class="note">No projects yet. Use Import folder.</p>')
        for p in projs:
            cls = "proj-btn current" if p["name"] == self.project else "proj-btn"
            b = ui.element("button").classes(cls).on("click", lambda _, n=p["name"]: self.select(n))
            with b:
                raw(f'<span>{esc(p["name"])}</span>', tag="span")
                raw(f'<span class="meta">{p["n_runs"]} runs</span>', tag="span")

    # Switch to another project and redraw.
    def select(self, name):
        self.project = name
        self.refresh_all()

    # Redraw every panel.
    def refresh_all(self):
        self.projects_panel.refresh()
        self.story_panel.refresh()
        self.chat_panel.refresh()

    # Middle panel: summary rows, gaps and the run timeline.
    @ui.refreshable_method
    def story_panel(self):
        raw(f"Storyline of {esc(self.project or '')}" if self.project else "Storyline").classes("panel-head")
        with ui.element("div").classes("panel-body"):
            p = db.get_project(self.project) if self.project else None
            if not p:
                raw('<div class="placeholder"><p><strong>No project selected.</strong></p>'
                    '<p>Use Import folder to scan a directory. The scan is read-only.</p></div>')
                return
            runs = self.runs()
            raw('<h2 class="sec">Summary</h2>')
            with ui.element("div").classes("summary"):
                for i, pt in enumerate(p["summary"]["points"]):
                    b = ui.element("button").classes("sumrow").on("click", lambda _, x=pt: self.open_point(x))
                    with b:
                        raw(f'<span class="lab">{esc(pt["label"])}</span>', tag="span")
                        raw(f'<span class="val">{etag(pt.get("state"))}{bold_value(pt["value"], pt.get("bold"))}</span>', tag="span")
                        raw('<span class="chev">&#x203A;</span>', tag="span")
            brs = p["summary"].get("branches") or []
            if brs:
                sel = brmod.find_branch(brs, self.branch.get(self.project)) or brs[0]
                self.branch[self.project] = sel["name"]
                raw('<h2 class="sec">Branches</h2>')
                svg = brmod.graph_svg(brs, runs, selected=sel["name"], event="altermag_branch")
                raw(svg).classes("graph-wrap")
                raw(f'<p class="legend">Node colour = phase: {brmod.graph_legend()}. Lines = parent links; '
                    'click a node or a lane to open its branch.</p>')
                self.branch_section(sel, runs, p)
            else:
                raw('<h2 class="sec">Runs in historical sequence</h2>')
                with ui.element("div").classes("timeline"):
                    for r in runs:
                        self.run_card(r)
            if p["summary"]["gaps"]:
                items = "".join(f"<li>{esc(g['gap'])}</li>" for g in p["summary"]["gaps"])
                raw(f'<strong>Gaps</strong><ul class="lines">{items}</ul>').classes("gaps")

    # Select one branch from a graph click and redraw the storyline.
    def select_branch(self, args):
        name = args[0] if isinstance(args, list) and args else args
        if isinstance(name, dict):
            name = name.get("branch")
        if name:
            self.branch[self.project] = str(name)
            self.story_panel.refresh()

    # Branch section: header with the 3D view of its first run, its runs in order, then the suggestions.
    def branch_section(self, b, runs, p):
        by = {r["name"]: r for r in runs}
        mine = [by[n] for n in b["runs"] if n in by]
        soc = "SOC (LSORBIT = .TRUE.)" if b["soc"] else "non-SOC (no LSORBIT)"
        with ui.element("div").classes("brhead"):
            raw(f'<h3>Branch {esc(b["name"])}</h3><p>Magnetic realization: {esc(b["realization"])}, {esc(b["pattern"])} '
                f'(from {esc(b.get("moment_source", ""))}); {esc(soc)}; {len(mine)} runs.</p>')
            first = next((r for r in mine if r["details"].get("structure_file")), None)
            if first:
                self.view3d(first, f"Branch {b['name']}: structure of {first['name']}")
        raw(f'<h2 class="sec">Runs of {esc(b["name"])} in historical sequence</h2>')
        with ui.element("div").classes("timeline"):
            for r in mine:
                self.run_card(r)
        sugg = b.get("suggestions")
        if sugg:
            raw(f'<h2 class="sec">Suggestions for the next DFT calculations</h2>')
            cls = {brmod.DONE: "st-done", brmod.PARTIAL: "st-partial", brmod.MISSING: "st-missing", brmod.NA: "st-na"}
            parts = []
            self.refresh_band_match(sugg)
            for st in sugg["steps"]:
                items = "".join(f"<li>{etag(PARSED) if e.startswith('band match:') else ''}{esc(e)}</li>" for e in st["evidence"])
                act = f'<li><strong>Action:</strong> {esc(st["action"])}</li>' if st.get("action") else ""
                parts.append(f'<h4>({st["step"]}) {esc(st["title"])}<span class="st {cls.get(st["state"], "st-na")}">'
                             f'{esc(st["state"])}</span></h4><ul class="lines">{items}{act}</ul>')
            parts.append(f'<p class="next"><strong>Recommended next calculation.</strong> {esc(sugg["next"])}</p>')
            raw("".join(parts)).classes("sugg")
            self.band_match_block(sugg)
            sb = next((x for x in sugg["steps"] if x["step"] == "b"), None)
            if sb and sb.get("action"):
                self.define_msg_block(b, runs, p)
        self.literature_block(p["summary"], runs, b)

    # Replace the stored step (e) band-match line with the current cached result from the index.
    def refresh_band_match(self, sugg):
        from .wannier_compare import band_match, band_match_text
        for st in sugg["steps"]:
            bm = st.get("band_match")
            if not bm:
                continue
            try:
                res = band_match(self.project, bm["run"], compute=False)
            except Exception as exc:
                res = {"status": "error", "reason": f"{type(exc).__name__}: {exc}"}
            bm.update(result=res, pending=res is None)
            st["evidence"] = [band_match_text(res) if e.startswith(("band match", "preliminary: band match")) else e for e in st["evidence"]]

    # Step (e) band match: a Compare button for the figure, and the lazy first computation, cached in the index.
    def band_match_block(self, sugg):
        st = next((x for x in sugg["steps"] if x.get("band_match")), None)
        if not st:
            return
        bm = st["band_match"]
        with ui.element("div").classes("files"):
            with ui.element("button").classes("fbtn").on("click.stop", lambda _, n=bm["run"]: self.open_compare(n)):
                raw(f'<span class="fl">Compare</span><span class="fn">{esc(bm["run"])} band match figure</span>', tag="span")
        if bm["pending"]:
            from .wannier_compare import band_match
            project = self.project

            # Compute the band match once in a worker thread, then redraw the storyline with the cached result.
            async def compute():
                await ng_run.io_bound(band_match, project, bm["run"], True)
                self.story_panel.refresh()

            ui.timer(0.05, compute, once=True)

    # "define MSG" button of step (b): spglib on the branch's first converged run, findsym.log comparison, findsym.in written.
    def define_msg_block(self, b, runs, p):
        project = self.project
        area = ui.element("div").classes("sugg")

        # Compute the label off the event loop and list it, tagged "parsed, pending validation".
        async def run_define():
            area.clear()
            with area:
                raw('<p class="note">Computing the magnetic space group with spglib...</p>')
            msg = await ng_run.io_bound(brmod.define_msg, b, runs, project, p.get("root"), (p.get("details") or {}).get("msg_logs"))
            area.clear()
            with area:
                items = "".join(f"<li>{esc(e)}</li>" for e in brmod.msg_evidence(msg))
                raw(f'<h4>(b) define MSG for {esc(b["name"])}</h4><ul class="lines">{items}</ul>')

        with area:
            ui.button("define MSG", on_click=run_define).props("outline no-caps")
            raw('<p class="note">Runs spglib on the first converged run of this branch and writes findsym.in under the index '
                'folder; the label is parsed, pending validation.</p>')

    # Interactive 3Dmol.js view of a run's structure as the main view, with the PNG as a fallback; drawn on first use.
    def view3d(self, r, caption):
        project, name = self.project, r["name"]
        area = ui.element("div").classes("w-full")

        # Show the iframe with the 3D page and the PNG fallback under it.
        def show(png):
            area.clear()
            html_p = Path(png).with_suffix(".html")
            with area:
                raw(f'<p class="caption">{esc(caption)}: input moments from INCAR MAGMOM</p>')
                if html_p.exists():
                    url = "/figures/" + html_p.relative_to(bands.FIG_DIR).as_posix()
                    raw(f'<iframe class="view3d" src="{esc(url)}" title="3D structure of {esc(name)}" loading="lazy"></iframe>'
                        f'<p class="note"><a class="link" href="{esc(url)}" target="_blank" rel="noopener">Open the 3D view in a new tab</a>'
                        ' (3Dmol.js needs internet)</p>')
                data = base64.b64encode(Path(png).read_bytes()).decode()
                raw(f'<details class="fallback"{"" if html_p.exists() else " open"}><summary>Static PNG fallback</summary>'
                    f'<div class="figure"><img src="data:image/png;base64,{data}" alt="structure of {esc(name)}"></div></details>')

        # Draw the structure figure off the event loop, then show it.
        async def draw():
            try:
                res = await ng_run.io_bound(structure.plot_structure, r["path"], project=project, run_name=name)
            except Exception as exc:
                area.clear()
                with area:
                    raw(f'<p class="note">Structure view failed: {esc(f"{type(exc).__name__}: {exc}")}</p>')
                return
            show(res["paths"]["png"])

        png = structure.existing_figure(r["path"], project, name)
        if png:
            show(png)
        else:
            with area:
                raw('<p class="note">Drawing the structure from POSCAR and INCAR...</p>')
            ui.timer(0.05, draw, once=True)

    # Run dialog parameter line: INCAR changes versus the parent run (parsers_extra.incar_diff), else the five key tags.
    def card_params(self, r):
        d = r["details"]
        par = self.find_run(r["parent"]) if r.get("parent") else None
        try:
            from .parsers_extra import incar_diff
        except ImportError:
            incar_diff = None
        if incar_diff and par and (Path(r["path"]) / "INCAR").exists() and (Path(par["path"]) / "INCAR").exists():
            try:
                diff = incar_diff(Path(r["path"]) / "INCAR", Path(par["path"]) / "INCAR")
                if diff and not diff.get("error"):
                    return f"INCAR vs parent {par['name']}: {diff['summary']}"
            except Exception:
                pass
        inc = d["incar"]
        keys = ["LSORBIT", "ICHARG", "ISTART", "ENCUT", "EDIFF"]
        parts = []
        for k in keys:
            if k in inc:
                v = inc[k]
                parts.append(f"{k}={'.TRUE.' if v is True else '.FALSE.' if v is False else v}")
        line = " ".join(parts)
        if d.get("error"):
            line += f' error: {d["error"]["line"][:40]}'
        return line

    # Literature section under the suggestions of the selected branch: published papers found by search, filled off the event loop.
    def literature_block(self, summary, runs, branch=None):
        raw(f'<h2 class="sec">{esc(literature.HEADER)}</h2>')
        area = ui.element("div").classes("diag lit")
        with area:
            raw(f'<p class="note">{esc(literature.HEADER_NOTE)}</p><p class="note">Searching...</p>')

        # Run the search in a worker thread and show the papers with links, or why the literature is unavailable.
        async def fill():
            sec = await ng_run.io_bound(literature.section, summary, runs, 3, branch)
            parts = [f'<p class="note">{esc(sec["note"])}</p>']
            if sec["unavailable"]:
                parts.append(f'<p class="note">{esc(sec["unavailable"])}</p>')
            for q in sec["queries"]:
                items = "".join(f'<li><a href="{esc(x["url"])}" target="_blank" rel="noopener">{esc(x["title"])}</a> '
                                f'({esc(x["first_author"])}, {esc(x["year"])})</li>' for x in q["papers"]) or "<li>no papers found</li>"
                parts.append(f'<h3>Search "{esc(q["query"])}" (relates to: {esc(q["topic"])})</h3><ul class="lines">{items}</ul>')
            area.clear()
            with area:
                raw("".join(parts))

        ui.timer(0.05, fill, once=True)

    # One timeline entry: clickable card plus its row of pill links.
    def run_card(self, r):
        d = r["details"]
        with ui.element("div").classes(f'tl-item {STATUS_DOT[r["status"]]}'):
            with ui.element("div").classes("card").on("click", lambda _, n=r["name"]: self.open_run(n)):
                raw(f'<span class="name">{r["seq"]}. {esc(r["name"])}</span><span class="date">{day(r["first_time"])}</span>'
                    f'{run_badges(r)}<span class="fig">{esc(phase_text(d))}</span>').classes("top")
                conv = short_convergence(r)
                if conv:
                    raw(esc(conv)).classes("params")
                with ui.element("div").classes("parent"):
                    raw("Parent:", tag="span")
                    if r["parent"]:
                        button(r["parent"], "link", lambda _, n=r["parent"]: self.open_run(n), stop=True)
                    else:
                        raw(f"<span>{esc(parent_short(r))}</span>", tag="span")
                    if r["seq"] > 1:
                        self.parent_selector(r)
            self.pill_row(r)

    # Small "change parent" selector; the choice is written into the index under ~/.altermagagent only, never into data.
    def parent_selector(self, r):
        earlier = [x["name"] for x in self.runs() if x["seq"] < r["seq"]]
        opts = {"__auto__": "automatic", "__unknown__": "parent unknown"} | {n: n for n in earlier}
        cur = "__auto__" if not r.get("parent_user") else (r["parent"] or "__unknown__")
        project, run_name = self.project, r["name"]

        # Store the selected parent and redraw the storyline.
        def changed(e):
            v = e.value
            db.set_parent_choice(project, run_name, None if v == "__auto__" else "" if v == "__unknown__" else v)
            ui.notify(f"Parent of {run_name} stored in the index; no data folder was changed")
            self.story_panel.refresh()

        with ui.element("span").classes("psel").on("click.stop", lambda: None):
            raw("change parent", tag="span")
            ui.select(opts, value=cur, on_change=changed).props("dense dark outlined options-dense").classes("min-w-[9rem]") \
                .props(f'aria-label="change parent of {esc(run_name)}"')

    # Row of file and folder pills for a run.
    def pill_row(self, r):
        with ui.element("div").classes("files"):
            for pl in files.pills(r):
                cls = "fbtn ref" if pl["action"] == "ref" else "fbtn"
                with ui.element("button").classes(cls).on("click.stop", lambda _, x=pl, n=r["name"]: self.pill_action(n, x)):
                    raw(f'<span class="fl">{esc(pl["label"])}</span><span class="fn">{esc(pl["mono"])}</span>', tag="span")

    # Dispatch a pill click to the folder, plot or file dialog.
    def pill_action(self, run_name, pl):
        if pl["action"] == "folder":
            self.open_folder(run_name)
        elif pl["action"] == "plot":
            self.open_bands(run_name)
        elif pl["action"] == "structure":
            self.open_structure(run_name)
        elif pl["action"] == "compare":
            self.open_compare(run_name)
        else:
            self.open_file(run_name, pl["file"], pl["action"])

    # Create a dialog with the mockup's header and return (dialog, body container).
    def dialog(self, title):
        dlg = ui.dialog().props("transition-show=fade transition-hide=fade")
        dlg.on("hide", lambda: dlg.delete())
        with dlg, ui.card().classes("dlg"):
            with ui.element("div").classes("dlg-head"):
                raw(f"<h3>{esc(title)}</h3>")
                button("✕", "close", dlg.close)
            body = ui.element("div").classes("dlg-body")
        return dlg, body

    # Path line with a Copy path button.
    def path_block(self, path):
        raw(f'<p><strong>Path:</strong> <span class="mono path">{esc(files.display_path(path))}</span></p>')
        button("Copy path", "btn", lambda: (ui.clipboard.write(path), ui.notify("Path copied")))

    # Summary point dialog: value and evidence list.
    def open_point(self, pt):
        dlg, body = self.dialog(pt["label"])
        with body:
            raw(f"<p>{etag(pt.get('state'))}{esc(pt['value'])}</p><p class=\"note\">Evidence state: {esc(pt.get('state_note', ''))}</p>"
                "<h4>Evidence</h4>")
            items = "".join(f'<li class="mono">{esc(e)}</li>' for e in pt["evidence"]) or "<li>none recorded</li>"
            raw(f'<ul class="lines">{items}</ul>')
        dlg.open()

    # Run dialog: status, files, parameters, energy, moments and figure.
    def open_run(self, name):
        r = self.find_run(name)
        if not r:
            return
        self.context[self.project] = name
        self.chat_panel.refresh()
        d = r["details"]
        dlg, body = self.dialog(f"Run {name}")
        with body:
            par = "parent " + parent_short(r)
            raw(f'<p>{run_badges(r)} <span class="note">{day(r["first_time"])} to {day(r["last_time"])}, {esc(par)}</span></p>'
                f'<p class="mono path">{esc(files.display_path(r["path"]))}</p>')
            if parent_reason(r):
                raw(f'<p class="note"><strong>Why this parent:</strong> {esc(parent_reason(r))}</p>')
            if d.get("structure_file"):
                raw("<h4>Structure</h4>")
                self.view3d(r, f"Structure of {name}")
            raw("<h4>Files</h4>")
            self.pill_row(r)
            rows = [(k, ".TRUE." if v is True else ".FALSE." if v is False else v) for k, v in d["incar"].items()]
            raw("<h4>Parameters</h4>" + kv_table(rows or [("INCAR", "no tracked tags")], "Tag", "Value"))
            raw(f'<p class="note">{esc(self.card_params(r))}</p>')
            if d.get("phase") == "wannier" or d.get("has_win"):
                wi = brmod.win_info(r)
                if wi:
                    raw(f'<p class="note">wannier90.win: {esc(wi)}</p>')
            en = f'{d["energy"]:.4f} eV (from {d["energy_source"]})' if d.get("energy") is not None else "not available"
            st = f'{d.get("formula")}, {d.get("spacegroup")} ({d.get("spacegroup_number")}), from {d.get("structure_file")}' \
                if d.get("formula") else "not readable"
            raw("<h4>Energy and structure</h4>" + kv_table([("Total energy", en), ("Structure", st)], "Quantity", "Value"))
            mom = d.get("moments")
            if mom:
                mrows = [(f'{s["element"]}{s["site"]}', f'{s["moment"]:+.3f} μB') for s in mom["sites"]]
                mrows.append(("net", f'{mom["net"]:+.3f} μB'))
                raw(f'<h4>Moments, OUTCAR magnetization ({mom["component"]})</h4>' + kv_table(mrows, "Site", "Moment"))
            else:
                raw("<h4>Moments</h4><p>Not written in OUTCAR.</p>")
            raw(f'<h4>Status</h4><p>{esc(STATUS_TEXT[r["status"]])}. {esc(d.get("status_note", ""))}</p>')
            conv = convergence_line(r)
            body, _, tag = conv.rpartition(" [")
            raw(f'<p>{etag(tag.rstrip("]")) if body else ""}{esc(body or conv)}</p>')
            raw(f'<h4>Phase</h4><p>{esc(phase_text(d))}: {esc(d.get("phase_reason", ""))}</p>')
            if d.get("checks"):
                word = {"ok": "matches the notes", "warn": "differs from the notes", "info": "information"}
                rows = "".join(f'<tr><td>{esc(word.get(lvl, lvl))}</td><td>{etag(st)}</td><td>{esc(t)}</td></tr>'
                               for lvl, t, st in map(check_parts, d["checks"]))
                raw(f"<h4>Local procedure checks (Diana's notes), pending her confirmation</h4>"
                    '<p class="note">These compare the run with a local procedure, not with physics; they are not physics errors.</p>'
                    '<div class="tbl-wrap"><table><thead><tr><th>Result</th><th>Evidence</th><th>Check</th></tr></thead>'
                    f"<tbody>{rows}</tbody></table></div>")
            if d.get("kind") == "bands" and (Path(r["path"]) / "vasprun.xml").exists():
                raw("<h4>Band figure</h4>")
                self.bands_block(r)
                raw("<h4>Figure in the run folder</h4>")
            else:
                raw("<h4>Figure</h4>")
            self.figure(r)
        dlg.open()

    # Band figure dialog opened from the Plot pill.
    def open_bands(self, name):
        r = self.find_run(name)
        if not r:
            return
        self.context[self.project] = name
        self.chat_panel.refresh()
        dlg, body = self.dialog(f"Band figure of {name}")
        with body:
            self.bands_block(r)
        dlg.open()

    # Band figure with energy window inputs and Redraw; draws the default window when no figure exists yet.
    def bands_block(self, r):
        project, name = self.project, r["name"]
        png = bands.existing_figure(r["path"], project, name)
        info = bands.figure_info(png) if png else None
        emin0 = info["emin"] if info else -3.0
        emax0 = info["emax"] if info else 3.0
        with ui.element("div").classes("flex flex-wrap items-end gap-2"):
            lo = ui.number("Lowest energy (eV)", value=emin0, step=0.5, format="%.2f").props("dark outlined dense").classes("w-40")
            hi = ui.number("Highest energy (eV)", value=emax0, step=0.5, format="%.2f").props("dark outlined dense").classes("w-40")
            button("Redraw", "btn primary", lambda: draw(lo.value, hi.value))
        area = ui.element("div").classes("w-full mt-2")

        # Show one PNG with its notes and file paths.
        def show(png_path, res):
            area.clear()
            with area:
                data = base64.b64encode(Path(png_path).read_bytes()).decode()
                raw(f'<img src="data:image/png;base64,{data}" alt="band structure of {esc(name)}">').classes("figure")
                if res:
                    split = res.get("largest_splitting_eV")
                    lines = list(res.get("notes", []))
                    if split is not None:
                        at = res["largest_splitting_at"]
                        path = "sampled along the k-path" if res.get("labels_found") else "sampled"
                        head = "none detected; " if split < SPLIT_NONE_EV else ""
                        lines.insert(0, f"{head}Largest up-down difference among the sampled k-points and bands in the window: "
                                        f"{split:.3f} eV, band {at['band']}, k-point {at['kpoint_index']} "
                                        f"{at['kpoint_label'] or at['kpoint_frac']} ({res.get('n_kpoints')} k-points {path}).")
                    if not res.get("labels_found"):
                        lines.insert(0, "No high-symmetry labels found; the x axis is the k-point index.")
                    raw('<ul class="lines">' + "".join(f"<li>{esc(t)}</li>" for t in lines) + "</ul>")
                    for kind, path in res["paths"].items():
                        raw(f'<p class="note">{kind.upper()}: <span class="mono path">{esc(files.display_path(path))}</span></p>')

        # Draw the figure off the event loop for the given window and show it.
        async def draw(emin, emax):
            if emin is None or emax is None or emin >= emax:
                ui.notify("The lowest energy must be below the highest energy.")
                return
            area.clear()
            with area:
                raw('<p class="note">Drawing the band figure from vasprun.xml...</p>')
            try:
                res = await ng_run.io_bound(bands.plot_bands, r["path"], emin=emin, emax=emax, project=project, run_name=name)
            except Exception as exc:
                area.clear()
                with area:
                    raw(f'<p class="note">Band figure failed: {esc(f"{type(exc).__name__}: {exc}")}</p>')
                return
            show(res["paths"]["png"], res)

        if png:
            show(png, info)
        else:
            ui.timer(0.05, lambda: draw(-3.0, 3.0), once=True)

    # Wannier comparison dialog opened from the Compare pill: cached overlay figure, or computed on first open.
    def open_compare(self, name):
        from . import wannier_compare as wc
        project = self.project
        dlg, body = self.dialog(f"Wannier bands versus VASP bands: {name}")
        with body:
            area = ui.element("div").classes("w-full")

        # Show the band-match line, the figure and its paths.
        def show(res):
            area.clear()
            with area:
                raw(f"<p>{etag(PARSED)}{esc(wc.band_match_text(res))}</p>")
                if res and res.get("png") and Path(res["png"]).exists():
                    data = base64.b64encode(Path(res["png"]).read_bytes()).decode()
                    raw(f'<img src="data:image/png;base64,{data}" alt="Wannier versus VASP bands of {esc(name)}">').classes("figure")
                    raw(f'<p class="note">{esc(res.get("evidence_state", ""))}</p>')
                    for note in res.get("fallbacks", []):
                        raw(f'<p class="note">{esc(note)}</p>')
                    raw(f'<p class="note">PNG: <span class="mono path">{esc(files.display_path(res["png"]))}</span></p>')

        # Compute the comparison off the event loop (it reads the whole Hamiltonian), then show it.
        async def compute():
            res = await ng_run.io_bound(wc.band_match, project, name, True)
            show(res)

        cached = wc.band_match(project, name, compute=False)
        if cached is not None:
            show(cached)
        else:
            with area:
                raw('<p class="note">Reading wannier90_hr.dat and interpolating along the ancestor k-path; '
                    'large Hamiltonians take about a minute...</p>')
            ui.timer(0.05, compute, once=True)
        dlg.open()

    # Structure dialog opened from the Structure pill; draws the figure on first open.
    def open_structure(self, name):
        r = self.find_run(name)
        if not r:
            return
        self.context[self.project] = name
        self.chat_panel.refresh()
        dlg, body = self.dialog(f"Structure of {name}")
        with body:
            area = ui.element("div").classes("w-full")
        project = self.project

        # Show the PNG, the moment table, the notes and a link to the interactive HTML view.
        def show(res):
            area.clear()
            with area:
                png = Path(res["paths"]["png"])
                data = base64.b64encode(png.read_bytes()).decode()
                url = "/figures/" + Path(res["paths"]["html"]).relative_to(bands.FIG_DIR).as_posix()
                raw(f'<p class="caption">Input magnetic configuration from INCAR MAGMOM (starting moments, '
                    'not the converged OUTCAR moments)</p>')
                raw(f'<iframe class="view3d" src="{esc(url)}" title="3D structure of {esc(name)}"></iframe>'
                    f'<p class="note"><a class="link" href="{esc(url)}" target="_blank" rel="noopener">Open the 3D view in a new tab</a>'
                    ' (3Dmol.js needs internet)</p>')
                raw(f'<details class="fallback"><summary>Static PNG fallback</summary><div class="figure">'
                    f'<img src="data:image/png;base64,{data}" alt="structure of {esc(name)}, input magnetic configuration from INCAR MAGMOM">'
                    '</div></details>')
                rows = [(f'{m["element"]}{m["site"]}', f'{m["magmom"]}', " ".join(f"{x:+.2f}" for x in m["vector"]))
                        for m in res["moments"]]
                rows = [(a, f"MAGMOM {b}, vector {c} μB") for a, b, c in rows]
                rows.append(("net", " ".join(f"{x:+.2f}" for x in res["net_moment"]) + " μB"))
                raw(f'<h4>{esc(res["formula"])}, {esc(res["spacegroup"])} ({res["spacegroup_number"]})</h4>'
                    + kv_table(rows, "Site", "Moment from INCAR"))
                raw('<ul class="lines">' + "".join(f"<li>{esc(t)}</li>" for t in res["notes"]) + "</ul>")
                for kind, path in res["paths"].items():
                    raw(f'<p class="note">{kind.upper()}: <span class="mono path">{esc(files.display_path(path))}</span></p>')

        # Draw the figure off the event loop and show it.
        async def draw():
            with area:
                raw('<p class="note">Drawing the structure from POSCAR and INCAR...</p>')
            try:
                res = await ng_run.io_bound(structure.plot_structure, r["path"], project=project, run_name=name)
            except Exception as exc:
                area.clear()
                with area:
                    raw(f'<p class="note">Structure figure failed: {esc(f"{type(exc).__name__}: {exc}")}</p>')
                return
            show(res)

        png = structure.existing_figure(r["path"], project, name)
        info = bands.figure_info(png) if png else None
        if info:
            show(info)
        else:
            ui.timer(0.05, draw, once=True)
        dlg.open()

    # Show the first figure file of a run inline, or a note when there is none.
    def figure(self, r):
        figs = r["details"].get("figures") or []
        img = next((f for f in figs if f.lower().endswith((".png", ".svg"))), None)
        if not img:
            raw(f'<p class="note">{"PDF figure " + esc(figs[0]) + ", open it from the folder." if figs else "No figure file in this run."}</p>')
            return
        p = Path(r["path"]) / img
        mime = "image/svg+xml" if img.lower().endswith(".svg") else "image/png"
        data = base64.b64encode(p.read_bytes()).decode()
        raw(f'<img src="data:{mime};base64,{data}" alt="{esc(img)}"><p class="note" style="color:#5b626b">{esc(img)}</p>').classes("figure")

    # File dialog: name, path, size, date, copy button and excerpt.
    def open_file(self, run_name, fname, action):
        r = self.find_run(run_name)
        if not r:
            return
        path = Path(r["path"]) / fname
        info = files.excerpt(path)
        dlg, body = self.dialog("File")
        with body:
            raw(f'<p><strong>File:</strong> <span class="mono">{esc(fname)}</span></p>'
                f'<p><strong>Run:</strong> <span class="mono">{esc(run_name)}</span></p>')
            self.path_block(str(path))
            raw(f'<p><strong>Size and date:</strong> {info["size_text"]}, {info["mtime_text"]}</p>')
            if action == "ref" or not info["opened"]:
                raw(f'<p class="note">{esc(info.get("note", "Reference only"))}</p>')
            else:
                if action == "tail":
                    raw('<p class="note">Last lines of a file that is still changing.</p>')
                raw("<h4>Excerpt</h4>")
                raw(f'<pre class="excerpt"><code>{esc(files.excerpt_text(info))}</code></pre>')
        dlg.open()

    # Folder dialog: path, copy button and file listing.
    def open_folder(self, run_name):
        r = self.find_run(run_name)
        if not r:
            return
        dlg, body = self.dialog("Folder")
        with body:
            self.path_block(r["path"].rstrip("/") + "/")
            rows = "".join(f'<tr><td class="mono">{esc(f["name"])}</td><td>{f["size_text"]}</td><td>{f["mtime_text"]}</td></tr>'
                           for f in files.folder_listing(r["path"]))
            raw('<h4>Files</h4><div class="tbl-wrap"><table><thead><tr><th>Name</th><th>Size</th><th>Modified</th></tr>'
                f'</thead><tbody>{rows}</tbody></table></div>')
        dlg.open()

    # Import dialog: path field, read-only note, scan and refresh.
    def open_import(self):
        dlg, body = self.dialog("Import folder")
        with body:
            path = ui.input("Folder path").props("dark outlined").classes("w-full")
            name = ui.input("Project name (optional)").props("dark outlined").classes("w-full")
            raw('<p class="note">Read-only scan, nothing is moved.</p>')
            status = raw("")

            # Run the importer off the event loop and refresh the panels.
            async def scan():
                if not path.value.strip():
                    status.set_content('<p class="note">Enter a folder path.</p>')
                    return
                status.set_content('<p class="note">Scanning...</p>')
                from .importer import import_project
                try:
                    pname, n = await ng_run.io_bound(import_project, path.value.strip(), name.value.strip() or None)
                except Exception as exc:
                    status.set_content(f'<p class="note">Import failed: {esc(str(exc))}</p>')
                    return
                self.project = pname
                dlg.close()
                ui.notify(f"Imported {n} runs into {pname}")
                self.refresh_all()

            with ui.element("div").classes("flex gap-2 mt-3"):
                button("Scan", "btn primary", scan)
                button("Cancel", "btn", dlg.close)
        dlg.open()

    # Right panel: context line, messages, input and Send.
    @ui.refreshable_method
    def chat_panel(self):
        ctx = self.context.get(self.project, "project summary") if self.project else "no project"
        raw(f'Agent<span class="ctx">Context: {esc(self.project or "none")}, {esc(ctx)}</span>').classes("panel-head chat-head")
        msgs = self.chats.setdefault(self.project, [])
        with ui.element("div").classes("panel-body"):
            with ui.element("div").classes("msgs"):
                if not msgs:
                    raw('<span class="who">Agent</span>Ask about this project. Answers come only from the indexed files '
                        'and cite run and file names.').classes("msg agent")
                for m in msgs:
                    who = "You" if m["role"] == "user" else "Agent"
                    cls = "user" if m["role"] == "user" else "agent"
                    if cls == "user":
                        raw(f'<span class="who">{who}</span>{esc(m["content"])}').classes(f"msg {cls}")
                    else:
                        with ui.element("div").classes("msg agent mdmsg"):
                            raw('<span class="who">Agent</span>')
                            ui.markdown(m["content"], extras=["fenced-code-blocks", "tables", "cuddled-lists", "code-friendly"])
        with ui.element("div").classes("chat-form"):
            inp = ui.input(placeholder="Ask about this project or phase").props("dark outlined dense").classes("flex-1")

            # Send the question to the agent and show the answer.
            async def send():
                q = (inp.value or "").strip()
                if not q or not self.project:
                    return
                history = list(msgs)
                msgs.append({"role": "user", "content": q})
                self.chat_panel.refresh()
                ans = await ng_run.io_bound(agent.chat, self.project, q, history)
                msgs.append({"role": "assistant", "content": ans})
                self.chat_panel.refresh()

            inp.on("keydown.enter", send)
            button("Send", "btn primary", send)


# Register the page and start the server.
def run(port=8080):
    # Build one fresh page per browser tab.
    @ui.page("/")
    def index(project: str = None):
        Page(project)

    bands.FIG_DIR.mkdir(parents=True, exist_ok=True)
    app.add_static_files("/figures", str(bands.FIG_DIR))
    ui.run(port=port, title="AltermagAgent", dark=True, reload=False, show=False, show_welcome_message=False)
