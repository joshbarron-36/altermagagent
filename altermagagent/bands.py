"""Band-structure figures from a VASP band run, written under ~/.altermagagent/figures and never into the run folder."""

import json
import re
from pathlib import Path

import numpy as np

from . import db
from .files import display_path

FIG_DIR = db.HOME_DIR / "figures"
STYLES = ("paper", "screen")
UP_COLOUR, DOWN_COLOUR, ONE_COLOUR = "#1f4e9c", "#c0392b", "#1a1a1a"
GREEK = {"GAMMA": "Γ", "G": "Γ", "\\GAMMA": "Γ", "Γ": "Γ", "SIGMA": "Σ", "\\SIGMA": "Σ", "DELTA": "Δ", "\\DELTA": "Δ",
         "LAMBDA": "Λ", "\\LAMBDA": "Λ"}
SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
UNSUBSCRIPT = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")

HTML_PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
body {{ margin:0; padding:16px; font-family: Georgia, "Times New Roman", serif; font-size:16px; background:#fff; color:#111; }}
h1 {{ font-size:1.2rem; margin:0 0 8px; }} p {{ margin:4px 0 12px; color:#444; }}
.plot {{ width:100%; height:70vh; min-height:320px; }}
@media (max-width: 600px) {{ body {{ padding:8px; }} }}
</style></head><body>
<h1>{title}</h1><p>{subtitle}</p>
<div class="plot">{plot}</div>
</body></html>
"""


MAX_LABEL = 6
MERGE_FRACTION = 0.02
TICK_FONT_SIZE = 12


# Strip markup and the leading k-point weight that pymatgen or VASPKIT glue in front, for example "1    GAMMA" -> "GAMMA".
def clean_part(part):
    part = (part or "").strip().strip("$").replace("\\", "").replace("{", "").replace("}", "").strip()
    return re.sub(r"^\d+\s*(?=[^\d\s])", "", part).strip()


# Printed symbol of one cleaned label part: GAMMA -> Γ, SIGMA_0 -> Σ₀, Y1 -> Y₁; a bare number gives "" unless keep_numbers.
def pretty_part(part, keep_numbers=False):
    part = clean_part(part)
    if not part:
        return ""
    if part.isdigit():
        return part if keep_numbers else ""
    m = re.fullmatch(r"([^\d_\s]+?)_?(\d*)_?", part)
    if not m:
        return part[:MAX_LABEL]
    base, sub = m.group(1), m.group(2)
    base = GREEK.get(base.upper(), base)
    return (base + sub.translate(SUBSCRIPT))[:MAX_LABEL]


# Turn a KPOINTS or KLABELS label such as GAMMA, SIGMA_0 or "Y|1GAMMA" into its printed form, for example Γ, Σ₀ or Y|Γ.
def pretty_label(label, keep_numbers=False):
    out = []
    for part in str(label or "").split("|"):
        lab = pretty_part(part, keep_numbers)
        if lab and (not out or out[-1] != lab):
            out.append(lab)
    return "|".join(out)


# Join two printed labels with "|", dropping repeats so Γ and Γ give Γ and Y|Γ and Γ give Y|Γ.
def join_labels(a, b):
    parts = [p for p in a.split("|") if p]
    for p in b.split("|"):
        if p and (not parts or parts[-1] != p):
            parts.append(p)
    return "|".join(parts)


# Group (x, label) entries by position and merge ticks closer than 2 percent of the axis width into one "A|B" tick.
def build_ticks(entries, width):
    ticks = []
    tol = max(width, 1e-9) * 1e-4
    for x, lab in sorted(entries, key=lambda e: e[0]):
        if ticks and abs(x - ticks[-1][0]) <= tol:
            ticks[-1] = (ticks[-1][0], join_labels(ticks[-1][1], lab))
        else:
            ticks.append((float(x), lab))
    merged = []
    for x, lab in ticks:
        if merged and x - merged[-1][0] < MERGE_FRACTION * width:
            merged[-1] = (merged[-1][0], join_labels(merged[-1][1], lab))
        else:
            merged.append((x, lab))
    return [(x, lab) for x, lab in merged if lab]


# Read a VASPKIT KLABELS file into (x, raw label) pairs, or None when it is missing or unreadable.
def read_klabels(path):
    if not path.exists():
        return None
    rows = []
    try:
        lines = path.read_text(errors="replace").splitlines()[1:]
    except OSError:
        return None
    for line in lines:
        bits = line.split()
        if not bits:
            break
        try:
            x = float(bits[-1])
        except ValueError:
            break
        rows.append((x, " ".join(bits[:-1])))
    return rows or None


# Ticks from KLABELS, scaled to the plotted distance; truncated names such as SIGMA_ are completed from KPOINTS labels.
def klabels_ticks(rows, x_end, kp_entries, notes):
    scale = 1.0
    if rows[-1][0] > 0 and abs(rows[-1][0] - x_end) > 0.01 * x_end:
        scale = x_end / rows[-1][0]
        notes.append(f"KLABELS distances rescaled by {scale:.4f} to match the plotted path length.")
    tol = max(x_end, 1e-9) * 1e-3
    entries = []
    for x, raw in rows:
        x = x * scale
        near = [clean_part(p) for kx, kr in kp_entries if abs(kx - x) <= tol for p in str(kr).split("|")]
        parts = []
        for part in raw.split("|"):
            c = clean_part(part)
            full = next((n for n in near if len(n) > len(c) and n.startswith(c)), c) if c else c
            parts.append(full)
        entries.append((x, pretty_label("|".join(parts), keep_numbers=True)))
    return build_ticks(entries, x_end)


# Matplotlib form of a printed label, with subscript digits drawn as mathtext subscripts.
def mpl_label(label):
    return re.sub(r"([₀-₉]+)", lambda m: "$_{" + m.group(1).translate(UNSUBSCRIPT) + "}$", label)


# Safe file-name part from a run or project name.
def slug(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip("/")) or "run"


# Default output folder: ~/.altermagagent/figures/<project or run name>.
def default_out_dir(run_dir, project=None):
    return FIG_DIR / slug(project or Path(run_dir).name)


# File stem of one figure, built from run name, style and energy window.
def figure_stem(run_name, style, emin, emax):
    return f"{slug(run_name)}_bands_{style}_{emin:g}_{emax:g}"


# Most recent band figure PNG of a run in the output folder, or None.
def existing_figure(run_dir, project=None, run_name=None, out_dir=None):
    folder = Path(out_dir) if out_dir else default_out_dir(run_dir, project)
    found = sorted(folder.glob(f"{slug(run_name or Path(run_dir).name)}_bands_*.png"), key=lambda p: p.stat().st_mtime)
    return found[-1] if found else None


# Stored result dict next to a figure PNG, or None when there is none.
def figure_info(png_path):
    side = Path(png_path).with_suffix(".json")
    return json.loads(side.read_text()) if side.exists() else None


# Return True when path lies inside folder (or is the folder itself).
def inside(path, folder):
    path, folder = Path(path).resolve(), Path(folder).resolve()
    return path == folder or folder in path.parents


# Raise ValueError(message) when out lies inside any of the roots; {root} in the message names the folder.
def refuse_inside_root(out, roots, message):
    for root in roots:
        if inside(out, root):
            raise ValueError(message.format(root=root))


# Read KPOINTS and say whether it is line mode with labels, without raising.
def kpoints_line_labels(kpath):
    from pymatgen.io.vasp import Kpoints
    if not kpath.exists():
        return False, "No KPOINTS file in the run folder."
    try:
        kp = Kpoints.from_file(kpath)
    except Exception as exc:
        return False, f"KPOINTS could not be read ({type(exc).__name__})."
    if "line" not in str(kp.style).lower():
        return False, f"KPOINTS is not in line mode (style {kp.style})."
    if not any(kp.labels or []):
        return False, "KPOINTS is in line mode but has no labels."
    return True, ""


# Load energies, x positions, branches and tick labels from vasprun.xml and KPOINTS.
def load_bands(run_dir, notes):
    from pymatgen.electronic_structure.core import Spin
    from pymatgen.io.vasp import Vasprun
    run_dir = Path(run_dir)
    v = Vasprun(str(run_dir / "vasprun.xml"), parse_projected_eigen=False, parse_potcar_file=False)
    efermi = float(v.efermi)
    par = v.parameters
    soc = bool(par.get("LSORBIT") or par.get("LNONCOLLINEAR"))
    spin_pol = bool(v.is_spin) and not soc
    data = {"efermi": efermi, "soc": soc, "spin_polarised": spin_pol, "labels_found": False}
    ok, why = kpoints_line_labels(run_dir / "KPOINTS")
    bs = None
    if ok:
        try:
            bs = v.get_band_structure(kpoints_filename=str(run_dir / "KPOINTS"), line_mode=True, efermi=efermi)
            if not getattr(bs, "branches", None):
                raise ValueError("no branches in the band structure")
        except Exception as exc:
            bs, why = None, f"Line-mode band structure failed ({type(exc).__name__}: {exc})."
    if bs is not None:
        data["labels_found"] = True
        data["x"] = np.array(bs.distance, dtype=float)
        data["bands"] = {s: np.array(bs.bands[s]) - efermi for s in bs.bands}
        data["branches"] = [(b["start_index"], b["end_index"]) for b in bs.branches]
        data["kfrac"] = [k.frac_coords for k in bs.kpoints]
        data["klabel"] = [k.label for k in bs.kpoints]
        kp_entries = []
        for b in bs.branches:
            for i, name in ((b["start_index"], b["name"].split("-")[0]), (b["end_index"], b["name"].split("-")[-1])):
                kp_entries.append((float(bs.distance[i]), bs.kpoints[i].label or name))
        width = float(data["x"][-1] - data["x"][0])
        rows = read_klabels(run_dir / "KLABELS")
        if rows:
            data["ticks"] = klabels_ticks(rows, float(data["x"][-1]), kp_entries, notes)
            notes.append("labels from KLABELS (VASPKIT)")
        else:
            data["ticks"] = build_ticks([(x, pretty_label(lab)) for x, lab in kp_entries], width)
            notes.append("labels from KPOINTS line mode")
    else:
        notes.append(why + " Plotted against the plain k-point index instead of the high-symmetry path.")
        eig = v.eigenvalues
        data["bands"] = {s: eig[s][:, :, 0].T - efermi for s in eig}
        nk = next(iter(data["bands"].values())).shape[1]
        data["x"] = np.arange(nk, dtype=float)
        data["branches"] = [(0, nk - 1)]
        data["kfrac"] = [np.array(k) for k in v.actual_kpoints]
        data["klabel"] = [None] * nk
        data["ticks"] = []
    if soc or not spin_pol:
        data["bands"] = {Spin.up: data["bands"][Spin.up]}
    data["spin_up"], data["spin_down"] = Spin.up, Spin.down
    return data


# Largest up-down difference per band inside the window and the overall largest, with its 1-based k-point index.
def spin_splitting(data, emin, emax):
    up, dn = data["bands"][data["spin_up"]], data["bands"][data["spin_down"]]
    nb = min(up.shape[0], dn.shape[0])
    per_band, best = [], None
    for i in range(nb):
        mask = ((up[i] >= emin) & (up[i] <= emax)) | ((dn[i] >= emin) & (dn[i] <= emax))
        if not mask.any():
            continue
        diff = np.where(mask, np.abs(up[i] - dn[i]), -1.0)
        k = int(np.argmax(diff))
        rec = {"band": i + 1, "max_eV": round(float(diff[k]), 4), "kpoint_index": k + 1}
        per_band.append(rec)
        if best is None or rec["max_eV"] > best["max_eV"]:
            best = rec
    if best is None:
        return None, None, []
    k = best["kpoint_index"] - 1
    where = {"band": best["band"], "kpoint_index": k + 1,
             "kpoint_frac": [round(float(c), 4) for c in data["kfrac"][k]],
             "kpoint_label": pretty_label(data["klabel"][k]) if data["klabel"][k] else None,
             "E_up_eV": round(float(up[best["band"] - 1][k]), 4), "E_down_eV": round(float(dn[best["band"] - 1][k]), 4)}
    return best["max_eV"], where, per_band


# Matplotlib settings of one style: paper is white, serif, thin lines and 300 dpi.
def style_rc(style):
    if style == "paper":
        return {"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 9,
                "mathtext.fontset": "dejavuserif",
                "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"}, 0.7, 300, (3.4, 3.0)
    return {"font.family": "sans-serif", "font.size": 12, "figure.facecolor": "white", "axes.facecolor": "white",
            "savefig.facecolor": "white"}, 1.3, 150, (7.0, 5.0)


# Line segments of one band, split at path breaks so separate branches are not joined.
def segments(x, y, branches):
    return [(x[a:b + 1], y[a:b + 1]) for a, b in branches]


# Shared decoration of a band panel: tick lines, dotted zero line, k-point labels, x range, y range and energy label.
def decorate_band_axes(axis, data, ylim, zero_colour):
    for position, _ in data["ticks"]:
        axis.axvline(position, color="0.85", lw=0.6, zorder=0)
    axis.axhline(0, color=zero_colour, ls=":", lw=0.8)
    axis.set_xticks([x for x, _ in data["ticks"]], [mpl_label(label) for _, label in data["ticks"]])
    axis.set(xlim=(data["x"][0], data["x"][-1]), ylim=ylim, ylabel=r"$E-E_F$ (eV)")


# Return True when any two neighbouring x tick labels overlap on the drawn figure.
def ticks_overlap(fig, ax):
    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in ax.get_xticklabels() if t.get_text()]
    return any(a.x1 + 2 > b.x0 for a, b in zip(boxes, boxes[1:]))


# Widen the figure step by step up to 7 inches until the x tick labels no longer overlap, then lay it out.
def fit_tick_labels(fig, ax, size, max_width=7.0):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    FigureCanvasAgg(fig)
    width = size[0]
    fig.tight_layout()
    while ticks_overlap(fig, ax) and width < max_width:
        width = min(max_width, width * 1.2)
        fig.set_size_inches(width, size[1])
        fig.tight_layout()


# Draw the PNG and SVG with matplotlib.
def draw_static(data, emin, emax, style, title, png, svg):
    from matplotlib import rc_context
    from matplotlib.figure import Figure
    rc, lw, dpi, size = style_rc(style)
    with rc_context(rc):
        fig = Figure(figsize=size)
        ax = fig.add_subplot(111)
        x = data["x"]
        two = data["spin_polarised"]
        for s, arr in data["bands"].items():
            colour = (UP_COLOUR if s == data["spin_up"] else DOWN_COLOUR) if two else ONE_COLOUR
            for i in range(arr.shape[0]):
                for j, (xs, ys) in enumerate(segments(x, arr[i], data["branches"])):
                    lab = None
                    if two and i == 0 and j == 0:
                        lab = "spin up" if s == data["spin_up"] else "spin down"
                    ax.plot(xs, ys, color=colour, lw=lw, label=lab)
        ax.axhline(0.0, color="#555555", lw=lw * 0.8, ls="--")
        for tx, _ in data["ticks"]:
            ax.axvline(tx, color="#999999", lw=lw * 0.6)
        if data["ticks"]:
            font = {"fontsize": TICK_FONT_SIZE, **({"fontfamily": "DejaVu Serif"} if style == "paper" else {})}
            ax.set_xticks([t for t, _ in data["ticks"]], [mpl_label(lab) for _, lab in data["ticks"]], **font)
        else:
            ax.set_xlabel("k-point index")
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(emin, emax)
        ax.set_ylabel(r"$E - E_\mathrm{F}$ (eV)")
        ax.set_title(title)
        if two:
            ax.legend(loc="best", frameon=True, framealpha=0.9, edgecolor="none", fontsize="small")
        fit_tick_labels(fig, ax, size)
        fig.savefig(png, dpi=dpi)
        fig.savefig(svg)


# Write the interactive HTML with plotly; hover shows band index and energy.
def draw_html(data, emin, emax, title, subtitle, html_path):
    import plotly.graph_objects as go
    x = data["x"]
    two = data["spin_polarised"]
    fig = go.Figure()
    for s, arr in data["bands"].items():
        name = ("spin up" if s == data["spin_up"] else "spin down") if two else "band"
        colour = (UP_COLOUR if s == data["spin_up"] else DOWN_COLOUR) if two else ONE_COLOUR
        for i in range(arr.shape[0]):
            if arr[i].max() < emin - 1 or arr[i].min() > emax + 1:
                continue
            xs, ys = [], []
            for a, b in data["branches"]:
                xs += list(x[a:b + 1]) + [None]
                ys += list(arr[i][a:b + 1]) + [None]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line={"color": colour, "width": 1.2}, name=name,
                                     legendgroup=name, showlegend=two and i == 0,
                                     customdata=[i + 1] * len(xs),
                                     hovertemplate=f"{name}, band %{{customdata}}<br>E - E_F = %{{y:.4f}} eV<extra></extra>"))
    fig.add_hline(y=0.0, line={"dash": "dash", "color": "#555", "width": 1})
    for tx, _ in data["ticks"]:
        fig.add_vline(x=tx, line={"color": "#bbb", "width": 1})
    xaxis = {"range": [float(x[0]), float(x[-1])], "showgrid": False}
    if data["ticks"]:
        xaxis.update(tickmode="array", tickvals=[t for t, _ in data["ticks"]], ticktext=[lab for _, lab in data["ticks"]],
                     tickfont={"size": TICK_FONT_SIZE, "family": "DejaVu Serif, Georgia, serif"})
    else:
        xaxis["title"] = "k-point index"
    fig.update_layout(template="simple_white", font={"family": "Georgia, Times New Roman, serif"}, xaxis=xaxis,
                      yaxis={"range": [emin, emax], "title": "E - E_F (eV)"}, margin={"l": 60, "r": 20, "t": 20, "b": 40},
                      autosize=True, hovermode="closest")
    plot = fig.to_html(full_html=False, include_plotlyjs=True, config={"responsive": True},
                       default_width="100%", default_height="100%")
    Path(html_path).write_text(HTML_PAGE.format(title=title, subtitle=subtitle, plot=plot), encoding="utf-8")


# Draw the band structure of one run into PNG, SVG and HTML and return a result dict.
def plot_bands(run_dir, emin=-3.0, emax=3.0, style="paper", out_dir=None, project=None, run_name=None):
    run_dir = Path(run_dir)
    emin, emax = float(emin), float(emax)
    if emin >= emax:
        raise ValueError(f"emin ({emin}) must be lower than emax ({emax}).")
    if style not in STYLES:
        raise ValueError(f"Unknown style {style}; use one of {', '.join(STYLES)}.")
    if not (run_dir / "vasprun.xml").exists():
        raise FileNotFoundError(f"No vasprun.xml in {run_dir}.")
    run_name = run_name or run_dir.name
    out = Path(out_dir).expanduser() if out_dir else default_out_dir(run_dir, project)
    refuse_inside_root(out, [run_dir], "Refusing to write into the run folder {root}; choose an output folder outside it.")
    if project:
        p = db.get_project(project)
        refuse_inside_root(out, [p["root"]] if p else [], "Refusing to write inside the scanned project folder {root}.")
    notes = []
    data = load_bands(run_dir, notes)
    nb, nk = next(iter(data["bands"].values())).shape
    split, where, per_band = (None, None, [])
    if data["spin_polarised"]:
        split, where, per_band = spin_splitting(data, emin, emax)
        if split is None:
            notes.append("No band falls inside the energy window, so no up-down difference was measured.")
    if data["soc"]:
        notes.append("SOC run (LSORBIT or LNONCOLLINEAR): spin is not a good quantum number, one colour for all bands.")
    elif not data["spin_polarised"]:
        notes.append("Non-spin-polarised run: one colour for all bands.")
    notes.append(f"Energies relative to the Fermi level {data['efermi']:.4f} eV from vasprun.xml.")
    out.mkdir(parents=True, exist_ok=True)
    stem = str(out / figure_stem(run_name, style, emin, emax))
    paths = {ext: stem + "." + ext for ext in ("png", "svg", "html")}
    title = f"{project}, {run_name}" if project else run_name
    draw_static(data, emin, emax, style, title, paths["png"], paths["svg"])
    sub = (f"{'Spin up and spin down' if data['spin_polarised'] else 'All bands'}, window {emin:g} to {emax:g} eV, "
           f"source {display_path(run_dir / 'vasprun.xml')}")
    draw_html(data, emin, emax, title, sub, paths["html"])
    result = {"paths": paths, "emin": emin, "emax": emax, "style": style, "spin_polarised": data["spin_polarised"],
              "soc": data["soc"], "labels_found": data["labels_found"],
              "labels": [lab for _, lab in data["ticks"]], "n_kpoints": int(nk), "n_bands": int(nb),
              "efermi_eV": round(data["efermi"], 4), "largest_splitting_eV": split, "largest_splitting_at": where,
              "splitting_per_band": per_band, "source": str(run_dir / "vasprun.xml"), "notes": notes}
    Path(stem + ".json").write_text(json.dumps(result, indent=1, ensure_ascii=False), encoding="utf-8")
    return result
