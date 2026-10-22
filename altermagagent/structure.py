"""Crystal structure figures with magnetic moment arrows from POSCAR (or CONTCAR) and the INCAR MAGMOM line.

Figures go to ~/.altermagagent/figures/<project>/ and never into the run folder.
"""

import json
from pathlib import Path

import numpy as np

from . import db
from .bands import default_out_dir, inside, slug
from .files import display_path

FALLBACK_COLOURS = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c"]
ARROW_COLOUR = "#d62728"
CDN_3DMOL = "https://3Dmol.org/build/3Dmol-min.js"

HTML_PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<!-- 3Dmol.js is loaded from its public CDN, so this page needs an internet connection; offline it shows only the text. -->
<script src="{cdn}"></script>
<style>
body {{ margin:0; padding:16px; font-family: system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif;
  font-size:16px; background:#000; color:#e6e8ec; }}
h1 {{ font-size:1.2rem; margin:0 0 8px; }} p {{ margin:4px 0 12px; color:#9aa3ad; }}
#view {{ position:relative; width:100%; height:70vh; min-height:320px; border:1px solid #262a33; border-radius:8px; }}
ul {{ padding-left:1.2rem; }} li {{ margin:2px 0; }}
.key span {{ display:inline-block; width:0.9em; height:0.9em; border-radius:50%; margin:0 4px 0 12px; vertical-align:middle; }}
@media (max-width: 600px) {{ body {{ padding:8px; }} }}
</style></head><body>
<h1>{title}</h1><p>{subtitle}</p>
<p class="key">{key}</p>
<div id="view"></div>
<ul>{notes}</ul>
<p>This view loads 3Dmol.js from {cdn}; without an internet connection use the PNG or SVG next to this file.</p>
<script>
const data = {data};
if (window.$3Dmol) {{
  const v = $3Dmol.createViewer("view", {{ backgroundColor: "black" }});
  for (const a of data.atoms) {{
    v.addSphere({{ center: {{x:a.x, y:a.y, z:a.z}}, radius: a.r, color: a.colour }});
  }}
  for (const e of data.edges) {{
    v.addCylinder({{ start: {{x:e[0][0], y:e[0][1], z:e[0][2]}}, end: {{x:e[1][0], y:e[1][1], z:e[1][2]}},
                    radius: 0.03, color: "#9aa3ad", fromCap: 1, toCap: 1 }});
  }}
  for (const m of data.arrows) {{
    v.addArrow({{ start: {{x:m[0][0], y:m[0][1], z:m[0][2]}}, end: {{x:m[1][0], y:m[1][1], z:m[1][2]}},
                 radius: 0.08, radiusRatio: 2.2, mid: 0.75, color: "{arrow}" }});
  }}
  v.zoomTo(); v.render();
}} else {{
  document.getElementById("view").textContent = "3Dmol.js could not be loaded (no internet connection?).";
}}
</script>
</body></html>
"""


# Element colours from pymatgen's Jmol scheme, with a fixed fallback palette for unknown elements.
def element_colours(elements):
    try:
        import yaml
        import pymatgen.vis
        table = yaml.safe_load((Path(pymatgen.vis.__file__).parent / "ElementColorSchemes.yaml").read_text())["Jmol"]
    except Exception:
        table = {}
    out = {}
    for i, el in enumerate(elements):
        rgb = table.get(el)
        out[el] = "#%02x%02x%02x" % tuple(rgb) if rgb else FALLBACK_COLOURS[i % len(FALLBACK_COLOURS)]
    return out


# Read the structure of a run from POSCAR, or CONTCAR when POSCAR is missing or unreadable.
def read_structure(run_dir, notes):
    from pymatgen.core import Structure
    for name in ("POSCAR", "CONTCAR"):
        p = Path(run_dir) / name
        if p.exists() and p.stat().st_size > 0:
            try:
                return Structure.from_file(str(p)), name
            except Exception as exc:
                notes.append(f"{name} could not be read ({type(exc).__name__}: {exc}).")
    raise FileNotFoundError(f"No readable POSCAR or CONTCAR in {run_dir}.")


# Rotate a MAGMOM vector given in the SAXIS frame into Cartesian coordinates, following the VASP manual.
def saxis_to_cartesian(m, saxis):
    sx, sy, sz = saxis
    a = np.arctan2(sy, sx)
    b = np.arctan2(np.hypot(sx, sy), sz)
    m1, m2, m3 = m
    return np.array([np.cos(b) * np.cos(a) * m1 - np.sin(a) * m2 + np.sin(b) * np.cos(a) * m3,
                     np.cos(b) * np.sin(a) * m1 + np.cos(a) * m2 + np.sin(b) * np.sin(a) * m3,
                     -np.sin(b) * m1 + np.cos(b) * m3])


# Parse MAGMOM from INCAR into one Cartesian vector per site, with the raw value, SOC flag and SAXIS.
def read_moments(run_dir, n_sites, notes):
    from pymatgen.io.vasp import Incar
    path = Path(run_dir) / "INCAR"
    if not path.exists():
        notes.append("No INCAR in the run folder, so no moments are drawn.")
        return None, [], False, None
    incar = Incar.from_file(str(path))
    soc = bool(incar.get("LSORBIT")) or bool(incar.get("LNONCOLLINEAR"))
    ispin = int(incar.get("ISPIN", 1) or 1)
    saxis = incar.get("SAXIS")
    raw = incar.get("MAGMOM")
    if raw is None:
        if not soc and ispin == 1:
            notes.append("ISPIN = 1 and no MAGMOM: non-magnetic run, no arrows.")
            return None, [], soc, saxis
        raw = [[1.0, 1.0, 1.0]] * n_sites if soc else [1.0] * n_sites
        notes.append("No MAGMOM in INCAR; VASP's default of 1 μB per atom (1 1 1 for SOC) is drawn.")
    raw = list(raw)
    if soc and raw and not isinstance(raw[0], (list, tuple)):
        if len(raw) % 3:
            notes.append(f"SOC run but MAGMOM has {len(raw)} numbers, not a multiple of 3; no arrows drawn.")
            return None, raw, soc, saxis
        raw = [raw[i:i + 3] for i in range(0, len(raw), 3)]
    if len(raw) != n_sites:
        notes.append(f"MAGMOM has {len(raw)} entries but the structure has {n_sites} sites; no arrows drawn.")
        return None, raw, soc, saxis
    if soc:
        vecs = np.array(raw, dtype=float)
        if saxis is not None and not np.allclose(np.array(saxis, float) / np.linalg.norm(saxis), [0, 0, 1]):
            vecs = np.array([saxis_to_cartesian(v, np.array(saxis, float)) for v in vecs])
            notes.append(f"SOC run with SAXIS = {' '.join(f'{x:g}' for x in saxis)}: MAGMOM vectors rotated from the SAXIS frame "
                         "into Cartesian axes.")
        else:
            notes.append("SOC run: MAGMOM gives a 3-vector per atom, Cartesian because SAXIS is the default 0 0 1.")
    else:
        vecs = np.array([[0.0, 0.0, float(m)] for m in raw])
        notes.append("Collinear run: each MAGMOM value is drawn along z; without SOC the spin axis is not tied to the lattice.")
    notes.append("Arrows show the starting moments from INCAR MAGMOM, not the converged moments from OUTCAR.")
    return vecs, raw, soc, saxis


# Shortest distance between two different sites (periodic images included), used to scale the arrows.
def bond_length(structure):
    if len(structure) > 1:
        dm = structure.distance_matrix
        pos = dm[dm > 1e-6]
        if pos.size:
            return float(pos.min())
    return float(min(structure.lattice.abc))


# The twelve edges of the unit cell as pairs of Cartesian points.
def cell_edges(lattice):
    m = lattice.matrix
    corners = {(i, j, k): i * m[0] + j * m[1] + k * m[2] for i in (0, 1) for j in (0, 1) for k in (0, 1)}
    edges = []
    for c, p in corners.items():
        for ax in range(3):
            if c[ax] == 0:
                d = list(c)
                d[ax] = 1
                edges.append((p, corners[tuple(d)]))
    return edges


# Draw the matplotlib 3D figure with atoms, cell edges and moment arrows into PNG and SVG.
def draw_static(structure, vecs, scale, colours, radii, title, png, svg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig = plt.figure(figsize=(7, 7), dpi=150)
    ax = fig.add_subplot(111, projection="3d")
    for a, b in cell_edges(structure.lattice):
        ax.plot(*zip(a, b), color="#888888", lw=1)
    xyz = structure.cart_coords
    els = [s.specie.symbol for s in structure]
    for el in colours:
        idx = [i for i, e in enumerate(els) if e == el]
        ax.scatter(*xyz[idx].T, s=radii[el] ** 2 * 900, color=colours[el], edgecolors="#222222", linewidths=0.6,
                   depthshade=False, alpha=0.9)
    if vecs is not None:
        for p, v in zip(xyz, vecs):
            if np.linalg.norm(v) < 1e-6:
                continue
            d = v * scale
            ax.quiver(*(p - d / 2), *d, color=ARROW_COLOUR, linewidth=2.5, arrow_length_ratio=0.3)
    pts = np.vstack([xyz] + [np.array(e) for e in cell_edges(structure.lattice)])
    lo, hi = pts.min(0) - 0.5, pts.max(0) + 0.5
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(hi - lo)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    handles = [Line2D([], [], marker="o", ls="", markersize=10, markerfacecolor=c, markeredgecolor="#222222", label=el)
               for el, c in colours.items()]
    if vecs is not None:
        handles.append(Line2D([], [], color=ARROW_COLOUR, lw=2.5, label="moment (MAGMOM)"))
    ax.legend(handles=handles, loc="upper left")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(png)
    fig.savefig(svg)
    plt.close(fig)


# Write the interactive 3Dmol.js page with the same atoms, cell edges and arrows.
def draw_html(structure, vecs, scale, colours, radii, title, subtitle, notes, html_path):
    import html as h
    atoms = [{"x": float(p[0]), "y": float(p[1]), "z": float(p[2]), "r": round(radii[s.specie.symbol] * 0.5, 3),
              "colour": colours[s.specie.symbol]} for s, p in zip(structure, structure.cart_coords)]
    edges = [[list(map(float, a)), list(map(float, b))] for a, b in cell_edges(structure.lattice)]
    arrows = []
    if vecs is not None:
        for p, v in zip(structure.cart_coords, vecs):
            if np.linalg.norm(v) > 1e-6:
                d = v * scale
                arrows.append([list(map(float, p - d / 2)), list(map(float, p + d / 2))])
    key = "".join(f'<span style="background:{c}"></span>{h.escape(el)}' for el, c in colours.items())
    if arrows:
        key += f'<span style="background:{ARROW_COLOUR}; border-radius:2px"></span>moment (MAGMOM)'
    page = HTML_PAGE.format(title=h.escape(title), subtitle=h.escape(subtitle), key=key, cdn=CDN_3DMOL, arrow=ARROW_COLOUR,
                            notes="".join(f"<li>{h.escape(n)}</li>" for n in notes),
                            data=json.dumps({"atoms": atoms, "edges": edges, "arrows": arrows}))
    Path(html_path).write_text(page, encoding="utf-8")


# Most recent structure figure PNG of a run in the output folder, or None.
def existing_figure(run_dir, project=None, run_name=None, out_dir=None):
    folder = Path(out_dir) if out_dir else default_out_dir(run_dir, project)
    p = folder / f"{slug(run_name or Path(run_dir).name)}_structure.png"
    return p if p.exists() else None


# Draw the crystal structure of one run with moment arrows into PNG, SVG and HTML and return a result dict.
def plot_structure(run_dir, out_dir=None, supercell=(1, 1, 1), project=None, run_name=None):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    run_dir = Path(run_dir)
    run_name = run_name or run_dir.name
    out = Path(out_dir).expanduser() if out_dir else default_out_dir(run_dir, project)
    if inside(out, run_dir):
        raise ValueError(f"Refusing to write into the run folder {run_dir}; choose an output folder outside it.")
    if project:
        p = db.get_project(project)
        if p and inside(out, p["root"]):
            raise ValueError(f"Refusing to write inside the scanned project folder {p['root']}.")
    notes = []
    structure, sfile = read_structure(run_dir, notes)
    vecs, raw, soc, saxis = read_moments(run_dir, len(structure), notes)
    try:
        sga = SpacegroupAnalyzer(structure, symprec=0.1)
        sg, sgn = sga.get_space_group_symbol(), sga.get_space_group_number()
    except Exception:
        sg, sgn = "unknown", None
    notes.append("Space group from the atom positions only (symprec 0.1 Å); the moments are not included.")
    moments = []
    for i, s in enumerate(structure):
        v = vecs[i] if vecs is not None else np.zeros(3)
        moments.append({"site": i + 1, "element": s.specie.symbol, "magmom": raw[i] if vecs is not None else None,
                        "vector": [round(float(x), 4) for x in v], "size": round(float(np.linalg.norm(v)), 4)})
    net = vecs.sum(0) if vecs is not None else np.zeros(3)
    sc = tuple(int(x) for x in supercell)
    draw = structure.copy()
    dvecs = vecs
    if sc != (1, 1, 1):
        if vecs is not None:
            draw.add_site_property("mvec", [list(v) for v in vecs])
        draw.make_supercell(sc)
        dvecs = np.array(draw.site_properties["mvec"]) if vecs is not None else None
        notes.append(f"Drawn as a {sc[0]}x{sc[1]}x{sc[2]} supercell; moments repeat in every copy.")
    bond = bond_length(structure)
    biggest = float(np.linalg.norm(vecs, axis=1).max()) if vecs is not None else 0.0
    scale = bond / biggest if biggest > 1e-6 else 0.0
    if vecs is not None and biggest <= 1e-6:
        notes.append("All MAGMOM values are zero, so no arrows are drawn.")
    notes.append(f"Arrow length: the largest moment ({biggest:.2f} μB) is drawn as the shortest bond, {bond:.2f} Å.")
    from pymatgen.core.periodic_table import Element
    els = list(dict.fromkeys(s.specie.symbol for s in structure))
    colours = element_colours(els)
    radii = {}
    for el in els:
        r = Element(el).atomic_radius
        radii[el] = min(float(r) if r else 1.2, 0.35 * bond)
    out.mkdir(parents=True, exist_ok=True)
    tag = "" if sc == (1, 1, 1) else f"_{sc[0]}x{sc[1]}x{sc[2]}"
    stem = str(out / f"{slug(run_name)}_structure{tag}")
    paths = {ext: stem + "." + ext for ext in ("png", "svg", "html")}
    formula = structure.composition.reduced_formula
    head = f"{project}, {run_name}" if project else run_name
    title = f"{head}: {formula}, {sg} ({sgn})"
    draw_static(draw, dvecs, scale, colours, radii, title, paths["png"], paths["svg"])
    sub = f"Structure from {display_path(run_dir / sfile)}, moments from {display_path(run_dir / 'INCAR')} MAGMOM"
    draw_html(draw, dvecs, scale, colours, radii, title, sub, notes, paths["html"])
    result = {"paths": paths, "formula": formula, "spacegroup": sg, "spacegroup_number": sgn,
              "structure_file": str(run_dir / sfile), "soc": soc, "saxis": list(saxis) if saxis is not None else None,
              "supercell": list(sc), "moments": moments, "net_moment": [round(float(x), 4) for x in net],
              "net_moment_size": round(float(np.linalg.norm(net)), 4), "notes": notes}
    Path(stem + ".json").write_text(json.dumps(result, indent=1, ensure_ascii=False), encoding="utf-8")
    return result
