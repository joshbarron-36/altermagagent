"""Magnetic space group, spin space group and wannhr_symm operations of one run, read-only.

Every label here is "parsed, pending validation" until a physicist confirms it.
Nothing is written into the run folder; only CONTCAR or POSCAR, INCAR, OUTCAR and the small wannier_sym text files are read.
"""

import re
from pathlib import Path

import numpy as np

from .importer import expand_magmom, iter_lines, parse_incar, read_moments, tag_value

PENDING = "parsed, pending validation"
SYMPREC = 0.01          # spglib position tolerance, Angstrom
MAG_SYMPREC = 0.05      # moment tolerance, Bohr magneton
TRANS_TOL = 1e-3        # tolerance on fractional translations when comparing operation sets
ROT_WORDS = {(1, 3): "identity", (-1, -3): "inversion", (1, -1): "2-fold rotation", (1, 0): "3-fold rotation",
             (1, 1): "4-fold rotation", (1, 2): "6-fold rotation", (-1, 1): "mirror", (-1, 0): "-3 rotoinversion",
             (-1, -1): "-4 rotoinversion", (-1, -2): "-6 rotoinversion"}
ORDER = {"identity": 1, "inversion": 2, "2-fold rotation": 2, "3-fold rotation": 3, "4-fold rotation": 4,
         "6-fold rotation": 6, "mirror": 2, "-3 rotoinversion": 6, "-4 rotoinversion": 4, "-6 rotoinversion": 6}


# Read the structure from CONTCAR, or POSCAR when CONTCAR is missing or unreadable, and return (structure, file name).
def read_structure(run_dir):
    from pymatgen.core import Structure
    for name in ("CONTCAR", "POSCAR"):
        p = Path(run_dir) / name
        if p.exists() and p.stat().st_size > 0:
            try:
                return Structure.from_file(str(p)), name
            except Exception:
                continue
    raise FileNotFoundError(f"no readable CONTCAR or POSCAR in {run_dir}")


# Spin settings of a run from INCAR: (ISPIN, LSORBIT or LNONCOLLINEAR on).
def spin_settings(run_dir):
    p = Path(run_dir) / "INCAR"
    raw = parse_incar(p) if p.exists() else {}
    ispin = tag_value(raw.get("ISPIN")) if "ISPIN" in raw else 1
    soc = tag_value(raw.get("LSORBIT")) is True or tag_value(raw.get("LNONCOLLINEAR")) is True
    return int(ispin or 1), soc


# Converged moments from the last OUTCAR magnetization blocks: (n x 3 vectors, "vector") or (n scalars, "scalar").
def read_converged_moments(run_dir, n_sites):
    folder = Path(run_dir)
    mom = read_moments(folder, {"OUTCAR"} if (folder / "OUTCAR").exists() else set())
    if not mom:
        raise ValueError(f"no magnetization block in {folder / 'OUTCAR'}")
    if len(mom["sites"]) != n_sites:
        raise ValueError(f"OUTCAR lists {len(mom['sites'])} sites but the structure has {n_sites}")
    if len(mom.get("components", [])) > 1:
        return np.array([s["vector"] for s in mom["sites"]], dtype=float), "vector"
    return np.array([s["moment"] for s in mom["sites"]], dtype=float), "scalar"


# Starting moments from INCAR MAGMOM: (n x 3 Cartesian vectors, "vector") rotated from SAXIS, or (n scalars, "scalar").
def read_incar_moments(run_dir, n_sites):
    p = Path(run_dir) / "INCAR"
    raw = parse_incar(p) if p.exists() else {}
    vals = expand_magmom(raw.get("MAGMOM")) or []
    if len(vals) == n_sites:
        return np.array(vals, dtype=float), "scalar"
    if len(vals) == 3 * n_sites:
        from .structure import saxis_to_cartesian
        sax = [float(x) for x in str(raw.get("SAXIS", "0 0 1")).split()[:3]]
        sax = sax if len(sax) == 3 else [0.0, 0.0, 1.0]
        return np.array([saxis_to_cartesian(vals[3 * i:3 * i + 3], sax) for i in range(n_sites)], dtype=float), "vector"
    raise ValueError(f"INCAR MAGMOM has {len(vals)} values for {n_sites} sites")


# Point-group operation type of an integer rotation in the lattice basis, as words such as "mirror".
def rotation_type(rot):
    r = np.array(rot)
    return ROT_WORDS.get((int(round(np.linalg.det(r))), int(round(np.trace(r)))), "unknown rotation")


# Intrinsic (screw or glide) part of the translation, (1/n) sum of R^k t for k < n, reduced to [-0.5, 0.5).
def intrinsic_translation(rot, trans, order):
    r, t = np.array(rot, dtype=float), np.array(trans, dtype=float)
    acc, rk = np.zeros(3), np.eye(3)
    for _ in range(order):
        acc += rk @ t
        rk = r @ rk
    w = acc / order
    return w - np.round(w)


# One operation in words, such as "2-fold screw (2_1)", "glide (mirror with in-plane translation)" or "inversion".
def describe(rot, trans):
    kind = rotation_type(rot)
    t = np.array(trans, dtype=float)
    t = t - np.round(t)
    if kind == "identity":
        return "identity" if np.allclose(t, 0, atol=TRANS_TOL) else "pure translation"
    w = intrinsic_translation(rot, trans, ORDER.get(kind, 1))
    if np.allclose(w, 0, atol=TRANS_TOL):
        return kind
    if kind.endswith("-fold rotation"):
        return kind.replace("rotation", "screw")
    if kind == "mirror":
        return "glide"
    return kind + " with translation"


# Operation list of rotations, translations and time-reversal flags, with words and unitary or antiunitary marks.
def operation_list(rotations, translations, time_reversals):
    ops = []
    for i, (r, t, tr) in enumerate(zip(rotations, translations, time_reversals), 1):
        tt = np.remainder(np.array(t, dtype=float), 1.0)
        tt[np.isclose(tt, 1.0, atol=TRANS_TOL)] = 0.0
        ops.append({"index": i, "rotation": np.array(r).astype(int).tolist(), "translation": [round(float(x), 5) for x in tt],
                    "time_reversal": bool(tr), "kind": "antiunitary" if tr else "unitary", "words": describe(r, t)})
    return ops


# BNS and OG symbols for a BNS number from pymatgen's ISO-MAG table, or (None, None) when it is not found.
def msg_symbols(bns_number):
    try:
        import sqlite3
        from pymatgen.symmetry.maggroups import MAGSYMM_DATA
        a, b = (int(x) for x in str(bns_number).split("."))
        with sqlite3.connect(MAGSYMM_DATA) as con:
            row = con.execute("SELECT BNS_label, OG_label FROM space_groups WHERE BNS1=? AND BNS2=?", (a, b)).fetchone()
        return (row[0], row[1]) if row else (None, None)
    except Exception:
        return None, None


# Plain space group (moments ignored) from spglib at SYMPREC: symbol, number and operations.
def plain_space_group(structure):
    import spglib
    cell = (structure.lattice.matrix, structure.frac_coords, [s.specie.Z for s in structure])
    ds = spglib.get_symmetry_dataset(cell, symprec=SYMPREC)
    ops = operation_list(ds.rotations, ds.translations, [False] * len(ds.rotations))
    return {"symbol": ds.international, "number": int(ds.number), "n_operations": len(ops), "operations": ops}


# Magnetic space group of one run from its structure and its converged OUTCAR moments, through spglib.
def magnetic_space_group(run_dir):
    import spglib
    run_dir = Path(run_dir)
    structure, sfile = read_structure(run_dir)
    if (run_dir / "OUTCAR").exists():
        moms, form = read_converged_moments(run_dir, len(structure))
        msrc = f"the last magnetization block of {run_dir.name}/OUTCAR"
    else:
        moms, form = read_incar_moments(run_dir, len(structure))
        msrc = f"the starting moments in {run_dir.name}/INCAR MAGMOM (no OUTCAR)"
    ispin, soc = spin_settings(run_dir)
    cell = (structure.lattice.matrix, structure.frac_coords, [s.specie.Z for s in structure], moms)
    ds = spglib.get_magnetic_symmetry_dataset(cell, is_axial=(form == "vector"), symprec=SYMPREC, mag_symprec=MAG_SYMPREC)
    if ds is None:
        raise ValueError("spglib returned no magnetic symmetry dataset")
    t = spglib.get_magnetic_spacegroup_type(ds.uni_number)
    bns_sym, og_sym = msg_symbols(t.bns_number)
    ops = operation_list(ds.rotations, ds.translations, ds.time_reversals)
    n_anti = sum(o["time_reversal"] for o in ops)
    return {"state": PENDING, "run": run_dir.name, "structure_file": sfile, "moment_form": form, "ispin": ispin, "soc": soc,
            "moments": moms.tolist(), "uni_number": int(ds.uni_number), "bns_number": t.bns_number, "og_number": t.og_number,
            "bns_symbol": bns_sym, "og_symbol": og_sym, "litvin_number": t.litvin_number, "type": int(t.type),
            "n_operations": len(ops), "n_unitary": len(ops) - n_anti, "n_antiunitary": n_anti, "operations": ops,
            "space_group": plain_space_group(structure), "symprec": SYMPREC, "mag_symprec": MAG_SYMPREC,
            "moment_source": msrc, "source": f"{run_dir.name}/{sfile} and {msrc}, "
                      f"spglib {spglib.__version__} get_magnetic_symmetry_dataset (symprec {SYMPREC} Å, moment tolerance "
                      f"{MAG_SYMPREC} μB, {'axial vectors' if form == 'vector' else 'collinear scalars'})"}


# Point group symbol of a set of integer rotations through spglib, or "?".
def point_group_symbol(rotations):
    import spglib
    try:
        uniq = {tuple(np.array(r).astype(int).ravel()) for r in rotations}
        return spglib.get_pointgroup(np.array([np.array(u).reshape(3, 3) for u in uniq]))[0].strip()
    except Exception:
        return "?"


# Spin space group of a collinear non-SOC run through spinspg; "not computed" with the reason when it cannot run.
def spin_space_group(run_dir):
    run_dir = Path(run_dir)
    try:
        ispin, soc = spin_settings(run_dir)
        if ispin != 2 or soc:
            return {"state": "not computed", "error": f"spin space group only for collinear non-SOC runs (ISPIN = {ispin}, "
                                                      f"SOC {'on' if soc else 'off'})"}
        import spinspg
        structure, sfile = read_structure(run_dir)
        moms, form = read_converged_moments(run_dir, len(structure))
        vec = np.array([[0.0, 0.0, m] for m in moms]) if form == "scalar" else moms
        sog, rots, trans, spin_rots = spinspg.get_spin_symmetry(structure.lattice.matrix, structure.frac_coords,
                                                               [s.specie.Z for s in structure], vec, symprec=SYMPREC,
                                                               mag_symprec=MAG_SYMPREC)
        axis = np.array([0.0, 0.0, 1.0]) if form == "scalar" else vec[np.argmax(np.linalg.norm(vec, axis=1))]
        axis = axis / np.linalg.norm(axis)
        flips = [np.dot(sr @ axis, axis) < 0 for sr in spin_rots]
        flip_words = sorted({describe(r, t) for r, t, f in zip(rots, trans, flips) if f})
        keep_rots = [r for r, f in zip(rots, flips) if not f]
        g_pg, h_pg = point_group_symbol(rots), point_group_symbol(keep_rots)
        by_trans_inv = [w for w in flip_words if w in ("pure translation", "inversion")]
        relation = ("by translation or inversion (" + ", ".join(by_trans_inv) + ")" if by_trans_inv else
                    "by rotation, mirror, screw, glide or rotoinversion, not by translation or inversion (" + ", ".join(flip_words) + ")"
                    if flip_words else "no spin-flip operation found")
        return {"state": PENDING, "run": run_dir.name, "spin_only_type": str(sog.spin_only_group_type).lower(),
                "label": f"spatial point group {g_pg}, spin-preserving subgroup {h_pg}",
                "n_operations": len(rots), "n_spin_flip": int(sum(flips)), "flip_operations": flip_words,
                "opposite_spin_relation": relation,
                "source": f"{run_dir.name}/{sfile} and {run_dir.name}/OUTCAR moments, spinspg {spinspg.__version__} "
                          f"get_spin_symmetry (symprec {SYMPREC} Å, moment tolerance {MAG_SYMPREC} μB); spinspg gives "
                          "operations, not a label, so the label is built here as G (all) and H (spin-preserving)"}
    except Exception as exc:
        return {"state": "not computed", "error": f"{type(exc).__name__}: {exc}"}


# wannhr_symm output format, learned from MnTe 9.6-wannier/wannier_sym (Diana's run):
#   mysymmop.dat: one block per space-group operation, "  ---------------    N ---------------", then "rotation:" and
#     three rows of an integer 3x3 matrix in the lattice basis (acting on fractional coordinates, as spglib),
#     then "translation:" with one row of three fractional numbers, then "rotmap:" with lines "i -> j" (atom map).
#     The file has no time-reversal flag; it lists the operations of the crystal without the moments.
#   magmom_unprimed.log: header "part-I, no Time inversion (unprimed), delta(h)=1", then per operation
#     "----------------isymm iptr  N 0", the Cartesian rot_glb matrix, the moment checks per atom and "MagSymop ? True|False".
#   magmom_primed.log: same layout with "part-II, with Time inversion (primed), delta(h)=-1".
#   An operation is read as a magnetic symmetry without time reversal when the unprimed log says True for its index,
#   and with time reversal attached when the primed log says True; operation N in the logs is block N of mysymmop.dat.
SYMOP_HEAD = re.compile(r"\s*-{5,}\s+(\d+)\s+-{5,}")
LOG_HEAD = re.compile(r"-{5,}isymm iptr\s+(\d+)")


# Parse mysymmop.dat into operations with index, integer rotation and fractional translation.
def parse_mysymmop(path):
    ops, cur, mode, rows = [], None, None, []
    for line in iter_lines(path):
        m = SYMOP_HEAD.match(line)
        if m:
            cur = {"index": int(m.group(1)), "rotation": [], "translation": None}
            ops.append(cur)
            mode = None
            continue
        if cur is None:
            continue
        s = line.strip()
        if s.startswith("rotation:"):
            mode, rows = "rot", []
        elif s.startswith("translation:"):
            mode = "trans"
        elif s.startswith("rotmap:"):
            mode = None
        elif mode == "rot" and s:
            rows.append([int(round(float(x))) for x in s.split()[:3]])
            if len(rows) == 3:
                cur["rotation"], mode = rows, None
        elif mode == "trans" and s:
            cur["translation"] = [float(x) for x in s.split()[:3]]
            mode = None
    return [o for o in ops if len(o["rotation"]) == 3 and o["translation"] is not None]


# Indices marked "MagSymop ? True" in one magmom_*.log, or None when the file is missing.
def magsymop_true(path):
    if not Path(path).exists():
        return None
    out, cur = set(), None
    for line in iter_lines(path):
        m = LOG_HEAD.search(line)
        if m:
            cur = int(m.group(1))
        elif "MagSymop ?" in line and cur is not None and line.strip().endswith("True"):
            out.add(cur)
    return out


# Key of one operation for set comparison: rotation tuple, translation reduced mod 1 and rounded, time-reversal flag.
def op_key(rot, trans, tr=False):
    t = np.remainder(np.array(trans, dtype=float), 1.0)
    t[np.isclose(t, 1.0, atol=TRANS_TOL)] = 0.0
    return tuple(np.array(rot).astype(int).ravel()), tuple(np.round(t, 3)), bool(tr)


# wannhr_symm operations of a run with time-reversal flags from the logs, compared with spglib's groups.
def wannhr_symm_operations(run_dir, reference=None):
    run_dir = Path(run_dir)
    found = [run_dir / "mysymmop.dat"] + sorted(run_dir.glob("*/mysymmop.dat")) + sorted(run_dir.glob("*/*/mysymmop.dat"))
    folder = next((p.parent for p in found if p.exists()), None)
    if folder is None:
        return {"state": "not computed", "error": f"no mysymmop.dat in {run_dir} or two folder levels below it"}
    ops = parse_mysymmop(folder / "mysymmop.dat")
    unp, pri = magsymop_true(folder / "magmom_unprimed.log"), magsymop_true(folder / "magmom_primed.log")
    for o in ops:
        o["words"] = describe(o["rotation"], o["translation"])
        o["magnetic_unitary"] = None if unp is None else o["index"] in unp
        o["magnetic_antiunitary"] = None if pri is None else o["index"] in pri
    rel = folder.relative_to(run_dir.parent).as_posix()
    res = {"state": PENDING, "folder": rel, "n_operations": len(ops), "operations": ops,
           "n_magnetic_unitary": None if unp is None else sum(bool(o["magnetic_unitary"]) for o in ops),
           "n_magnetic_antiunitary": None if pri is None else sum(bool(o["magnetic_antiunitary"]) for o in ops),
           "logs": [n for n in ("magmom_unprimed.log", "magmom_primed.log") if (folder / n).exists()]}
    if reference is None:
        try:
            reference = magnetic_space_group(run_dir)
        except Exception as exc:
            res["comparison"] = f"not compared: no spglib result for {run_dir.name} ({type(exc).__name__}: {exc})"
            return res
    res["comparison"] = compare(res, reference)
    return res


# Compare wannhr_symm operations with spglib's plain and magnetic groups: "consistent" or "differs: ...".
def compare(w, ref):
    diffs = []
    sg = ref["space_group"]
    wset = {op_key(o["rotation"], o["translation"]) for o in w["operations"]}
    wset_t = {op_key(np.array(o["rotation"]).T, o["translation"]) for o in w["operations"]}
    sset = {op_key(o["rotation"], o["translation"]) for o in sg["operations"]}
    if w["n_operations"] != sg["n_operations"]:
        diffs.append(f'{w["n_operations"]} operations in mysymmop.dat, {sg["n_operations"]} in spglib\'s plain space group '
                     f'{sg["symbol"]} ({sg["number"]})')
    if wset != sset:
        diffs.append("operation set differs from spglib's plain space group" +
                     (" (it matches with transposed rotations)" if wset_t == sset else ""))
    if w.get("n_magnetic_unitary") is not None and w.get("n_magnetic_antiunitary") is not None:
        mag = {op_key(o["rotation"], o["translation"], False) for o in w["operations"] if o["magnetic_unitary"]} | \
              {op_key(o["rotation"], o["translation"], True) for o in w["operations"] if o["magnetic_antiunitary"]}
        rset = {op_key(o["rotation"], o["translation"], o["time_reversal"]) for o in ref["operations"]}
        if mag != rset:
            only_r = sorted({o["words"] + (" with time reversal" if o["time_reversal"] else "") for o in ref["operations"]
                             if op_key(o["rotation"], o["translation"], o["time_reversal"]) not in mag})
            only_w = sorted({o["words"] for o in w["operations"]
                             if (o["magnetic_unitary"] and op_key(o["rotation"], o["translation"], False) not in rset) or
                             (o["magnetic_antiunitary"] and op_key(o["rotation"], o["translation"], True) not in rset)})
            diffs.append(f'magnetic operations: wannhr_symm logs keep {w["n_magnetic_unitary"]} unitary and '
                         f'{w["n_magnetic_antiunitary"]} antiunitary, spglib keeps {ref["n_unitary"]} unitary and '
                         f'{ref["n_antiunitary"]} antiunitary ({ref["run"]} moments)' +
                         (f'; only in spglib: {", ".join(only_r)}' if only_r else "") +
                         (f'; only in wannhr_symm: {", ".join(only_w)}' if only_w else ""))
    return "consistent" if not diffs else "differs: " + "; ".join(diffs)


# Magnetic group, spin group and wannhr_symm results of one run for the index; each part fails alone, never the import.
def analyse(run_dir, has_outcar=True, has_symm=False):
    out = {}
    if has_outcar:
        try:
            out["magnetic"] = magnetic_space_group(run_dir)
        except Exception as exc:
            out["magnetic"] = {"state": "not computed", "error": f"{type(exc).__name__}: {exc}"}
        out["spin"] = spin_space_group(run_dir)
    if has_symm:
        try:
            ref = out.get("magnetic") if out.get("magnetic", {}).get("state") == PENDING else None
            out["wannhr"] = wannhr_symm_operations(run_dir, ref)
        except Exception as exc:
            out["wannhr"] = {"state": "not computed", "error": f"{type(exc).__name__}: {exc}"}
    return out
