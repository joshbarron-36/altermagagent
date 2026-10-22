"""Pure read-only parsers for FINDSYM logs, wannier90.win files and INCAR differences.

Nothing here writes inside a scanned folder; only write_findsym_input writes, and only to the path it is given.

FINDSYM log format (version 4.2, plain text, as written by `findsym < findsym.in > findsym.log`):

    FINDSYM, Version 4.2, June 2014              <- program name and version
    ...authors and university...
    CrSe                                         <- title line copied from findsym.in
    Tolerance:    0.00010
    Lattice parameters, a,b,c,alpha,beta,gamma:  <- the input cell, values on the next line
    Centering: P
    Number of atoms in unit cell:                <- count on the next line
    Type of each atom:                           <- for example "6*Cr 6*Se"
    Position and magnetic moment of each atom (dimensionless coordinates)
       1   x y z   mx my mz                      <- one line per input atom
    ------------------------------------------   <- separator between input echo and result
    Space Group 189.225 P-62'm'                  <- result: BNS number and BNS label for a magnetic search,
                                                    or "Space Group 194 D6h-4 P6_3/mmc" without moments
    Origin at    x y z
    Vectors a,b,c:                               <- three lines, the new basis
    Values of a,b,c,alpha,beta,gamma:            <- the conventional cell, values on the next line
    Atomic positions and magnetic moments in terms of a,b,c:
    Wyckoff position i, x = ..., mx = ...        <- one header per orbit, then its atoms

Newer FINDSYM versions write CIF tags instead, such as `_space_group_magn.number_BNS`, `_space_group_magn.name_BNS`,
`_space_group_magn.number_OG`, `_space_group_magn.name_OG`, `_space_group.IT_number` and a loop of
`_space_group_symop_magn_operation.xyz` lines; those tags are read too.
"""

import re
import warnings
from pathlib import Path

_FLOAT = r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?"
_SUBSCRIPTS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


# Turn a FINDSYM ASCII group label such as P-62'm' or P6_3/mmc into its Unicode form P6̄2′m′ or P6₃/mmc.
def unicode_label(label):
    if not label:
        return ""
    out = re.sub(r"_(\d)", lambda m: m.group(1).translate(_SUBSCRIPTS), label)
    out = re.sub(r"-(\d)", "\\1\u0305", out)
    return out.replace("'", "′")


# Read a list of floats from a text line, returning only the numeric tokens.
def _floats(line):
    vals = []
    for tok in line.split():
        try:
            vals.append(float(tok.replace("d", "e").replace("D", "e")))
        except ValueError:
            pass
    return vals


# Strip quotes FINDSYM CIF output puts around labels.
def _unquote(text):
    return text.strip().strip("\"")


# Look up the OG label, OG number and operation count for a BNS number in the pymatgen magnetic group table.
def _magnetic_table(bns_number):
    try:
        from pymatgen.symmetry.maggroups import MagneticSpaceGroup
        parts = [int(p) for p in bns_number.split(".")]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grp = MagneticSpaceGroup(parts)
        return {"og_label": grp._data["og_label"], "og_number": ".".join(str(n) for n in grp._data["og_number"]),
                "bns_label": grp._data["bns_label"], "num_operations": len(grp.symmetry_ops)}
    except Exception:
        return None


# Look up the Hermann-Mauguin symbol of an ordinary space group number with pymatgen.
def _space_group_symbol(number):
    try:
        from pymatgen.symmetry.groups import SpaceGroup
        return SpaceGroup.from_int_number(int(number)).symbol
    except Exception:
        return None


# Parse a FINDSYM log into space group, magnetic group labels, operations, lattice and notes; never raises.
def parse_findsym_log(path):
    try:
        text = Path(path).read_text(errors="replace")
    except Exception as exc:
        return {"error": f"cannot read {path}: {exc}"}
    try:
        return _parse_findsym_text(text)
    except Exception as exc:
        return {"error": f"cannot parse {path}: {exc}"}


# Parse the text of a FINDSYM log; the caller catches any exception.
def _parse_findsym_text(text):
    lines = text.splitlines()
    res = {"program": None, "title": None, "tolerance": None, "space_group_number": None,
           "space_group_symbol": None, "bns_number": None, "bns_label": None, "bns_label_unicode": None,
           "og_number": None, "og_label": None, "og_label_unicode": None, "magnetic": False,
           "num_operations": None, "lattice_input": None, "lattice": None, "centering": None,
           "num_atoms": None, "origin": None, "wyckoff": [], "notes": []}
    notes = res["notes"]
    ops = []
    in_ops = False
    for i, raw in enumerate(lines):
        line = raw.strip()
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        low = line.lower()
        if in_ops:
            if line and not line.startswith(("_", "loop_", "#")):
                ops.append(line)
                continue
            in_ops = False
        if line.startswith("FINDSYM") and res["program"] is None:
            res["program"] = line
            if i + 4 < len(lines):
                res["title"] = lines[i + 4].strip() or None
        elif low.startswith("tolerance:"):
            vals = _floats(line.split(":", 1)[1])
            res["tolerance"] = vals[0] if vals else None
        elif low.startswith("lattice parameters, a,b,c"):
            vals = _floats(nxt)
            res["lattice_input"] = vals[:6] if len(vals) >= 6 else None
        elif low.startswith("values of a,b,c,alpha,beta,gamma"):
            vals = _floats(nxt)
            res["lattice"] = vals[:6] if len(vals) >= 6 else None
        elif low.startswith("centering:"):
            res["centering"] = line.split(":", 1)[1].strip()
        elif low.startswith("number of atoms in unit cell"):
            vals = _floats(nxt)
            res["num_atoms"] = int(vals[0]) if vals else None
        elif low.startswith("origin at"):
            res["origin"] = _floats(line)[:3]
        elif low.startswith("wyckoff position"):
            m = re.match(r"wyckoff position\s+(\w+)", line, re.I)
            if m:
                res["wyckoff"].append(m.group(1))
        elif re.match(r"space group\s+\d", line, re.I):
            _read_space_group_line(line, res)
        elif re.search(r"number of (symmetry )?operations", low):
            vals = _floats(line.split(":", 1)[-1]) or _floats(nxt)
            if vals:
                res["num_operations"] = int(vals[0])
        elif line.startswith("_"):
            in_ops = _read_cif_tag(line, res)
        elif re.search(r"\bBNS\b", line) and res["bns_label"] is None:
            m = re.search(r"BNS\D*?(\d+\.\d+)\s+(\S+)", line)
            if m:
                res["bns_number"], res["bns_label"], res["magnetic"] = m.group(1), m.group(2), True
        if re.search(r"\bOG\b", line) and not line.startswith("_") and res["og_label"] is None:
            m = re.search(r"OG\D*?(\d+\.\d+\.\d+)\s+(\S+)", line)
            if m:
                res["og_number"], res["og_label"] = m.group(1), m.group(2)
    if ops and res["num_operations"] is None:
        res["num_operations"] = len(ops)
    _finish_findsym(res, notes)
    return res


# Read the "Space Group ..." result line of a FINDSYM 4.x log into the result dict.
def _read_space_group_line(line, res):
    toks = line.split()[2:]
    if not toks:
        return
    num = toks[0]
    if "." in num:
        res["bns_number"] = num
        res["bns_label"] = toks[1] if len(toks) > 1 else None
        res["magnetic"] = True
        res["space_group_number"] = int(num.split(".")[0])
    else:
        res["space_group_number"] = int(num)
        # Non-magnetic lines carry a Schoenflies symbol before the Hermann-Mauguin symbol.
        res["space_group_symbol"] = toks[-1] if len(toks) > 1 else None
        if len(toks) > 2:
            res["schoenflies"] = toks[1]


# Read one CIF tag line written by newer FINDSYM versions; returns True when an operation loop starts.
def _read_cif_tag(line, res):
    parts = line.split(None, 1)
    tag = parts[0].lower()
    val = _unquote(parts[1]) if len(parts) > 1 else ""
    if tag == "_space_group_magn.number_bns" and val:
        res["bns_number"], res["magnetic"] = val, True
    elif tag == "_space_group_magn.name_bns" and val:
        res["bns_label"], res["magnetic"] = val, True
    elif tag == "_space_group_magn.number_og" and val:
        res["og_number"] = val
    elif tag == "_space_group_magn.name_og" and val:
        res["og_label"] = val
    elif tag in ("_space_group.it_number", "_symmetry_int_tables_number") and val:
        res["space_group_number"] = int(float(val))
    elif tag in ("_space_group.name_h-m_alt", "_symmetry_space_group_name_h-m") and val:
        res["space_group_symbol"] = val
    elif tag in ("_space_group_symop_magn_operation.xyz", "_space_group_symop.operation_xyz",
                 "_symmetry_equiv_pos_as_xyz"):
        return True
    return False


# Fill derived fields (group number, symbol, Unicode labels, table lookups) and the notes list.
def _finish_findsym(res, notes):
    if res["bns_number"] and res["space_group_number"] is None:
        res["space_group_number"] = int(res["bns_number"].split(".")[0])
    if res["space_group_number"] and not res["space_group_symbol"]:
        sym = _space_group_symbol(res["space_group_number"])
        if sym:
            res["space_group_symbol"] = sym
            notes.append("space group symbol taken from the pymatgen table for the printed group number")
    res["space_group_symbol_unicode"] = unicode_label(res["space_group_symbol"])
    res["bns_label_unicode"] = unicode_label(res["bns_label"]) or None
    if res["bns_number"]:
        table = _magnetic_table(res["bns_number"])
        if table:
            res["num_operations_table"] = table["num_operations"]
            if table["bns_label"] != res["bns_label"]:
                notes.append(f"pymatgen table gives BNS label {table['bns_label']} for {res['bns_number']}")
            if res["og_label"] is None:
                res["og_label"], res["og_number"] = table["og_label"], table["og_number"]
                notes.append("OG label not printed by FINDSYM; taken from the pymatgen magnetic group table")
        if res["og_label"] is None:
            notes.append("no OG label printed")
    res["og_label_unicode"] = unicode_label(res["og_label"]) or None
    if res["space_group_number"] is None:
        notes.append("no space group line found; FINDSYM may have failed or the log is incomplete")
    if res["num_operations"] is None:
        notes.append("number of operations not listed in the log")
    if res["lattice"] is None and res["lattice_input"] is not None:
        notes.append("only the input lattice is listed")


# Write a findsym.in in the physicist's layout for a pymatgen-readable structure and one moment per atom.
def write_findsym_input(structure_path, magmom_vectors, out_path, title=None, tolerance=0.0001):
    from pymatgen.core import Structure
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = Structure.from_file(str(structure_path))
    moms = []
    for m in magmom_vectors:
        if isinstance(m, (int, float)):
            moms.append((0.0, 0.0, float(m)))
        else:
            vec = [float(x) for x in m]
            if len(vec) != 3:
                raise ValueError(f"magnetic moment {m!r} is not a 3-vector")
            moms.append(tuple(vec))
    if len(moms) != len(st):
        raise ValueError(f"{len(moms)} moments given for {len(st)} atoms")
    lat = st.lattice
    out = [title or st.composition.reduced_formula, f"{tolerance:g}", "2",
           f"{lat.a:.5f} {lat.b:.5f} {lat.c:.5f}  {lat.alpha:.4f}  {lat.beta:.4f} {lat.gamma:.4f}",
           "2", "P", str(len(st)), " ".join(site.specie.symbol for site in st), "magnetic"]
    for site, mom in zip(st, moms):
        x, y, z = site.frac_coords
        out.append(f"{x:.6f} {y:.6f} {z:.6f}     {mom[0]:.4f} {mom[1]:.4f} {mom[2]:.4f}")
    Path(out_path).write_text("\n".join(out) + "\n")
    return str(out_path)


# Parse a wannier90.win file into counts, windows, projections, flags and a one-line summary; never raises.
def parse_wannier_win(path):
    try:
        text = Path(path).read_text(errors="replace")
    except Exception as exc:
        return {"error": f"cannot read {path}: {exc}"}
    try:
        return parse_wannier_win_text(text)
    except Exception as exc:
        return {"error": f"cannot parse {path}: {exc}"}


# Turn a wannier90 logical value (true, .true., T, F) into a Python bool, or None when it is not logical.
def _win_bool(val):
    v = val.strip().strip(".").lower()
    if v in ("true", "t"):
        return True
    if v in ("false", "f"):
        return False
    return None


# Parse the text of a wannier90.win file; comments start with "#" or "!" and keys are case-insensitive.
def parse_wannier_win_text(text):
    keys, blocks, current = {}, {}, None
    for raw in text.splitlines():
        line = re.split(r"[#!]", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        low = line.lower()
        m = re.match(r"begin\s+(\w+)", low)
        if m:
            current = m.group(1)
            blocks[current] = []
            continue
        if re.match(r"end\s+\w+", low):
            current = None
            continue
        if current is not None:
            blocks[current].append(line)
            continue
        m = re.match(r"([A-Za-z_][\w]*)\s*(?:[=:]\s*|\s+)(.*)$", line)
        if m:
            keys[m.group(1).lower()] = m.group(2).strip()
    res = {}
    for k in ("num_wann", "num_bands"):
        vals = _floats(keys.get(k, ""))
        res[k] = int(vals[0]) if vals else None
    for k in ("dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max"):
        vals = _floats(keys.get(k, ""))
        res[k] = vals[0] if vals else None
    res["projections"] = blocks.get("projections", [])
    for k in ("spinors", "bands_plot", "write_hr"):
        res[k] = _win_bool(keys[k]) if k in keys else False
    res["kpoint_path"] = bool(blocks.get("kpoint_path"))
    res["fermi_energy"] = (_floats(keys["fermi_energy"]) or [None])[0] if "fermi_energy" in keys else None
    res["summary"] = _win_summary(res)
    return res


# Format a number for a summary line without trailing zeros (10.5, -1, 8.337).
def _num(x):
    return f"{x:.3f}".rstrip("0").rstrip(".") if isinstance(x, float) else str(x)


# Build the one-line summary of a parsed wannier90.win.
def _win_summary(res):
    parts = []
    if res["num_wann"] is not None:
        parts.append(f"num_wann {res['num_wann']}")
    if res["num_bands"] is not None:
        parts.append(f"bands {res['num_bands']}")
    for name, lo, hi in (("outer window", "dis_win_min", "dis_win_max"), ("frozen", "dis_froz_min", "dis_froz_max")):
        a, b = res[lo], res[hi]
        if a is None and b is None:
            parts.append(f"{name}: no window set")
        else:
            parts.append(f"{name} {_num(a) if a is not None else 'default'} to {_num(b) if b is not None else 'default'} eV")
    if res["projections"]:
        parts.append("projections: " + ", ".join(re.sub(r"\s+", "", p) for p in res["projections"]))
    else:
        parts.append("no projections block")
    flags = [k for k in ("spinors", "write_hr", "bands_plot") if res[k]]
    if res["kpoint_path"]:
        flags.append("kpoint_path")
    if flags:
        parts.append("flags: " + " ".join(flags))
    return ", ".join(parts)


# Read an INCAR with pymatgen and return a plain dict of tag to value with MAGMOM flattened.
def _read_incar(path):
    from pymatgen.io.vasp import Incar
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inc = Incar.from_file(str(path))
    out = {}
    for k, v in inc.items():
        out[k.upper()] = _flatten(v) if isinstance(v, list) else v
    return out


# Flatten a nested list such as a non-collinear MAGMOM into one list.
def _flatten(v):
    flat = []
    for x in v:
        if isinstance(x, (list, tuple)):
            flat.extend(_flatten(x))
        else:
            flat.append(x)
    return flat


# Normalise one INCAR value so that .TRUE., T and True compare equal and numbers compare as numbers.
def _norm(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        b = _win_bool(v)
        if b is not None:
            return b
        try:
            return float(v)
        except ValueError:
            return v.strip().lower()
    return v


# Decide whether two INCAR values are equal, lists element-wise with a 0.01 tolerance.
def _same(a, b):
    a, b = _norm(a), _norm(b)
    if isinstance(a, list) or isinstance(b, list):
        if not (isinstance(a, list) and isinstance(b, list)) or len(a) != len(b):
            return False
        return all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= 0.01 if isinstance(a, float) or isinstance(b, float) else a == b
    return a == b


# Show an INCAR value in VASP style for a summary line.
def _show(v):
    if isinstance(v, str) and _win_bool(v) is not None:
        v = _win_bool(v)
    if isinstance(v, bool):
        return ".TRUE." if v else ".FALSE."
    if isinstance(v, float):
        return f"{v:g}"
    if isinstance(v, list):
        return " ".join(_show(x) for x in v)
    return str(v)


# Tags whose values are long lists and are summarised by name only.
_LONG_TAGS = {"MAGMOM", "M_CONSTR", "RWIGS", "LDAUU", "LDAUJ", "LDAUL"}


# Compare an INCAR with its parent INCAR and return changed, added and removed tags with a one-line summary.
def incar_diff(incar_path, parent_incar_path):
    try:
        this, parent = _read_incar(incar_path), _read_incar(parent_incar_path)
    except Exception as exc:
        return {"error": f"cannot read INCAR: {exc}"}
    return diff_incar_dicts(this, parent)


# Compare two INCAR tag dicts (this run and its parent) and build the diff result.
def diff_incar_dicts(this, parent):
    this = {k.upper(): (_flatten(v) if isinstance(v, list) else v) for k, v in this.items()}
    parent = {k.upper(): (_flatten(v) if isinstance(v, list) else v) for k, v in parent.items()}
    changed, added, removed = {}, {}, {}
    for k in sorted(set(this) | set(parent)):
        if k in this and k not in parent:
            added[k] = this[k]
        elif k in parent and k not in this:
            removed[k] = parent[k]
        elif not _same(this[k], parent[k]):
            changed[k] = (parent[k], this[k])
    items = []
    for k, (old, new) in changed.items():
        items.append(f"{k} changed" if k in _LONG_TAGS or k == "SYSTEM" else f"{k} {_show(old)} -> {_show(new)}")
    for k, v in added.items():
        items.append(f"{k} added" if k in _LONG_TAGS else f"{k} added {_show(v)}")
    for k in removed:
        items.append(f"{k} removed")
    # Physics tags first, the SYSTEM title last, so the six shown are the informative ones.
    items.sort(key=lambda s: s.startswith("SYSTEM"))
    summary = ", ".join(items[:6]) if items else "no changes"
    if len(items) > 6:
        summary += f", and {len(items) - 6} more"
    return {"changed": changed, "added": added, "removed": removed, "summary": summary}
