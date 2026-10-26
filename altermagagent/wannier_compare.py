"""Compare an indexed Wannier Hamiltonian with bands from its recorded calculation lineage."""

import csv
import hashlib
import json
import re
from collections import namedtuple
from pathlib import Path
from uuid import uuid4

import numpy as np

from . import bands, db
from .files import display_path

STRUCTURE_TOL = 1e-6
KPOINT_TOL = 1e-7
HERMITIAN_TOL_EV = 1e-5
WS_FALLBACK_NOTE = ("no seedname_wsvec.dat and no use_ws_distance=false declaration: "
                    "direct-HR convention used as a fallback, not the Wannier90 default")


# Hash an input in bounded memory and reject a file changed during the read.
def file_hash(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Input changed while being checked: {path}")
    return digest.hexdigest()


# Require the same lattice basis, ordered species and periodic coordinates without rescaling.
def require_same_structure(first, second):
    if list(first.species) != list(second.species):
        raise ValueError("POSCAR species or atom order differ.")
    if not np.allclose(first.lattice.matrix, second.lattice.matrix, atol=STRUCTURE_TOL, rtol=0):
        raise ValueError("POSCAR lattice bases differ; the same fractional k-path cannot be reused.")
    delta = np.asarray(first.frac_coords) - np.asarray(second.frac_coords)
    delta -= np.rint(delta)
    if np.max(np.linalg.norm(delta @ first.lattice.matrix, axis=1)) > STRUCTURE_TOL:
        raise ValueError("POSCAR atomic positions differ.")


# Require identical SOC, noncollinearity, spin axis and starting magnetic texture even for preliminary plots.
def require_same_magnetism(left, right):
    from pymatgen.io.vasp import Incar

    first, second = (Incar.from_file(Path(folder) / "INCAR") for folder in (left, right))
    for key in ("LSORBIT", "LNONCOLLINEAR"):
        if bool(first.get(key, False)) != bool(second.get(key, False)):
            raise ValueError(f"{key} differs; reference and Hamiltonian must have the same magnetic treatment.")
    for key, default in (("SAXIS", [0, 0, 1]), ("MAGMOM", [])):
        a, b = np.asarray(first.get(key, default)), np.asarray(second.get(key, default))
        if a.shape != b.shape or not np.allclose(a, b, atol=1e-6, rtol=0):
            raise ValueError(f"{key} differs; the same magnetic configuration is required.")
    return first


# Check structure and required input equality, retaining a complete audit for preliminary overlays.
def check_inputs(wannier_run, reference_run, structure_only=False):
    from pymatgen.io.vasp import Incar, Poscar

    left, right = Path(wannier_run["path"]), Path(reference_run["path"])
    structures = [Poscar.from_file(p / "POSCAR").structure for p in (left, right)]
    require_same_structure(*structures)
    incar = require_same_magnetism(left, right)
    checks = {"structure": {"status": "matched", "tolerance_angstrom": STRUCTURE_TOL},
              "magnetic_treatment": {"soc": bool(incar.get("LSORBIT", False)),
                                     "noncollinear": bool(incar.get("LNONCOLLINEAR", False))}}
    for name in ("POSCAR", "INCAR", "CHGCAR", "CHG"):
        paths = [p / name for p in (left, right)]
        hashes = [file_hash(p) if p.is_file() else None for p in paths]
        same = all(hashes) and hashes[0] == hashes[1]
        checks[name] = {"status": "matched" if same else "missing" if None in hashes else "different",
                        "paths": list(map(str, paths)), "sha256": hashes}
    if all((p / "INCAR").is_file() for p in (left, right)):
        inc = [dict(Incar.from_file(p / "INCAR")) for p in (left, right)]
        differences = {k: [inc[0].get(k), inc[1].get(k)] for k in sorted(inc[0].keys() | inc[1].keys())
                       if inc[0].get(k) != inc[1].get(k)}
        checks["incar_differences"] = differences
    failed = [name for name in ("POSCAR", "INCAR", "CHGCAR", "CHG") if checks[name]["status"] != "matched"]
    if failed and not structure_only:
        detail = ", ".join(f"{name}: {checks[name]['status']}" for name in failed)
        raise ValueError(f"Required input equality failed ({detail}). "
                         "Use --structure-only only for an explicitly preliminary overlay.")
    return structures[0], checks, failed


# Follow recorded parents with cycle, missing-parent and ambiguity checks.
def ancestry(run, by_name):
    chain, seen = [], {run["name"]}
    current = run
    while current.get("parent"):
        name = current["parent"]
        if name in seen:
            raise ValueError(f"Parent cycle at {name}; correct the indexed lineage.")
        if name not in by_name:
            raise ValueError(f"Indexed parent {name} is missing; re-import or correct the lineage.")
        seen.add(name)
        current = by_name[name]
        chain.append(current)
    return chain


# No band run qualifies as the reference of a Wannier run.
class NoReference(ValueError):
    """Raised when no same-structure band run can serve as the comparison reference."""


# Identify band runs without guessing from directory names.
def is_band(run):
    details = run.get("details") or {}
    return details.get("phase") == "band" or details.get("kind") == "bands"


# Select a literal band ancestor or an unambiguous earlier band run sharing an SCF ancestor; optionally the latest one.
# Returns the reference, the lineage, the Wannier run's ancestor chain and its parsed POSCAR structure.
def select_reference(wannier_run, runs, reference_run=None, latest_if_ambiguous=False):
    from pymatgen.io.vasp import Poscar

    by_name = {r["name"]: r for r in runs}
    chain = ancestry(wannier_run, by_name)
    direct = [r for r in chain if is_band(r)]
    candidates = list(direct)
    scf_names = {r["name"] for r in chain if r.get("details", {}).get("phase") == "scf"}
    if not direct or reference_run:
        for run in runs:
            if not is_band(run) or run in candidates or run["seq"] >= wannier_run["seq"]:
                continue
            if scf_names.intersection(r["name"] for r in ancestry(run, by_name)):
                candidates.append(run)
    structure = Poscar.from_file(Path(wannier_run["path"]) / "POSCAR").structure
    compatible = []
    for run in candidates:
        try:
            require_same_structure(structure, Poscar.from_file(Path(run["path"]) / "POSCAR").structure)
            require_same_magnetism(wannier_run["path"], run["path"])
        except (OSError, ValueError):
            continue
        compatible.append(run)
    note = None
    if reference_run:
        reference = next((r for r in compatible if r["name"] == reference_run), None)
        if reference is None:
            raise ValueError(f"{reference_run} is not a same-structure band ancestor or earlier shared-SCF band run.")
    elif direct and direct[0] in compatible:
        reference = direct[0]
    elif len(compatible) == 1:
        reference = compatible[0]
    elif compatible and latest_if_ambiguous:
        reference = max(compatible, key=lambda r: r["seq"])
        note = f"several eligible band runs; the latest, {reference['name']}, was used"
    else:
        names = ", ".join(r["name"] for r in compatible) or "none"
        raise NoReference(f"No unique same-structure reference (eligible runs: {names}). "
                          "Use --reference-run to resolve candidates; correct unknown parents with altermag parent.")
    ref_chain = ancestry(reference, by_name)
    common = next((r for r in chain if r["name"] in {x["name"] for x in ref_chain}), None)
    return reference, {"relation": "band ancestor" if reference in direct else "shared SCF ancestor",
                       "wannier_ancestors": [r["name"] for r in chain],
                       "reference_ancestors": [r["name"] for r in ref_chain],
                       "common_ancestor": common["name"] if common else None,
                       "parent_notes": {r["name"]: r.get("parent_note") for r in [wannier_run, *chain, reference, *ref_chain]},
                       "selection_note": note}, chain, structure


# Find the requested Hamiltonian and its most specific indexed owner without a latest-run guess.
def locate_hamiltonian(project, runs, hr_file=None, run_name=None):
    root = Path(project["root"]).resolve()
    if hr_file:
        path = Path(hr_file).expanduser()
        path = (path if path.is_absolute() else root / path).resolve()
        if not bands.inside(path, root):
            raise ValueError("The Hamiltonian must belong to this indexed project's root.")
        owners = [r for r in runs if bands.inside(path, r["path"])]
        if not owners:
            raise ValueError("No indexed run owns this Hamiltonian; import its calculation first.")
        run = max(owners, key=lambda r: len(Path(r["path"]).parts))
        if run_name and run["name"] != run_name:
            raise ValueError("The Hamiltonian does not belong to the requested run.")
    else:
        candidates = [r for r in runs if (not run_name or r["name"] == run_name)
                      and (Path(r["path"]) / "wannier90_hr.dat").is_file()]
        if len(candidates) != 1:
            raise ValueError("Specify hr_file or run_name; no unique indexed wannier90_hr.dat was found.")
        run = candidates[0]
        path = (Path(run["path"]) / "wannier90_hr.dat").resolve()
        if not bands.inside(path, root):
            raise ValueError("The Hamiltonian resolves outside the indexed project.")
    if not path.is_file() or "_hr" not in path.name or ".dat" not in path.name:
        raise ValueError(f"Expected an existing seedname_hr.dat file or a variant such as seedname_hr.dat_*: {path}")
    return path, run


# Parse real-space Hamiltonian blocks with explicit indices and Wigner-Seitz degeneracies.
def read_hr(path):
    with Path(path).open() as stream:
        header = stream.readline().strip()
        nw, nr = int(stream.readline()), int(stream.readline())
        if nw < 1 or nr < 1:
            raise ValueError("Hamiltonian dimensions must be positive.")
        degeneracies = []
        while len(degeneracies) < nr:
            line = stream.readline()
            if not line:
                raise ValueError("Truncated Hamiltonian degeneracy block.")
            degeneracies.extend(map(int, line.split()))
        if len(degeneracies) != nr or min(degeneracies) <= 0:
            raise ValueError("Invalid Hamiltonian degeneracies.")
        vectors, matrices = [], []
        for _ in range(nr):
            matrix, seen, vector = np.zeros((nw, nw), complex), set(), None
            for _ in range(nw * nw):
                row = stream.readline().split()
                if len(row) != 7:
                    raise ValueError("Truncated or malformed Hamiltonian matrix row.")
                r = tuple(map(int, row[:3]))
                m, n = int(row[3]) - 1, int(row[4]) - 1
                if not (0 <= m < nw and 0 <= n < nw) or (m, n) in seen or (vector is not None and r != vector):
                    raise ValueError("Inconsistent Hamiltonian indices or R blocks.")
                vector = r
                seen.add((m, n))
                matrix[m, n] = complex(*(float(x.replace("D", "E").replace("d", "e")) for x in row[5:]))
            vectors.append(vector)
            matrices.append(matrix)
        if any(line.strip() for line in stream) or len(set(vectors)) != nr:
            raise ValueError("Unexpected extra rows or repeated R blocks in Hamiltonian.")
    matrices = np.asarray(matrices)
    if not np.isfinite(matrices).all():
        raise ValueError("Non-finite Hamiltonian matrix elements.")
    return {"header": header, "vectors": np.asarray(vectors), "matrices": matrices,
            "degeneracies": np.asarray(degeneracies), "num_wann": nw}


# Read small Wannier control flags, including flags written in export headers.
def control_flags(text):
    return {key.lower(): value.lower().strip(".") in ("t", "true")
            for key, value in re.findall(r"(use_ws_distance|write_ndegen_applied)\s*[=:]\s*(\.?\w+\.?)", text, re.I)}


# The seedname.win and seedname_wsvec.dat paths next to a seedname_hr.dat file.
def companion_files(path):
    seed = path.name.split("_hr", 1)[0]
    return path.with_name(seed + ".win"), path.with_name(seed + "_wsvec.dat")


# Fold pair-dependent superlattice translations once into effective real-space matrices.
def interpolation_terms(path, hr, use_ws_distance=None):
    from .importer import read_win

    win, wsvec = companion_files(path)
    flags = {}
    if win.is_file():
        settings = read_win(win)
        flags.update(control_flags("\n".join(f"{k}={v}" for k, v in settings.items())))
        if "num_wann" in settings and int(settings["num_wann"]) != hr["num_wann"]:
            raise ValueError("num_wann in the .win file disagrees with the Hamiltonian.")
    flags.update(control_flags(hr["header"]))
    wsheader = ""
    if wsvec.is_file():
        with wsvec.open() as stream:
            wsheader = stream.readline()
        header_flags = control_flags(wsheader)
        for key, val in header_flags.items():
            if key in flags and flags[key] != val:
                raise ValueError(f"Contradictory {key} in Wannier inputs and export headers.")
        flags.update(header_flags)
    applied = flags.get("write_ndegen_applied", False)
    if applied and np.any(hr["degeneracies"] != 1):
        raise ValueError("An export with applied degeneracies must have unit degeneracy entries.")
    declared_false = flags.get("use_ws_distance") is False
    if use_ws_distance is not None:
        flags["use_ws_distance"] = bool(use_ws_distance)
    use_ws = flags.get("use_ws_distance", True) and not applied
    if use_ws and not wsvec.is_file():
        raise ValueError("Missing seedname_wsvec.dat for Wigner-Seitz interpolation. "
                         "Provide it, or a .win/export header explicitly declaring use_ws_distance=false.")
    if use_ws:
        nw = hr["num_wann"]
        lookup = {tuple(r): i for i, r in enumerate(hr["vectors"])}
        seen = np.zeros(hr["matrices"].shape, dtype=bool)
        expanded = {}
        with wsvec.open() as stream:
            next(stream)
            lines = iter(line for line in stream if line.strip())
            try:
                for line in lines:
                    key = tuple(map(int, line.split()))
                    if len(key) != 5 or key[:3] not in lookup:
                        raise ValueError("Malformed wsvec pair or unknown R vector.")
                    r, m, n = key[:3], key[3] - 1, key[4] - 1
                    ir = lookup[r]
                    if not (0 <= m < nw and 0 <= n < nw) or seen[ir, m, n]:
                        raise ValueError("Invalid or repeated wsvec orbital indices.")
                    seen[ir, m, n] = True
                    count = int(next(lines))
                    if count < 1:
                        raise ValueError("Invalid wsvec translation count.")
                    value = hr["matrices"][ir, m, n] / (hr["degeneracies"][ir] * count)
                    for _ in range(count):
                        translation = tuple(map(int, next(lines).split()))
                        if len(translation) != 3:
                            raise ValueError("Malformed wsvec translation.")
                        vector = tuple(a + b for a, b in zip(r, translation))
                        if vector not in expanded:
                            expanded[vector] = np.zeros((nw, nw), complex)
                        expanded[vector][m, n] += value
            except StopIteration:
                raise ValueError("Truncated wsvec file.") from None
        if not seen.all():
            raise ValueError("wsvec is missing Hamiltonian pairs.")
        hr = dict(hr, vectors=np.asarray(list(expanded)), matrices=np.asarray(list(expanded.values())),
                  degeneracies=np.ones(len(expanded)))
    return hr, {"use_ws_distance": use_ws, "write_ndegen_applied": applied,
                "win": str(win) if win.is_file() else None, "wsvec": str(wsvec) if use_ws else None,
                "effective_R_vectors": len(hr["vectors"]),
                "convention_note": None if use_ws or applied or declared_false else WS_FALLBACK_NOTE}


# Call visit(index, Hermitian part of H(k)) at each fractional k-point; return the largest Hermiticity residual.
def hamiltonians(hr, kfrac, visit):
    kfrac = np.asarray(kfrac, dtype=float)
    if kfrac.ndim != 2 or kfrac.shape[1] != 3 or not np.isfinite(kfrac).all():
        raise ValueError("Expected finite fractional k-points with shape (n, 3).")
    nw = hr["num_wann"]
    weighted = hr["matrices"] / hr["degeneracies"][:, None, None]
    max_residual = 0.0
    for start in range(0, len(kfrac), 32):
        phases = np.exp(2j * np.pi * (kfrac[start:start + 32] @ hr["vectors"].T))
        block = (phases @ weighted.reshape(len(weighted), -1)).reshape(-1, nw, nw)
        for local, hk in enumerate(block):
            residual = float(np.max(np.abs(hk - hk.conj().T)))
            if residual > HERMITIAN_TOL_EV:
                raise ValueError(f"Non-Hermitian H(k) at k-point {start + local + 1}: residual {residual:g} eV.")
            max_residual = max(max_residual, residual)
            visit(start + local, (hk + hk.conj().T) / 2)
    return max_residual


# Fourier transform and diagonalize at the exact reference k-points, rejecting non-Hermitian results.
def interpolate(hr, kfrac):
    result = np.empty((len(kfrac), hr["num_wann"]))

    # Store the eigenvalues of one H(k) in its row.
    def store(ik, hk):
        result[ik] = np.linalg.eigvalsh(hk)

    residual = hamiltonians(hr, kfrac, store)
    return result.T, residual


# Expand reciprocal line-mode KPOINTS and preserve every branch boundary and duplicate endpoint.
def reference_path(folder, structure):
    from pymatgen.io.vasp import Kpoints

    kp = Kpoints.from_file(folder / "KPOINTS")
    if "line" not in str(kp.style).lower() or not str(kp.coord_type).lower().startswith("r"):
        raise ValueError("Comparison requires reciprocal-coordinate line-mode KPOINTS.")
    count = kp.num_kpts
    if count < 2 or len(kp.kpts) % 2 or not kp.kpts:
        raise ValueError("Malformed line-mode KPOINTS endpoints or subdivision count.")
    frac, distance, branches, ticks = [], [], [], []
    offset = 0.0
    labels = kp.labels or [""] * len(kp.kpts)
    for j in range(0, len(kp.kpts), 2):
        start, stop = np.asarray(kp.kpts[j]), np.asarray(kp.kpts[j + 1])
        segment = np.linspace(start, stop, count)
        length = np.linalg.norm((stop - start) @ structure.lattice.reciprocal_lattice.matrix)
        a = len(frac)
        frac.extend(segment)
        distance.extend(np.linspace(offset, offset + length, count))
        branches.append((a, a + count - 1))
        ticks.extend([(offset, bands.pretty_label(labels[j])), (offset + length, bands.pretty_label(labels[j + 1]))])
        offset += length
    if offset <= 0:
        raise ValueError("The reference k-path has zero length.")
    return {"kfrac": np.asarray(frac), "x": np.asarray(distance), "branches": branches,
            "ticks": bands.build_ticks(ticks, offset)}


# True for a spin-polarised collinear INCAR mapping: ISPIN 2 without LSORBIT and without LNONCOLLINEAR.
def collinear_polarised(incar):
    return incar.get("ISPIN", 1) == 2 and not (incar.get("LSORBIT") or incar.get("LNONCOLLINEAR"))


# Read absolute VASP eigenvalues or undo the recorded VASPKIT energy shift.
def load_reference(folder, structure, spin=None):
    from pymatgen.electronic_structure.core import Spin
    from pymatgen.io.vasp.outputs import Eigenval, Vasprun

    data = reference_path(folder, structure)
    fermi = None
    if (folder / "EIGENVAL").is_file():
        run = Eigenval(folder / "EIGENVAL")
        kfrac, values = run.kpoints, run.eigenvalues
        source = folder / "EIGENVAL"
    elif (folder / "vasprun.xml").is_file():
        run = Vasprun(folder / "vasprun.xml", parse_projected_eigen=False, parse_potcar_file=False)
        require_same_structure(structure, run.initial_structure)
        kfrac, values = run.actual_kpoints, run.eigenvalues
        fermi = run.efermi
        source = folder / "vasprun.xml"
    else:
        source = next((folder / name for name in ("REFORMATTED_BAND.dat", "REFORMATED_BAND.dat")
                       if (folder / name).is_file()), None)
        if source is None:
            raise ValueError("The ancestor needs vasprun.xml, EIGENVAL or REFORMATTED_BAND.dat; a figure alone is insufficient.")
        fermi = read_fermi(folder / "FERMI_ENERGY")
        table = np.loadtxt(source, comments="#", ndmin=2)
        if not np.isfinite(table).all():
            raise ValueError("Non-finite VASPKIT energies or path distances.")
        if table.shape[1] < 2 or len(table) != len(data["x"]):
            raise ValueError("VASPKIT bands do not have the expected line-mode k-point count.")
        if table[-1, 0] <= table[0, 0] or not np.allclose(
                (table[:, 0] - table[0, 0]) / (table[-1, 0] - table[0, 0]),
                data["x"] / data["x"][-1], atol=1e-5, rtol=0):
            raise ValueError("VASPKIT distances disagree with the ancestor KPOINTS path.")
        from pymatgen.io.vasp import Incar
        if collinear_polarised(Incar.from_file(folder / "INCAR")):
            raise ValueError("Spin-polarised VASPKIT tables need EIGENVAL or vasprun.xml to identify spin columns.")
        if spin == "down":
            raise ValueError("The reference table contains a single channel; spin down is unavailable.")
        data.update(energies=table[:, 1:].T + fermi, source=str(source), fermi=fermi,
                    kpath_evidence="KPOINTS expansion checked against VASPKIT distances; no output k-coordinates")
        return data
    actual = np.asarray(kfrac)
    if actual.shape != data["kfrac"].shape or not np.allclose(actual, data["kfrac"], atol=KPOINT_TOL, rtol=0):
        raise ValueError("Output k-points disagree with line-mode KPOINTS; refusing to overlay different paths.")
    if len(values) > 1 and spin is None:
        raise ValueError("The reference has two spin channels; specify --spin up or --spin down for this Hamiltonian.")
    selected = Spin.down if spin == "down" else Spin.up
    if selected not in values:
        raise ValueError("The requested spin channel is absent from the reference output.")
    energies = np.asarray(values[selected])[:, :, 0].T
    if not np.isfinite(energies).all():
        raise ValueError("Non-finite reference energies.")
    data.update(energies=energies, source=str(source), fermi=fermi, kpath_evidence="output coordinates checked against KPOINTS")
    return data


# Read a scalar Fermi-energy file without assuming a project-specific value.
def read_fermi(path):
    for line in path.read_text().splitlines():
        text = line.split("#", 1)[0].strip()
        if text:
            value = float(text.split()[0])
            if np.isfinite(value):
                return value
    raise ValueError(f"No finite Fermi energy in {path}")


# Prefer a common SCF energy reference and otherwise explicitly use the band calculation's value.
def common_fermi(chain, lineage, data, reference, explicit=None):
    if explicit is not None:
        if not np.isfinite(float(explicit)):
            raise ValueError("Fermi energy must be finite.")
        return float(explicit), "explicit --fermi-eV"
    shared = set(lineage["reference_ancestors"])
    for run in chain:
        if run["name"] not in shared or run.get("details", {}).get("phase") != "scf":
            continue
        folder = Path(run["path"])
        if (folder / "OUTCAR").is_file():
            value = None
            with (folder / "OUTCAR").open(errors="replace") as stream:
                for line in stream:
                    match = re.search(r"E-fermi\s*:\s*([-+\d.eEdD]+)", line)
                    if match:
                        value = float(match[1].replace("D", "E"))
            if value is not None and np.isfinite(value):
                return value, str(folder / "OUTCAR")
        if (folder / "FERMI_ENERGY").is_file():
            return read_fermi(folder / "FERMI_ENERGY"), str(folder / "FERMI_ENERGY")
    value = data["fermi"]
    if value is None and (Path(reference["path"]) / "FERMI_ENERGY").is_file():
        return read_fermi(Path(reference["path"]) / "FERMI_ENERGY"), str(Path(reference["path"]) / "FERMI_ENERGY")
    if value is None or not np.isfinite(value):
        raise ValueError("No usable Fermi energy found; supply --fermi-eV from the common SCF calculation.")
    return float(value), data["source"] + " (band calculation; common SCF Fermi energy unavailable)"


# Compute bandwise errors only for an explicitly supplied contiguous band correspondence.
def error_metrics(reference, wannier, band_start, emin, emax):
    if band_start is None:
        return None
    if isinstance(band_start, bool) or int(band_start) != band_start:
        raise ValueError("band_start must be a one-based integer.")
    start = int(band_start) - 1
    if start < 0 or start + len(wannier) > len(reference):
        raise ValueError("The supplied reference band range does not fit the available bands.")
    selected = reference[start:start + len(wannier)]
    mask = (selected >= emin) & (selected <= emax)
    errors = (wannier - selected)[mask]
    if not errors.size:
        raise ValueError("No mapped reference energies fall inside the requested window.")
    absolute = np.abs(errors) * 1000
    return {"mapping": f"Wannier bands 1-{len(wannier)} to reference bands {start + 1}-{start + len(wannier)}",
            "samples": int(errors.size), "mae_meV": float(absolute.mean()),
            "rmse_meV": float(np.sqrt(np.mean(absolute ** 2))), "p95_meV": float(np.quantile(absolute, 0.95)),
            "max_meV": float(absolute.max()), "samples_over_20_meV": int(np.count_nonzero(absolute > 20))}


# Text of the error box: MAE, RMSE, 95th percentile, max and count above 20 meV, the window, the mapping and its source.
def error_box_text(result):
    m = result["metrics"]
    if m is None:
        return "no band mapping:\nno errors computed"
    return (f"MAE {m['mae_meV']:.3f} meV\nRMSE {m['rmse_meV']:.3f} meV\n95th percentile {m['p95_meV']:.3f} meV\n"
            f"max {m['max_meV']:.3f} meV\n{m['samples_over_20_meV']} of {m['samples']} samples above 20 meV\n"
            f"window {result['metric_emin']:g} to {result['metric_emax']:g} eV\n{m['mapping']}\n{m['mapping_source']}")


# Draw new static overlays while keeping disconnected path segments separate: full window left, metric window right.
def draw_overlay(data, wannier, result, paths):
    from matplotlib.figure import Figure

    figure = Figure(figsize=(11, 5))
    axes = figure.subplots(1, 2, sharex=True)
    limits = ((result["emin"], result["emax"]), (result["metric_emin"], result["metric_emax"]))
    for axis, ylim in zip(axes, limits):
        for energies, colour, style, label in ((data["energies"], "#202020", "-", "VASP"),
                                                (wannier, "#d95f02", "--", "Wannier Hamiltonian")):
            first = True
            for band in energies:
                for xs, ys in bands.segments(data["x"], band, data["branches"]):
                    axis.plot(xs, ys, color=colour, ls=style, lw=0.8, label=label if first else None)
                    first = False
        bands.decorate_band_axes(axis, data, ylim, zero_colour="#2166ac")
        axis.set_xlabel("Wave vector")
    axes[0].set_title(f"{result['wannier_run']} versus {result['reference_run']}", fontsize=9)
    axes[1].set_title(f"error window {result['metric_emin']:g} to {result['metric_emax']:g} eV", fontsize=9)
    axes[0].legend(loc="lower right", fontsize=8)
    axes[1].text(0.02, 0.98, error_box_text(result), transform=axes[1].transAxes, va="top", ha="left", fontsize=8,
                 family="monospace", bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.5", "alpha": 0.9})
    footer = "\n".join(filter(None, [result["evidence_state"], result["interpolation"].get("convention_note")]))
    figure.text(0.5, 0.01, footer, ha="center", va="bottom", fontsize=9)
    figure.tight_layout(rect=(0, 0.07, 1, 1))
    for kind in ("png", "pdf", "svg"):
        figure.savefig(paths[kind], dpi=250)


# First reference band of the Wannier manifold from a leading exclude_bands block (1-12 gives 13), or None.
def excluded_offset(hr_path):
    win = hr_path.parent / "wannier90.win"
    if not win.is_file():
        win = next((p / "wannier90.win" for p in hr_path.parents if (p / "wannier90.win").is_file()), None)
    if win is None:
        return None, None
    text = win.read_text(errors="replace")
    m = re.search(r"^\s*exclude_bands\s*[=:]\s*([^!#\n]+)", text, re.I | re.M)
    if not m:
        return None, None
    excluded = set()
    for part in re.split(r"[,\s]+", m.group(1).strip()):
        a, _, b = part.partition("-")
        if a.isdigit():
            excluded.update(range(int(a), int(b or a) + 1))
    first = 1
    while first in excluded:
        first += 1
    return first, str(win)


# Choose the contiguous band mapping: from exclude_bands when it fits, else the offset with the smallest MAE.
def auto_band_start(reference, wannier, hr_path, emin, emax):
    first, win = excluded_offset(hr_path)
    fits = []
    for start in range(1, len(reference) - len(wannier) + 2):
        try:
            fits.append((error_metrics(reference, wannier, start, emin, emax)["mae_meV"], start))
        except ValueError:
            continue
    if not fits:
        return None, "no contiguous band mapping has reference energies in the window"
    best = min(fits)[1]
    if first is not None and any(start == first for _, start in fits):
        how = f"band offset {first} from exclude_bands in {win}, not verified"
        if first != best:
            how += f" (smallest-error offset would be {best})"
        return first, how
    return best, f"band offset {best} chosen by smallest error, not verified"


Prepared = namedtuple("Prepared", "hr_path run reference lineage chain structure checks failed data spin fermi fermi_source "
                                   "hr convention out emin emax notes")


# Spin channel "up" with a note when the reference INCAR is spin-polarised and collinear (ISPIN 2, no SOC, no LNONCOLLINEAR).
def collinear_reference_spin(reference):
    from pymatgen.io.vasp import Incar

    incar_path = Path(reference["path"]) / "INCAR"
    if incar_path.is_file() and collinear_polarised(Incar.from_file(incar_path)):
        return "up", "spin-up channel of the reference used for the k-path and Fermi level"
    return None, None


# Shared setup of a Wannier comparison or projection: Hamiltonian, output guard, reference bands, Fermi level and H(R).
def prepare(project, hr_file=None, run_name=None, reference_run=None, emin=-3.0, emax=3.0, out_dir=None, fermi_eV=None,
            latest_if_ambiguous=False, spin=None, structure_only=None, use_ws_distance=None, ws_fallback=False):
    project_info = db.get_project(project)
    if project_info is None:
        raise ValueError(f"Project {project} is not indexed; run altermag import first.")
    emin, emax = float(emin), float(emax)
    if not np.isfinite([emin, emax]).all() or emin >= emax:
        raise ValueError("The energy window must have finite emin < emax.")
    runs = db.get_runs(project)
    hr_path, run = locate_hamiltonian(project_info, runs, hr_file, run_name)
    out = Path(out_dir).expanduser() if out_dir else bands.default_out_dir(run["path"], project)
    bands.refuse_inside_root(out, [p["root"] for p in db.list_projects()],
                             "Refusing to write inside a scanned data root; choose --out-dir outside all indexed projects.")
    reference, lineage, chain, structure = select_reference(run, runs, reference_run, latest_if_ambiguous)
    notes = [lineage["selection_note"]] if lineage["selection_note"] else []
    checks, failed = None, []
    if structure_only is not None:
        structure, checks, failed = check_inputs(run, reference, structure_only)
    if spin in (None, "auto"):
        spin, note = collinear_reference_spin(reference)
        notes += [note] if note else []
    data = load_reference(Path(reference["path"]), structure, spin)
    fermi, fermi_source = common_fermi(chain, lineage, data, reference, fermi_eV)
    hr = read_hr(hr_path)
    fallback = ws_fallback and use_ws_distance is None
    if fallback and not companion_files(hr_path)[1].is_file():
        use_ws_distance = False
    try:
        terms, convention = interpolation_terms(hr_path, hr, use_ws_distance)
    except ValueError as exc:
        if not (fallback and use_ws_distance is None and "wsvec" in str(exc)):
            raise
        terms, convention = interpolation_terms(hr_path, hr, False)
        notes.append(f"seedname_wsvec.dat unusable ({exc}); retried without Wigner-Seitz distances")
        convention["convention_note"] = convention["convention_note"] or WS_FALLBACK_NOTE
    hr = terms
    if convention["convention_note"]:
        notes.append(convention["convention_note"])
    return Prepared(hr_path, run, reference, lineage, chain, structure, checks, failed, data, spin, fermi, fermi_source,
                    hr, convention, out, emin, emax, notes)


# Validate lineage and inputs, calculate the overlay, and write only outside every indexed data root.
def compare_wannier(project, hr_file=None, run_name=None, reference_run=None, emin=-3.0, emax=3.0,
                    out_dir=None, structure_only=False, fermi_eV=None, spin=None, band_start=None,
                    window=None, use_ws_distance=None, latest_if_ambiguous=False, ws_fallback=False):
    if spin not in (None, "up", "down"):
        raise ValueError("spin must be up or down.")
    if not isinstance(structure_only, bool):
        raise ValueError("structure_only must be a boolean.")
    prep = prepare(project, hr_file, run_name, reference_run, emin, emax, out_dir, fermi_eV, latest_if_ambiguous, spin,
                   structure_only, use_ws_distance, ws_fallback)
    run, reference, data, hr, convention = prep.run, prep.reference, prep.data, prep.hr, prep.convention
    hr_path, failed, fermi, emin, emax = prep.hr_path, prep.failed, prep.fermi, prep.emin, prep.emax
    wannier, residual = interpolate(hr, data["kfrac"])
    data["energies"] = data["energies"] - fermi
    wannier -= fermi
    metric_emin, metric_emax = (-float(window), float(window)) if window else (emin, emax)
    mapping_note = f"band offset {band_start} supplied with --band-start"
    if band_start == "auto":
        band_start, mapping_note = auto_band_start(data["energies"], wannier, hr_path, metric_emin, metric_emax)
    metrics = error_metrics(data["energies"], wannier, band_start, metric_emin, metric_emax)
    if metrics is not None:
        metrics["mapping_source"] = mapping_note
    originals = bands.existing_figure(reference["path"], project, reference["name"])
    source_figures = [str(originals)] if originals else []
    source_figures += [str(p) for p in sorted(Path(reference["path"]).iterdir())
                       if p.is_file() and p.suffix.lower() in (".png", ".jpg", ".jpeg", ".svg", ".pdf")]
    notes = ["Absolute VASP and Wannier energies share one energy zero; no fitted energy shift was applied.",
             "Parent links are indexed evidence and may be inferred; inspect lineage.parent_notes.",
             "Input equality and visual agreement do not establish physical convergence."]
    if metrics is None:
        notes.append("No bandwise errors: supply band_start only when the band correspondence is known.")
    else:
        notes.append(metrics["mapping_source"])
    if failed:
        notes.append("Preliminary overlay; input checks failed: " + ", ".join(failed))
    notes += prep.notes
    result = {"project": project, "wannier_run": run["name"], "reference_run": reference["name"],
              "hamiltonian": str(hr_path), "reference_source": data["source"], "original_figures": source_figures,
              "evidence_state": "preliminary: structure matched; required inputs differ or are missing" if failed
              else "input-matched comparison; interpolation quality requires researcher review",
              "lineage": prep.lineage, "input_checks": prep.checks, "interpolation": convention,
              "kpath_evidence": data["kpath_evidence"], "n_kpoints": len(data["x"]), "num_wann": hr["num_wann"],
              "spin": prep.spin, "fermi_eV": fermi, "fermi_source": prep.fermi_source, "emin": emin, "emax": emax,
              "metric_emin": metric_emin, "metric_emax": metric_emax,
              "max_hermiticity_residual_eV": residual, "metrics": metrics, "notes": notes}
    folder = prep.out / (bands.slug(run["name"]) + "_wannier_compare_" + uuid4().hex)
    folder.mkdir(parents=True, exist_ok=False)
    paths = {kind: str(folder / ("overlay." + kind)) for kind in ("png", "pdf", "svg", "json", "csv")}
    result["paths"] = paths
    draw_overlay(data, wannier, result, paths)
    with Path(paths["csv"]).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["k_index", "kx_frac", "ky_frac", "kz_frac", "distance_inv_angstrom"]
                        + [f"vasp_{i + 1}_eV" for i in range(len(data["energies"]))]
                        + [f"wannier_{i + 1}_eV" for i in range(len(wannier))])
        for ik, k in enumerate(data["kfrac"]):
            writer.writerow([ik + 1, *k, data["x"][ik], *data["energies"][:, ik], *wannier[:, ik]])
    Path(paths["json"]).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    return result


BAND_MATCH_WINDOW_EV = 1.0


# Raw or symmetrised: a file under wannier_sym or with a "Symmed by" header is the symmetrised Hamiltonian.
def hr_kind(path):
    with Path(path).open(errors="replace") as stream:
        header = stream.readline()
    return "symmetrised" if "wannier_sym" in Path(path).parts or "symm" in header.lower() else "raw"


# Run compare_wannier with the direct-HR fallback when wsvec is absent, and as a preliminary overlay only if allowed.
def run_comparison(project, allow_preliminary=False, **kw):
    notes = []
    for _ in range(2):
        try:
            res = compare_wannier(project, ws_fallback=True, **kw)
            break
        except ValueError as exc:
            text = str(exc)
            if "Required input equality failed" in text and allow_preliminary and not kw.get("structure_only"):
                kw["structure_only"] = True
                notes.append(text.split(" Use --structure-only")[0])
            else:
                raise
    res["fallbacks"] = notes
    res["hr_kind"] = hr_kind(res["hamiltonian"])
    return res


# Cache key of one Hamiltonian file: path, size and modification time.
def hr_key(path):
    st = Path(path).stat()
    return f"v2|{path}|{st.st_size}|{st.st_mtime_ns}"


# Band match of one Wannier run from the index cache; computes and caches it when compute is True, else returns None if not cached.
def band_match(project, run_name, compute=True, window=BAND_MATCH_WINDOW_EV):
    runs = db.get_runs(project)
    run = next((r for r in runs if r["name"] == run_name), None)
    hr = Path(run["path"]) / "wannier90_hr.dat" if run else None
    if hr is None or not hr.is_file():
        return {"status": "no_hr"}
    key = hr_key(hr)
    cached = db.get_band_match(project, run_name, key)
    if cached is not None or not compute:
        return cached
    try:
        res = run_comparison(project, allow_preliminary=True, run_name=run_name, latest_if_ambiguous=True,
                             band_start="auto", window=window)
    except NoReference as exc:
        out = {"status": "no_reference", "reason": str(exc)}
    except Exception as exc:
        out = {"status": "error", "reason": f"{type(exc).__name__}: {exc}"}
    else:
        m = res["metrics"]
        out = {"status": "ok" if m else "no_mapping", "reference_run": res["reference_run"], "metrics": m,
               "window": [res["metric_emin"], res["metric_emax"]], "png": res["paths"]["png"],
               "json": res["paths"]["json"], "evidence_state": res["evidence_state"], "hr_kind": res["hr_kind"],
               "fallbacks": res["fallbacks"] + [n for n in [res["lineage"]["selection_note"]] if n]}
    db.set_band_match(project, run_name, key, out)
    return out


# True when the band match read wannier90_hr.dat: status "ok" or "no_mapping".
def hr_parsed(res):
    return bool(res) and res.get("status") in ("ok", "no_mapping")


# One-line text of a band-match result for the suggestions step (e).
def band_match_text(res):
    if res is None:
        return "band match: not computed yet (computed on first view)"
    st = res.get("status")
    prefix = "preliminary: " if str(res.get("evidence_state", "")).startswith("preliminary") else ""
    if st == "no_reference":
        return "band match: no ancestor band run with eigenvalues"
    if st == "error":
        return f"band match: comparison failed ({res['reason']})"
    if st == "no_mapping":
        return f"{prefix}band match: no band mapping found against {res['reference_run']} (figure {display_path(res['png'])})"
    if st != "ok":
        return "band match: no wannier90_hr.dat"
    m, (lo, hi) = res["metrics"], res["window"]
    return (f"{prefix}band match: MAE {m['mae_meV']:.3f} meV, RMSE {m['rmse_meV']:.3f} meV, "
            f"95th percentile {m['p95_meV']:.3f} meV, max {m['max_meV']:.3f} meV, "
            f"{m['samples_over_20_meV']} of {m['samples']} samples above 20 meV "
            f"in the window {lo:g} to +{hi:g} eV against {res['reference_run']}, {res['hr_kind']} Hamiltonian; "
            f"{m.get('mapping_source', 'band offset source not recorded')} (figure {display_path(res['png'])})")
