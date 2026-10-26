"""Spin projection of Wannier bands: expectation values of sigma_x, sigma_y and sigma_z along the band k-path.

Based on Diana's reconstruct_spin_sz_bands.py. The colours assume a diagonal spin operator in the Wannier basis;
they are not interpolated ab initio spin matrices (.spn).
"""

import json
import re
from pathlib import Path
from uuid import uuid4

import numpy as np

from . import bands
from . import wannier_compare as wc
from .branches import noncollinear
from .files import display_path
from .parsers_extra import parse_wannier_win

DEGENERACY_TOL_EV = 1e-5
ORDERS = ("udud", "uudd")
COMPONENTS = ("x", "y", "z")
PNG_DPI = 120
NOT_CHECKED = "spin projections: not checked"
ROTATION = "each component diagonalised inside degenerate groups"
NO_WIN_NOTE = "wannier90.win not found; spinor basis assumed"


# Spin order stated in the first header line ("permuted uudd->udud") or in the file name (_udud, _uudd), else None.
def spin_order_from_hamiltonian(path):
    path = Path(path)
    with path.open(errors="replace") as stream:
        header = stream.readline().lower()
    stated = re.findall(r"->\s*(udud|uudd)\b", header) or re.findall(r"\b(udud|uudd)\b", header)
    if stated:
        return stated[-1]
    named = re.findall(r"(?:^|[_.\-])(udud|uudd)(?=[_.\-]|$)", path.name.lower())
    return named[-1] if named else None


# Spin order and its source: a given order beats the file, a stated order beats the udud assumption.
def resolve_order(path, given=None):
    detected = spin_order_from_hamiltonian(path)
    if given is not None:
        return given, "given" if detected in (None, given) else f"given; the file states {detected}"
    if detected is not None:
        return detected, "detected"
    return "udud", "assumed"


# Row indices of spin up and spin down for each orbital: udud is up, down per orbital; uudd is all up, then all down.
def spin_rows(nw, order):
    if order not in ORDERS:
        raise ValueError(f"Spin order must be udud or uudd, not {order!r}.")
    if nw < 2 or nw % 2:
        raise ValueError(f"A spinor Hamiltonian needs an even number of Wannier functions, not {nw}.")
    if order == "udud":
        return np.arange(0, nw, 2), np.arange(1, nw, 2)
    return np.arange(nw // 2), np.arange(nw // 2, nw)


# Reference Pauli matrices in the Wannier basis (nw x nw), placed on the rows given by spin_rows; used by tests only.
def spin_operators(nw, order):
    up, down = spin_rows(nw, order)
    ops = {c: np.zeros((nw, nw), complex) for c in COMPONENTS}
    ops["z"][up, up], ops["z"][down, down] = 1, -1
    ops["x"][up, down], ops["x"][down, up] = 1, 1
    ops["y"][up, down], ops["y"][down, up] = -1j, 1j
    return ops


# Pauli operator of one component times vectors (nw x g), as a signed row permutation of the up and down rows.
def apply_operator(vectors, rows, component):
    up, down = rows
    out = np.empty_like(vectors, dtype=complex)
    if component == "z":
        out[up], out[down] = vectors[up], -vectors[down]
    elif component == "x":
        out[up], out[down] = vectors[down], vectors[up]
    elif component == "y":
        out[up], out[down] = -1j * vectors[down], 1j * vectors[up]
    else:
        raise ValueError(f"Spin component must be x, y or z, not {component!r}.")
    return out


# Real expectation value of one Pauli component for every column of vectors.
def expectation(vectors, rows, component):
    return np.real(np.sum(vectors.conj() * apply_operator(vectors, rows, component), axis=0))


# Degenerate groups as (first, last) index ranges: each group spans at most tol above its first state.
def degenerate_groups(values, tol=DEGENERACY_TOL_EV):
    groups, first = [], 0
    while first < len(values):
        last = first + 1
        while last < len(values) and values[last] - values[first] <= tol:
            last += 1
        groups.append((first, last))
        first = last
    return groups


# Diagonalise one component inside a degenerate group (columns sub): its sorted expectations and the rotation.
def group_rotation(sub, rows, component):
    return np.linalg.eigh(sub.conj().T @ apply_operator(sub, rows, component))


# Eigenvalues and <sigma_x>, <sigma_y>, <sigma_z> of every band of one Hermitian H(k); each component is
# diagonalised on its own inside degenerate groups (spread at most tol), so every component is basis-independent.
def expectations(hk, rows, tol=DEGENERACY_TOL_EV):
    values, vectors = np.linalg.eigh(hk)
    spins = {c: expectation(vectors, rows, c) for c in COMPONENTS}
    if np.all(np.diff(values) > tol):
        return values, spins
    for first, last in degenerate_groups(values, tol):
        if last - first > 1:
            for c in COMPONENTS:
                spins[c][first:last] = group_rotation(vectors[:, first:last], rows, c)[0]
    return values, spins


# Energies (bands x k) and the three expectation arrays along the k-points, with the largest Hermiticity residual.
def project_path(hr, kfrac, order, tol=DEGENERACY_TOL_EV):
    nw, nk = hr["num_wann"], len(kfrac)
    rows = spin_rows(nw, order)
    energies = np.empty((nw, nk))
    spins = {c: np.empty((nw, nk)) for c in COMPONENTS}

    # Store the energies and the three expectations of one k-point.
    def store(ik, hk):
        energies[:, ik], values = expectations(hk, rows, tol)
        for c in COMPONENTS:
            spins[c][:, ik] = values[c]

    residual = wc.hamiltonians(hr, kfrac, store)
    return energies, spins, residual


# Caption under the figure: the basis assumption, the order and its source, the rotation, the SOC warning and basis notes.
def caption_text(order, order_source, soc, basis_note=None):
    text = f"Colours assume a diagonal spin operator in the Wannier basis, order {order} ({order_source}); {ROTATION}"
    return text + ("; with SOC these are not the physical spin" if soc else "") + (f"; {basis_note}" if basis_note else "")


# Refuse a non-spinor Hamiltonian: spinors false in wannier90.win or an odd num_wann; returns a note when no .win exists.
def spinor_check(hr_path, run, num_wann):
    if num_wann % 2:
        raise ValueError(f"num_wann {num_wann} is odd, so {hr_path.name} is not a spinor Hamiltonian; "
                         "spin projection needs spinor Wannier functions.")
    candidates = [wc.companion_files(hr_path)[0], Path(run["path"]) / "wannier90.win"]
    win = next((path for path in candidates if path.is_file()), None)
    if win is None:
        return NO_WIN_NOTE
    settings = parse_wannier_win(win)
    if "error" in settings:
        raise ValueError(f"wannier90.win could not be read ({settings['error']}); spinor basis unknown.")
    if settings["spinors"] is False:
        raise ValueError(f"{win.name} has spinors = false (the Wannier90 default when absent): a collinear basis, "
                         "so the spin projection would be meaningless; use a spinor Wannier run.")
    return None


# Mean |sigma|, min and max per component over the states inside the energy window, or None when the window is empty.
def window_stats(energies, spins, emin, emax):
    mask = (energies >= emin) & (energies <= emax)
    if not mask.any():
        return {c: None for c in COMPONENTS}, 0
    return {c: {"mean_abs": float(np.abs(spins[c][mask]).mean()), "min": float(spins[c][mask].min()),
                "max": float(spins[c][mask].max())} for c in COMPONENTS}, int(mask.sum())


# Draw the bands inside the window as grey lines with coloured segments, one panel per component, and save the PNG.
def draw(data, energies, spins, components, result, png):
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.figure import Figure
    import matplotlib.cm as cm

    figure = Figure(figsize=(7.5, 3.6 * len(components) + 0.6))
    axes = figure.subplots(len(components), 1, sharex=True, squeeze=False)[:, 0]
    norm = Normalize(-1, 1)
    keep = [i for i, band in enumerate(energies) if band.max() >= result["emin"] and band.min() <= result["emax"]]
    for axis, c in zip(axes, components):
        grey, pieces, weights = [], [], []
        for i in keep:
            for (xs, ys), (_, ws) in zip(bands.segments(data["x"], energies[i], data["branches"]),
                                         bands.segments(data["x"], spins[c][i], data["branches"])):
                points = np.column_stack((xs, ys))
                grey.append(points)
                pieces.append(np.stack((points[:-1], points[1:]), axis=1))
                weights.append(0.5 * (ws[:-1] + ws[1:]))
        if keep:
            axis.add_collection(LineCollection(grey, colors="0.75", linewidths=0.4, zorder=1))
            axis.add_collection(LineCollection(np.concatenate(pieces), cmap="coolwarm_r", norm=norm,
                                               array=np.concatenate(weights), linewidth=1.0, zorder=2))
        bands.decorate_band_axes(axis, data, (result["emin"], result["emax"]), zero_colour="0.45")
        bar = figure.colorbar(cm.ScalarMappable(norm=norm, cmap="coolwarm_r"), ax=axis, pad=0.02)
        bar.set_label(rf"$\langle\sigma_{c}\rangle$ (assumed Wannier basis)")
        bar.set_ticks([-1, 0, 1])
    axes[-1].set_xlabel("Wave vector")
    axes[0].set_title(f"{result['wannier_run']}: spin projection along the k-path of {result['reference_run']}", fontsize=9)
    figure.text(0.5, 0.01, result["caption"], ha="center", va="bottom", fontsize=8, wrap=True)
    figure.tight_layout(rect=(0, 0.05 / len(components), 1, 1))
    figure.savefig(png, dpi=PNG_DPI)


# Spin projection of an indexed Wannier Hamiltonian along its ancestor band run's k-path: PNG, JSON and a result dict.
def spin_projection(project, hr_file=None, run_name=None, reference_run=None, emin=-1.0, emax=1.0, component="z",
                    order=None, out_dir=None, fermi_eV=None):
    prep = wc.prepare(project, hr_file, run_name, reference_run, emin, emax, out_dir, fermi_eV,
                      latest_if_ambiguous=True, spin="auto", ws_fallback=True)
    run, reference, data, hr = prep.run, prep.reference, prep.data, prep.hr
    if component not in (*COMPONENTS, "all"):
        raise ValueError(f"Spin component must be x, y, z or all, not {component!r}.")
    basis_note = spinor_check(prep.hr_path, run, hr["num_wann"])
    notes = prep.notes + ([basis_note] if basis_note else [])
    order, order_source = resolve_order(prep.hr_path, order)
    energies, spins, residual = project_path(hr, data["kfrac"], order)
    energies -= prep.fermi
    stats, states = window_stats(energies, spins, prep.emin, prep.emax)
    details = run.get("details") or {}
    soc = bool(details.get("soc"))
    result = {"wannier_run": run["name"], "reference_run": reference["name"], "hamiltonian": str(prep.hr_path),
              "spin_order": order, "spin_order_source": order_source, "soc": soc,
              "noncollinear": bool(noncollinear(details)), "component": component,
              "degeneracy_tolerance_eV": DEGENERACY_TOL_EV, "degenerate_rotation": ROTATION, "fermi_eV": prep.fermi, "fermi_source": prep.fermi_source,
              "emin": prep.emin, "emax": prep.emax, "kpath_source": str(Path(reference["path"]) / "KPOINTS"),
              "n_kpoints": len(data["x"]), "num_wann": hr["num_wann"], "max_hermiticity_residual_eV": residual,
              "interpolation": prep.convention, "states_in_window": states, "window_stats": stats,
              "caption": caption_text(order, order_source, soc, basis_note), "notes": notes}
    folder = prep.out / f"{bands.slug(run['name'])}_spin_projection_{uuid4().hex}"
    folder.mkdir(parents=True, exist_ok=False)
    stem = f"{bands.slug(prep.hr_path.name)}_spin_{component}_{prep.emin:g}_{prep.emax:g}"
    result["paths"] = {"png": str(folder / (stem + ".png")), "json": str(folder / (stem + ".json"))}
    draw(data, energies, spins, COMPONENTS if component == "all" else (component,), result, result["paths"]["png"])
    Path(result["paths"]["json"]).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    return result


# A Fermi source for printing: its leading path shortened by display_path when it exists, other text unchanged.
def source_text(source):
    path, sep, rest = str(source).partition(" (")
    return display_path(path) + sep + rest if path and Path(path).exists() else str(source)


# Printed lines of a result: caption, order, window means of |sigma|, Fermi level, k-points, residual and notes.
def summary_lines(result):
    means = ", ".join(f"{c} {s['mean_abs']:.3f}" if s else f"{c} none" for c, s in result["window_stats"].items())
    return [result["caption"],
            f"Spin order {result['spin_order']} ({result['spin_order_source']}); window {result['emin']:g} to "
            f"{result['emax']:g} eV, {result['states_in_window']} states; mean |sigma|: {means}",
            f"Fermi level {result['fermi_eV']:.4f} eV from {source_text(result['fermi_source'])}; "
            f"{result['n_kpoints']} k-points of {result['reference_run']}, {result['num_wann']} Wannier functions",
            f"Largest Hermiticity residual {result['max_hermiticity_residual_eV']:.2e} eV", *result["notes"]]


# Step (e) suggestion after the band comparison when it parsed wannier90_hr.dat, else None.
def spin_line(run_name, band_match_result):
    if not wc.hr_parsed(band_match_result):
        return None
    return (f"Next: spin projection of {run_name}/wannier90_hr.dat (sigma_x, sigma_y, sigma_z along the band path); "
            "confirm the Wannier spin order first")


# Step (e) spin entry: the run, whether the projection is available, and its evidence line.
def spin_entry(run_name, band_match_result):
    line = spin_line(run_name, band_match_result)
    return {"run": run_name, "available": line is not None, "line": line or NOT_CHECKED}


# Refresh the spin entry of a stored step (e) from a new band match; steps stored without "spin" use the literal line.
def refresh_spin_step(step, run_name, band_match_result):
    old = (step.get("spin") or {}).get("line", NOT_CHECKED)
    new = spin_entry(run_name, band_match_result)
    step["evidence"] = [new["line"] if e == old else e for e in step["evidence"]]
    step["spin"] = new
    return step
