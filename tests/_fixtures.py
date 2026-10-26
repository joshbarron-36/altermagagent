"""Shared test setup: the synthetic demo tree in a temporary folder and a temporary index home per test."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

from altermagagent import bands, db

REPO = Path(__file__).resolve().parents[1]


# Write the demo tree into a new temporary folder; returns the TemporaryDirectory and the tree path.
def demo_tree():
    tmp = tempfile.TemporaryDirectory()
    demo = Path(tmp.name) / "CrSe_demo"
    subprocess.run([sys.executable, str(REPO / "demo" / "make_demo.py"), "--out", str(demo)],
                   check=True, capture_output=True)
    return tmp, demo


# Point the index and default figure folder at a fresh home under parent for the running test, undone at cleanup; returns the home.
def temp_index(test, parent):
    home = Path(tempfile.mkdtemp(dir=parent)) / "home"
    for p in (patch.object(db, "HOME_DIR", home), patch.object(db, "DB_PATH", home / "index.sqlite"),
              patch.object(bands, "FIG_DIR", home / "figures"), patch.dict(os.environ, {"ALTERMAG_HOME": str(home)})):
        p.start()
        test.addCleanup(p.stop)
    return home


# Copy the demo tree to parent/name with a collinear (no SOC) Wannier run, its own Hamiltonian and .win, then import it.
# The Wannier run gets bands_nosoc (ISPIN 2, no SOC) as its parent; wsvec text, when given, is written next to the Hamiltonian.
def collinear_wannier_project(demo, parent, name, hr_text, win_text=None, wsvec_text=None):
    import shutil
    from altermagagent import importer

    root = Path(parent) / name
    shutil.copytree(demo, root)
    wannier = root / "wannier"
    shutil.copy(root / "bands_nosoc" / "POSCAR", wannier / "POSCAR")
    (wannier / "INCAR").write_text("SYSTEM = CrSe wannier no SOC\nISPIN = 2\nMAGMOM = 3 -3 0 0\nENCUT = 500\n"
                                   "NSW = 0\nIBRION = -1\nLWANNIER90 = .TRUE.\n")
    (wannier / "wannier90_hr.dat").write_text(hr_text)
    win = wannier / "wannier90.win"
    if win_text is None:
        win.unlink()
    else:
        win.write_text(win_text)
    if wsvec_text is not None:
        (wannier / "wannier90_wsvec.dat").write_text(wsvec_text)
    write_spin_polarised_bands(root / "bands_nosoc")
    importer.import_project(root, name)
    db.set_parent_choice(name, "wannier", "bands_nosoc")
    return root


# Replace a band run's path by a four-point G-M line with a two-channel EIGENVAL and FERMI_ENERGY 0, as ISPIN 2 writes it.
def write_spin_polarised_bands(folder):
    from pymatgen.io.vasp import Poscar
    from altermagagent.wannier_compare import reference_path

    (folder / "KPOINTS").write_text("G-M\n4\nLine-mode\nReciprocal\n 0.0 0.0 0.0 ! G\n 0.5 0.0 0.0 ! M\n")
    kfrac = reference_path(folder, Poscar.from_file(folder / "POSCAR").structure)["kfrac"]
    levels = [-3.0, -1.2, -0.8, 0.8, 1.2, 3.0, 4.0, 5.0]
    lines = ["    4    4    1    2", "  0.1E+02", "  0.0", "  CAR", " CrSe", f"  8  {len(kfrac)}  {len(levels)}"]
    for k in kfrac:
        lines += ["", "  " + "  ".join(f"{x:.7E}" for x in k) + "  0.2500000E+00"]
        lines += [f"{i + 1:5d}  {e:12.6f}  {e + 0.3:12.6f}  1.0  1.0" for i, e in enumerate(levels)]
    (folder / "EIGENVAL").write_text("\n".join(lines) + "\n")
    (folder / "FERMI_ENERGY").write_text("0.0\n")
