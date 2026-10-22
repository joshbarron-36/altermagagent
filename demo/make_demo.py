"""Write a synthetic CrSe VASP project tree for AltermagAgent demos (no VASP needed)."""

import math
import os
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "CrSe_demo"
NEGATIVE_ROOT = HERE / "CrSe_demo_without_nosoc_bands"

A_POS, C_POS = 3.700, 6.050      # starting lattice (POSCAR of relax)
A_REL, C_REL = 3.684, 6.019      # relaxed lattice (CONTCAR of relax and later POSCARs)
EFERMI = 5.0                     # Fermi level used in the bands vasprun.xml

# name, start, end, status kind, energy (eV), INCAR text
RUNS = [
    ("relax", "2026-08-02 09:14", "2026-08-02 13:52", "ok", -36.8814,
     "SYSTEM = CrSe relax\nISPIN = 2\nMAGMOM = 3 -3 0 0\nENCUT = 500\nISIF = 3\nIBRION = 2\nNSW = 40\nEDIFFG = -0.01\nISMEAR = 0\nSIGMA = 0.05\n"),
    ("scf_nosoc", "2026-08-05 10:03", "2026-08-05 11:26", "ok", -36.9027,
     "SYSTEM = CrSe scf\nISPIN = 2\nMAGMOM = 3 -3 0 0\nENCUT = 500\nNSW = 0\nIBRION = -1\nLORBIT = 11\nLCHARG = .TRUE.\nLWAVE = .TRUE.\nISMEAR = 0\nSIGMA = 0.05\n"),
    ("bands_nosoc", "2026-08-06 08:40", "2026-08-06 09:18", "ok", -36.9025,
     "SYSTEM = CrSe bands\nISPIN = 2\nMAGMOM = 3 -3 0 0\nENCUT = 500\nICHARG = 11\nNSW = 0\nIBRION = -1\nNBANDS = 8\nLORBIT = 11\nISMEAR = 0\nSIGMA = 0.05\n"),
    ("scf_soc", "2026-09-01 09:30", "2026-09-01 15:07", "ok", -37.0412,
     "SYSTEM = CrSe scf SOC\nISPIN = 2\nLSORBIT = .TRUE.\nMAGMOM = 0 0 3 0 0 -3 0 0 0 0 0 0\nSAXIS = 0 0 1\nISYM = 0\nENCUT = 500\nNSW = 0\nIBRION = -1\nLCHARG = .TRUE.\nLWAVE = .TRUE.\n"),
    ("bands_soc", "2026-09-03 10:11", "2026-09-03 10:24", "err", -37.0408,
     "SYSTEM = CrSe bands SOC\nISPIN = 2\nLSORBIT = .TRUE.\nMAGMOM = 0 0 3 0 0 -3 0 0 0 0 0 0\nSAXIS = 0 0 1\nISYM = 0\nENCUT = 500\nICHARG = 11\nIBRION = 2\nNSW = 10\n"),
    ("wannier", "2026-09-14 09:05", "2026-09-14 12:48", "ok", -37.0410,
     "SYSTEM = CrSe wannier\nISPIN = 2\nLSORBIT = .TRUE.\nMAGMOM = 0 0 3 0 0 -3 0 0 0 0 0 0\nSAXIS = 0 0 1\nISYM = 0\nENCUT = 500\nISTART = 1\nNSW = 0\nIBRION = -1\nNBANDS = 48\nLWANNIER90 = .TRUE.\n"),
    ("bands_soc_retry", None, None, "run", -37.0405,
     "SYSTEM = CrSe bands SOC retry\nISPIN = 2\nLSORBIT = .TRUE.\nMAGMOM = 0 0 3 0 0 -3 0 0 0 0 0 0\nSAXIS = 0 0 1\nISYM = 0\nENCUT = 500\nICHARG = 11\nIBRION = -1\nNSW = 0\n"),
]


# Return POSCAR text for hexagonal NiAs-type CrSe with lattice a, c.
def poscar(a, c, title="CrSe"):
    h = a * math.sqrt(3) / 2
    v = lambda x: f"{x:17.10f}"
    return "\n".join([
        title, "   1.00000000000000",
        v(a / 2) + v(-h) + v(0), v(a / 2) + v(h) + v(0), v(0) + v(0) + v(c),
        "   Cr   Se", "     2     2", "Direct",
        "  0.0000000000  0.0000000000  0.0000000000",
        "  0.0000000000  0.0000000000  0.5000000000",
        "  0.3333333333  0.6666666667  0.2500000000",
        "  0.6666666667  0.3333333333  0.7500000000", ""])


# Return a magnetization block in OUTCAR format for one component.
def mag_block(comp, m):
    d = m - 0.035
    row = lambda i, s, p, dd, el: f"{i:6d}{s:12.3f}{p:8.3f}{dd:8.3f}{s + p + dd:8.3f}"
    return "\n".join([
        f" magnetization ({comp})", "", " # of ion       s       p       d       tot",
        " ------------------------------------------",
        row(1, 0.011 * (m != 0), 0.024 * (m != 0), d if m else 0, "Cr"),
        row(2, -0.011 * (m != 0), -0.024 * (m != 0), -d if m else 0, "Cr"),
        row(3, -0.001, 0.003, 0, "Se"), row(4, 0.001, -0.003, 0, "Se"),
        " --------------------------------------------------",
        f"tot{0:15.3f}{0:8.3f}{0:8.3f}{0:8.3f}", ""])


# Return OUTCAR text for a run of the given kind.
def outcar(name, kind, energy, soc):
    head = [" vasp.6.4.2 20Jul23 (build Aug 10 2026 14:02:11) complex", "",
            " POSCAR = CrSe", " running on   32 total cores", ""]
    if kind == "err":
        return "\n".join(head + [
            " ------------------------ aborting loop because EDIFF is reached ----------------------------------------",
            f"  free  energy   TOTEN  =  {energy:16.8f} eV",
            " BRION: g(F)=  0.218E-04 g(S)=  0.000E+00 retain N=  1 mean eig= 0.35",
            " trial: gam= 0.00000 g(F)=  0.221E-04 g(S)=  0.000E+00 ort =-0.211E-05 (trialstep = 0.100E+01)",
            " bond charge predicted", ""])
    if kind == "run":
        return "\n".join(head + [
            " charge density read from file: CHGCAR",
            " ----------------------------------------- Iteration    1(  31)  ---------------------------------------",
            "    POTLOK:  cpu time    0.6021: real time    0.6034",
            "    EDDAV:   cpu time   41.2210: real time   41.3018", ""])
    body = []
    steps = [energy + 0.0213, energy] if name == "relax" else [energy]
    for i, e in enumerate(steps):
        m = 3.08 if name == "relax" else 3.12 if not soc else 3.10
        if i < len(steps) - 1:
            m -= 0.05
        body += ["  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)", "  ---------------------------------------------------",
                 f"  free  energy   TOTEN  =  {e:16.8f} eV", f"  energy  without entropy=  {e + 0.0011:16.8f}  energy(sigma->0) =  {e + 0.0005:16.8f}", ""]
        if soc:
            body += [mag_block("x", 0), mag_block("y", 0), mag_block("z", m)]
        else:
            body += [mag_block("x", m)]
    tail = ["", " General timing and accounting informations for this job:",
            " ========================================================", "",
            "                  Total CPU time used (sec):     4870.112",
            "                            User time (sec):     4801.220",
            "                          Elapsed time (sec):     4902.551", ""]
    return "\n".join(head + body + tail)


# Format a number in the VASP OSZICAR F= style, for example -.36881400E+02.
def vf(x):
    ex = math.floor(math.log10(abs(x))) + 1
    s = f"{x / 10 ** ex:.8f}E{'-' if ex < 0 else '+'}{abs(ex):02d}"
    return s.replace("0.", ".", 1)


# Return OSZICAR text; running runs have no final F= line yet.
def oszicar(kind, energy, soc):
    lines = ["       N       E                     dE             d eps       ncg     rms          rms(c)"]
    for n in range(1, 13 if kind in ("ok", "err") else 7):
        lines.append(f"DAV:{n:4d}    {energy + 0.5 / n ** 3:+.6E}   {-0.5 / n ** 4:+.2E}   {-0.05 / n ** 4:+.2E}  2208   {0.2 / n ** 2:.2E}")
    mag = "mag=     0.0000     0.0000     0.0003" if soc else "mag=     0.0003"
    if kind in ("ok", "err"):
        lines.append(f"   1 F= {vf(energy)} E0= {vf(energy + 0.0005)}  d E =-.98000000E-05  {mag}")
    return "\n".join(lines) + "\n"


# Return spin up and down eigenvalues per k-point for the bands run.
def eigenvalues(nk=6, nb=8):
    centres = [-4.0, -2.5, -1.2, 4.3, 5.6, 6.6, 8.0, 9.5]
    splits = [0.05, 0.08, 0.10, 0.42, 0.18, 0.12, 0.09, 0.04]
    up, dn = [], []
    for k in range(nk):
        t = k / (nk - 1)
        shape = (0.0, 0.3, 0.7, 1.0, 0.6, 0.0)[k] if nk == 6 else math.sin(math.pi * t) ** 2  # zero at the path ends
        ku, kd = [], []
        for b in range(nb):
            e = EFERMI + centres[b] - 5.0 + 0.6 * math.cos(2 * math.pi * t + b)
            s = splits[b] * shape / 2
            ku.append(e + s)
            kd.append(e - s)
        up.append(ku)
        dn.append(kd)
    return up, dn


# Return a small but valid spin-polarised vasprun.xml for the bands run.
def vasprun(a, c):
    h = a * math.sqrt(3) / 2
    kpts = [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.2, 0.05, 0.0), (0.333, 0.1, 0.0), (0.4, 0.05, 0.1), (0.5, 0.0, 0.0)]
    up, dn = eigenvalues(len(kpts))
    frac = [(0, 0, 0), (0, 0, 0.5), (1 / 3, 2 / 3, 0.25), (2 / 3, 1 / 3, 0.75)]
    rec = 1 / a

    # Return a <structure> block, named or anonymous.
    def structure(nm=None):
        tag = f' <structure name="{nm}">' if nm else "  <structure>"
        return "\n".join([
            tag, "  <crystal>", '   <varray name="basis" >',
            f"    <v> {a / 2:16.8f} {-h:16.8f} {0:16.8f} </v>", f"    <v> {a / 2:16.8f} {h:16.8f} {0:16.8f} </v>",
            f"    <v> {0:16.8f} {0:16.8f} {c:16.8f} </v>", "   </varray>",
            f'   <i name="volume">{a * a * h * c:16.8f} </i>', '   <varray name="rec_basis" >',
            f"    <v> {rec:16.8f} {-rec / math.sqrt(3):16.8f} 0.0 </v>", f"    <v> {rec:16.8f} {rec / math.sqrt(3):16.8f} 0.0 </v>",
            f"    <v> 0.0 0.0 {1 / c:16.8f} </v>", "   </varray>", "  </crystal>", '  <varray name="positions" >']
            + [f"   <v> {x:16.8f} {y:16.8f} {z:16.8f} </v>" for x, y, z in frac] + ["  </varray>", " </structure>"])

    # Return the eigenvalue <set> lines of one spin channel.
    def eig_set(vals, spin):
        out = [f'      <set comment="spin {spin}">']
        for k, ev in enumerate(vals):
            out.append(f'       <set comment="kpoint {k + 1}">')
            out += [f"        <r> {e:10.4f} {1.0 if e < EFERMI else 0.0:8.4f} </r>" for e in ev]
            out.append("       </set>")
        out.append("      </set>")
        return out

    lines = ['<?xml version="1.0" encoding="ISO-8859-1"?>', "<modeling>", " <generator>",
             '  <i name="program" type="string">vasp </i>', '  <i name="version" type="string">6.4.2  </i>',
             '  <i name="subversion" type="string">20Jul23 complex  parallel </i>',
             '  <i name="platform" type="string">LinuxIFC </i>', '  <i name="date" type="string">2026 08 06 </i>',
             '  <i name="time" type="string">09:18:00 </i>', " </generator>", " <incar>",
             '  <i type="string" name="SYSTEM">CrSe bands</i>', '  <i type="int" name="ISPIN">     2</i>',
             '  <i type="int" name="ICHARG">    11</i>', '  <i type="int" name="NBANDS">     8</i>',
             '  <i name="ENCUT">    500.00000000</i>', " </incar>", " <kpoints>",
             '  <generation param="listgenerated">', '   <i name="divisions" type="int">       5 </i>',
             '   <v>       0.00000000       0.00000000       0.00000000 </v>',
             '   <v>       0.50000000       0.00000000       0.00000000 </v>', '  </generation>',
             '  <varray name="kpointlist" >'] + [f"   <v> {x:16.8f} {y:16.8f} {z:16.8f} </v>" for x, y, z in kpts] + [
             "  </varray>", '  <varray name="weights" >'] + ["   <v>       0.16666667 </v>" for _ in kpts] + [
             "  </varray>", " </kpoints>", " <parameters>", '  <separator name="electronic" >',
             '   <i type="string" name="PREC">accurate</i>', '   <i name="ENMAX">    500.00000000</i>',
             '   <i name="NELECT">     24.00000000</i>', '   <i type="int" name="NELM">    60</i>', '   <i name="EDIFF">      0.00000100</i>', '   <i type="int" name="NBANDS">     8</i>',
             '   <separator name="electronic spin" >', '    <i type="int" name="ISPIN">     2</i>',
             '    <i type="logical" name="LNONCOLLINEAR"> F  </i>',
             '    <v name="MAGMOM">      3.00000000     -3.00000000      0.00000000      0.00000000</v>',
             '    <i type="logical" name="LSORBIT"> F  </i>', "   </separator>", "  </separator>",
             '  <separator name="ionic" >', '   <i type="int" name="NSW">     0</i>', '   <i type="int" name="IBRION">    -1</i>',
             "  </separator>", '  <i type="int" name="ICHARG">    11</i>', '  <i type="logical" name="LEPSILON"> F  </i>',
             " </parameters>", " <atominfo>", "  <atoms>       4 </atoms>", "  <types>       2 </types>",
             '  <array name="atoms" >', '   <dimension dim="1">ion</dimension>',
             '   <field type="string">element</field>', '   <field type="int">atomtype</field>', "   <set>",
             "    <rc><c>Cr</c><c>   1</c></rc>", "    <rc><c>Cr</c><c>   1</c></rc>",
             "    <rc><c>Se</c><c>   2</c></rc>", "    <rc><c>Se</c><c>   2</c></rc>", "   </set>", "  </array>",
             '  <array name="atomtypes" >', '   <dimension dim="1">type</dimension>',
             '   <field type="int">atomspertype</field>', '   <field type="string">element</field>',
             '   <field>mass</field>', '   <field>valence</field>', '   <field type="string">pseudopotential</field>',
             "   <set>", "    <rc><c>   2</c><c>Cr</c><c>     52.00000000</c><c>      6.00000000</c><c>  PAW_PBE Cr 06Sep2000                  </c></rc>",
             "    <rc><c>   2</c><c>Se</c><c>     78.96000000</c><c>      6.00000000</c><c>  PAW_PBE Se 06Sep2000                  </c></rc>",
             "   </set>", "  </array>", " </atominfo>", structure("initialpos"), " <calculation>",
             "  <scstep>", "   <energy>", '    <i name="e_fr_energy">    -36.90250000 </i>',
             '    <i name="e_wo_entrp">    -36.90160000 </i>', '    <i name="e_0_energy">    -36.90205000 </i>', "   </energy>",
             "  </scstep>", structure(), "  <energy>",
             '   <i name="e_fr_energy">    -36.90250000 </i>', '   <i name="e_wo_entrp">    -36.90160000 </i>',
             '   <i name="e_0_energy">    -36.90205000 </i>', "  </energy>", "  <eigenvalues>", "   <array>",
             '    <dimension dim="1">band</dimension>', '    <dimension dim="2">kpoint</dimension>',
             '    <dimension dim="3">spin</dimension>', '    <field>eigene</field>', '    <field>occ</field>', "    <set>"]
    lines += eig_set(up, 1) + eig_set(dn, 2)
    lines += ["    </set>", "   </array>", "  </eigenvalues>", "  <dos>", f'   <i name="efermi">     {EFERMI:.8f} </i>'] + total_dos() + [
              "  </dos>", " </calculation>", structure("finalpos"), "</modeling>", ""]
    return "\n".join(l for l in lines if l != "")


# Return a coarse total DOS block (two spins, five energies) so parsers find E_F.
def total_dos():
    out = ["   <total>", "    <array>", '     <dimension dim="1">gridpoints</dimension>', '     <dimension dim="2">spin</dimension>',
           "     <field>energy</field>", "     <field>total</field>", "     <field>integrated</field>", "     <set>"]
    for spin in (1, 2):
        out.append(f'      <set comment="spin {spin}">')
        out += [f"       <r> {EFERMI + e:10.4f} {1.2:10.4f} {6.0 * (i + 1):10.4f} </r>" for i, e in enumerate((-2, -1, 0, 1, 2))]
        out.append("      </set>")
    return out + ["     </set>", "    </array>", "   </total>"]


# Draw a bands figure for bands_nosoc with matplotlib (shipped with pymatgen).
def bands_png(path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    up, dn = eigenvalues(6)
    fig, ax = plt.subplots(figsize=(6, 3.2), dpi=100)
    for b in range(len(up[0])):
        ax.plot(range(6), [k[b] - EFERMI for k in up], color="#dc2626", lw=1.5)
        ax.plot(range(6), [k[b] - EFERMI for k in dn], color="#2563eb", lw=1.5)
    ax.axhline(0, color="#5b626b", ls="--", lw=1)
    ax.set_ylim(-3, 3)
    ax.set_ylabel("E - E_F (eV)")
    ax.set_xticks([0, 5], ["Γ", "M"])
    ax.set_title("CrSe bands without SOC, spin up red, spin down blue", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return True


# Write one text file.
def put(folder, name, text):
    (folder / name).write_text(text)


# Write every file of one run and return the list of (input?, path) pairs for time stamping.
def write_run(name, kind, energy, incar, root=ROOT):
    d = root / name
    d.mkdir(parents=True)
    soc = "LSORBIT" in incar
    inputs, outputs = ["INCAR", "KPOINTS", "POSCAR"], []
    put(d, "INCAR", incar)
    if name.startswith("bands"):
        put(d, "KPOINTS", "Line mode G-M-K-G-A\n60\nLine-mode\nReciprocal\n 0.0 0.0 0.0 ! G\n 0.5 0.0 0.0 ! M\n")
    else:
        put(d, "KPOINTS", "Automatic mesh\n0\nGamma\n 16 16 10\n 0 0 0\n")
    put(d, "POSCAR", poscar(A_POS, C_POS) if name == "relax" else poscar(A_REL, C_REL))
    if kind != "run":
        put(d, "CONTCAR", poscar(A_REL, C_REL))
        outputs.append("CONTCAR")
    put(d, "OUTCAR", outcar(name, kind, energy, soc))
    put(d, "OSZICAR", oszicar(kind, energy, soc))
    outputs += ["OUTCAR", "OSZICAR"]
    log = [" running on   32 total cores", " vasp.6.4.2 20Jul23 (build Aug 10 2026 14:02:11) complex",
           " POSCAR found :  2 types and       4 ions", " LDA part: xc-table for Pade appr. of Perdew"]
    if kind == "err":
        log += ["   1 F= " + vf(energy), " ZBRENT: fatal error in bracketing",
                "     please rerun with smaller EDIFF, or copy CONTCAR", "     to POSCAR and continue"]
    elif kind == "run":
        log += [" charge density read from file: CHGCAR", " entering main loop",
                "DAV:  31    -0.370405E+02   -0.17E-03   -0.31E-02  8064   0.19E-01"]
    else:
        log += ["   1 F= " + vf(energy), " writing wavefunctions"]
    put(d, "vasp.log", "\n".join(log) + "\n")
    outputs.append("vasp.log")
    if name.startswith("scf"):
        put(d, "CHGCAR", "placeholder charge density, reference only\n")
        put(d, "WAVECAR", "placeholder wavefunctions\n")
        outputs += ["CHGCAR", "WAVECAR"]
    if name.startswith("bands"):
        put(d, "CHGCAR", "placeholder charge density copied from the parent run\n")
        inputs.append("CHGCAR")
    if name == "bands_nosoc":
        put(d, "vasprun.xml", vasprun(A_REL, C_REL))
        outputs.append("vasprun.xml")
        if bands_png(d / "bands.png"):
            outputs.append("bands.png")
    if name == "wannier":
        put(d, "WAVECAR", "placeholder wavefunctions copied from the parent run\n")
        inputs += ["WAVECAR", "wannier90.win"]
        put(d, "wannier90.win", "num_wann        = 22\nnum_bands       = 48\nnum_iter        = 200\ndis_win_max     = 8.0\n"
            "dis_froz_max    = 2.0\nwrite_hr        = .true.\n\nbegin projections\nCr: d\nSe: p\nend projections\nmp_grid         = 8 8 6\n")
        put(d, "wannier90.wout", "            <<< Wannierisation convergence criteria satisfied >>>\n\n Final State\n"
            "    Final Spread (Ang^2)       Omega Total  =    38.451236700\n All done: wannier90 exiting\n")
        hr = [" written on 14Sep2026 at 12:48:03", "          22", "         397",
              "    4    4    4    2    2    2    1    1    1    1    1    1    1    1    1",
              "   -3   -3   -2    1    1   -0.000412    0.000000", "   -3   -3   -2    2    1    0.000037    0.000000",
              "    0    0    0    1    1    3.214765    0.000000", "    0    0    0    2    1   -0.018342    0.000000"]
        put(d, "wannier90_hr.dat", "\n".join(hr) + "\n")
        outputs += ["wannier90.wout", "wannier90_hr.dat"]
    return d, inputs, outputs


# Turn a "YYYY-MM-DD HH:MM" string into a POSIX time stamp.
def ts(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M").timestamp()


# Build the whole demo tree, replacing an earlier copy; --without-nosoc-bands writes the negative fixture without bands_nosoc.
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Write the synthetic CrSe demo tree.")
    ap.add_argument("--without-nosoc-bands", action="store_true",
                    help=f"leave out the bands_nosoc run (negative fixture), default folder {NEGATIVE_ROOT.name}")
    ap.add_argument("--out", help="output folder, default demo/CrSe_demo or demo/CrSe_demo_without_nosoc_bands")
    args = ap.parse_args(argv)
    root = Path(args.out).resolve() if args.out else NEGATIVE_ROOT if args.without_nosoc_bands else ROOT
    runs = [r for r in RUNS if not (args.without_nosoc_bands and r[0] == "bands_nosoc")]
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    now = time.time()
    for name, start, end, kind, energy, incar in runs:
        d, inputs, outputs = write_run(name, kind, energy, incar, root)
        t0 = ts(start) if start else now - 2 * 3600
        t1 = ts(end) if end else now - 10 * 60
        for f in inputs:
            os.utime(d / f, (t0, t0))
        for f in outputs:
            os.utime(d / f, (t1, t1))
        os.utime(d, (t1, t1))
    print(f"Demo tree written to {root} with {len(runs)} runs")


if __name__ == "__main__":
    sys.exit(main())
