"""Builder of small synthetic vasprun.xml files with known eigenvalues and projection weights.

Run `python tests/fixtures/orbital_vasprun_builder.py` to rewrite the files in tests/fixtures/orbital/.
Every value is given by a formula below, so a test can decode which spin, k-point, band, ion and orbital a number came from.
Run `python tests/fixtures/orbital_vasprun_builder.py --trim <vasprun.xml>` to rewrite the fixture cut from a real band run.
The builder imports the parser's component label table, so the synthetic fixtures are not an independent check of the layout.
The fixture trimmed from a real VASP band run is that independent check.
"""

from __future__ import annotations

import sys
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from vasp_orbital_plotter.parser import (COMPONENT_LABELS, LAYOUT_COLLINEAR, LAYOUT_NONCOLLINEAR,  # noqa: E402
                                         eigen_spin_blocks)

ORBITALS_9 = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]
ORBITALS_16 = ORBITALS_9 + ["fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
SPECIES = ["Cr", "Se", "Cr"]
KPOINTS = [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.5, 0.25, 0.0)]
NBANDS = 3
EFERMI = 1.23456789
LATTICE = [(3.6, 0.0, 0.0), (-1.8, 3.11769145, 0.0), (0.0, 0.0, 5.9)]
RECIPROCAL = [(0.27777778, 0.16037507, 0.0), (0.0, 0.32075015, 0.0), (0.0, 0.0, 0.16949153)]
POSITIONS = [(0.0, 0.0, 0.0), (0.33333333, 0.66666667, 0.25), (0.0, 0.0, 0.5)]
OUT_DIR = Path(__file__).resolve().parent / "orbital"
TRIMMED_NAME = "mnte_9.2-band_trimmed_vasprun.xml"


# Known eigenvalue in eV of spin block s, k-point k and band b (all 0-based).
def eigenvalue(s: int, k: int, b: int) -> float:
    return round(-2.0 + 1.5 * b + 0.125 * k + 0.5 * s, 4)


# Known occupation of band b: the lowest band is full, the others empty.
def occupation(b: int) -> float:
    return 1.0 if b == 0 else 0.0


# Known weight of component c, k-point k, band b, ion i, orbital o (0-based); components 1..3 of noncollinear runs are signed.
def weight(c: int, k: int, b: int, i: int, o: int, signed: bool = False) -> float:
    value = 0.1 * b + 0.01 * (i + 1) + 0.001 * (o + 1) + 0.0001 * k + (0.3 if c == 1 and not signed else 0.0)
    if signed and c > 0:
        value = (-1) ** c * (value + 0.0002 * c)
    return round(value, 4)


# Text of one <i> or <v> parameter line in VASP's style.
def _param(name: str, value: object, kind: str | None = None) -> str:
    kind_attr = f' type="{kind}"' if kind else ""
    return f'   <i{kind_attr} name="{name}">{value}</i>'


# Text of a complete vasprun.xml for the given layout; options break it on purpose for the error tests.
def build(layout: str, orbitals: list[str] = ORBITALS_9, projected: bool = True, drop_band_at: tuple[int, int] | None = None,
          extra_projected_block: bool = False, labels: bool = False, band_outside_kpoint: bool = False) -> str:
    nonc = layout == LAYOUT_NONCOLLINEAR
    ispin = 2 if layout == LAYOUT_COLLINEAR else 1
    eig_blocks = eigen_spin_blocks(layout)
    proj_blocks = len(COMPONENT_LABELS[layout]) + (1 if extra_projected_block else 0)
    logical = lambda flag: " T  " if flag else " F  "
    out = ['<?xml version="1.0" encoding="ISO-8859-1"?>', "<modeling>", " <generator>",
           '  <i name="program" type="string">vasp </i>', '  <i name="version" type="string">6.3.2  </i>', " </generator>",
           " <incar>", '  <i type="int" name="ISPIN">     2</i>' if nonc else f'  <i type="int" name="ISPIN">     {ispin}</i>',
           f'  <i type="logical" name="LSORBIT">{logical(nonc)}</i>', f'  <i type="int" name="NBANDS">    {NBANDS}</i>',
           f'  <i type="int" name="LORBIT">    {11 if projected else 10}</i>', " </incar>", " <kpoints>"]
    out.append('  <generation param="listgenerated">')
    out.append('   <i name="divisions" type="int">       2 </i>')
    out.extend(f"   <v> {k[0]:16.8f} {k[1]:16.8f} {k[2]:16.8f} </v>" for k in (KPOINTS[0], KPOINTS[2], KPOINTS[3]))
    out.append("  </generation>")
    out.append('  <varray name="kpointlist" >')
    out.extend(f"   <v> {k[0]:16.8f} {k[1]:16.8f} {k[2]:16.8f} </v>" for k in KPOINTS)
    out.append("  </varray>")
    out.append('  <varray name="weights" >')
    out.extend(f"   <v> {0.25:16.8f} </v>" for _ in KPOINTS)
    out.append("  </varray>")
    if labels:
        out += ["  <kpoints_labels>", '   <i name="GAMMA" type="int">       1 </i>', '   <i name="M" type="int">       2 </i>',
                '   <i name="M" type="int">       3 </i>', '   <i name="K" type="int">       4 </i>', "  </kpoints_labels>"]
    out += [" </kpoints>", " <parameters>", '  <separator name="electronic" >', '   <separator name="electronic spin" >',
            _param("ISPIN", f"     {ispin}", "int"), _param("LNONCOLLINEAR", logical(nonc), "logical"),
            _param("LSORBIT", logical(nonc), "logical"), "   </separator>", _param("NBANDS", f"     {NBANDS}", "int"),
            _param("NELM", "    60", "int"), "  </separator>", '  <separator name="ionic" >', _param("NSW", "     0", "int"),
            _param("IBRION", "    -1", "int"), _param("ISIF", "     2", "int"), "  </separator>", '  <separator name="dos" >', _param("LORBIT", f"     {11 if projected else 10}", "int"),
            "  </separator>", " </parameters>", " <atominfo>", f"  <atoms>       {len(SPECIES)} </atoms>",
            '  <array name="atoms" >', '   <dimension dim="1">ion</dimension>', '   <field type="string">element</field>',
            '   <field type="int">atomtype</field>', "   <set>"]
    out.extend(f"    <rc><c>{name:2s}</c><c>   {1 if name == 'Cr' else 2}</c></rc>" for name in SPECIES)
    out += ["   </set>", "  </array>", " </atominfo>"]
    out += _structure("initialpos", " ")
    out += [" <calculation>", "  <scstep>", '   <energy><i name="e_fr_energy">     -1.0 </i></energy>', "  </scstep>"]
    out += _structure(None, "  ")
    out += _eigen_block("  ", eig_blocks)
    out += ["  <dos>", f'   <i name="efermi">   {EFERMI:14.8f} </i>', "   <total>", "    <array>",
            '     <dimension dim="1">gridpoints</dimension>', '     <dimension dim="2">spin</dimension>',
            "     <field>energy</field>", "     <field>total</field>", "     <field>integrated</field>", "     <set>"]
    for s in range(len(COMPONENT_LABELS[layout])):
        out += [f'      <set comment="spin {s + 1}">', "       <r>   -1.0000   0.5000   0.0000 </r>",
                "       <r>    0.0000   0.7000   0.6000 </r>", "      </set>"]
    out += ["     </set>", "    </array>", "   </total>", "   <partial>", "    <array>",
            '     <dimension dim="1">gridpoints</dimension>', '     <dimension dim="2">spin</dimension>',
            '     <dimension dim="3">ion</dimension>', "     <field>energy</field>"]
    out.extend(f"     <field>{name:>5s}</field>" for name in orbitals)
    out.append("     <set>")
    for i in range(len(SPECIES)):
        out.append(f'      <set comment="ion {i + 1}">')
        out += ['       <set comment="spin 1">', "        <r> -1.0 " + " 0.9999" * len(orbitals) + " </r>", "       </set>"]
        out.append("      </set>")
    out += ["     </set>", "    </array>", "   </partial>", "  </dos>"]
    if projected:
        out += ["  <projected>"]
        out += _eigen_block("   ", eig_blocks)
        # Dimension names from the outside in; band_outside_kpoint swaps them to test that the parser follows <dimension>.
        outer, middle = ("band", "kpoint") if band_outside_kpoint else ("kpoint", "band")
        counts = {"kpoint": len(KPOINTS), "band": NBANDS}
        out += ["   <array>", '    <dimension dim="1">ion</dimension>', f'    <dimension dim="2">{middle}</dimension>',
                f'    <dimension dim="3">{outer}</dimension>', '    <dimension dim="4">spin</dimension>']
        out.extend(f"    <field>{name:>5s}</field>" for name in orbitals)
        out.append("    <set>")
        for c in range(proj_blocks):
            out.append(f'     <set comment="spin{c + 1}">')
            for x in range(counts[outer]):
                out.append(f'      <set comment="{outer} {x + 1}">')
                for y in range(counts[middle]):
                    k, b = (y, x) if band_outside_kpoint else (x, y)
                    if drop_band_at == (k, b):
                        continue
                    out.append(f'       <set comment="{middle} {y + 1}">')
                    for i in range(len(SPECIES)):
                        row = " ".join(f"{weight(c, k, b, i, o, signed=nonc):7.4f}" for o in range(len(orbitals)))
                        out.append(f"        <r> {row} </r>")
                    out.append("       </set>")
                out.append("      </set>")
            out.append("     </set>")
        out += ["    </set>", "   </array>", "  </projected>"]
    out += [" </calculation>"]
    out += _structure("finalpos", " ")
    out.append("</modeling>")
    return "\n".join(out) + "\n"


# Text of an <eigenvalues> block with the known eigenvalues and occupations.
def _eigen_block(indent: str, blocks: int) -> list[str]:
    out = [f"{indent}<eigenvalues>", f"{indent} <array>", f'{indent}  <dimension dim="1">band</dimension>',
           f'{indent}  <dimension dim="2">kpoint</dimension>', f'{indent}  <dimension dim="3">spin</dimension>',
           f"{indent}  <field>eigene</field>", f"{indent}  <field>occ</field>", f"{indent}  <set>"]
    for s in range(blocks):
        out.append(f'{indent}   <set comment="spin {s + 1}">')
        for k in range(len(KPOINTS)):
            out.append(f'{indent}    <set comment="kpoint {k + 1}">')
            out.extend(f"{indent}     <r> {eigenvalue(s, k, b):9.4f} {occupation(b):9.4f} </r>" for b in range(NBANDS))
            out.append(f"{indent}    </set>")
        out.append(f"{indent}   </set>")
    out += [f"{indent}  </set>", f"{indent} </array>", f"{indent}</eigenvalues>"]
    return out


# Text of a <structure> block with the fixed lattice and positions.
def _structure(name: str | None, indent: str) -> list[str]:
    head = f'{indent}<structure name="{name}" >' if name else f"{indent}<structure>"
    out = [head, f"{indent} <crystal>", f'{indent}  <varray name="basis" >']
    out.extend(f"{indent}   <v> {a:16.8f} {b:16.8f} {c:16.8f} </v>" for a, b, c in LATTICE)
    out += [f"{indent}  </varray>", f'{indent}  <i name="volume">     66.18 </i>', f'{indent}  <varray name="rec_basis" >']
    out.extend(f"{indent}   <v> {a:16.8f} {b:16.8f} {c:16.8f} </v>" for a, b, c in RECIPROCAL)
    out += [f"{indent}  </varray>", f"{indent} </crystal>", f'{indent} <varray name="positions" >']
    out.extend(f"{indent}  <v> {a:16.8f} {b:16.8f} {c:16.8f} </v>" for a, b, c in POSITIONS)
    out += [f"{indent} </varray>", f"{indent}</structure>"]
    return out


# File name and builder arguments of every fixture written to tests/fixtures/orbital/.
FIXTURES = {
    "nonmagnetic_vasprun.xml": dict(layout="nonmagnetic"),
    "collinear_f_vasprun.xml": dict(layout="collinear", orbitals=ORBITALS_16),
    "soc_vasprun.xml": dict(layout="noncollinear", labels=True),
    "no_projected_vasprun.xml": dict(layout="nonmagnetic", projected=False),
    "missing_band_vasprun.xml": dict(layout="nonmagnetic", drop_band_at=(2, 1)),
    "extra_block_vasprun.xml": dict(layout="collinear", extra_projected_block=True),
    "permuted_dims_vasprun.xml": dict(layout="collinear", orbitals=ORBITALS_16, band_outside_kpoint=True),
}


# Text of every fixture, keyed by file name.
def all_fixtures() -> dict[str, str]:
    texts = {name: build(**args) for name, args in FIXTURES.items()}
    whole = texts["soc_vasprun.xml"]
    texts["truncated_vasprun.xml"] = whole[: whole.index('<set comment="spin3">')]
    return texts


# Write every fixture file into tests/fixtures/orbital/.
def write_all(out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, text in all_fixtures().items():
        (out_dir / name).write_text(text, encoding="latin-1")


# Keep the first nk k-points and nb bands of every band-data <array> below one element.
def _trim_band_arrays(parent: etree._Element, nk: int, nb: int) -> None:
    for array in parent.iter("array"):
        for spin in array.find("set").iterchildren("set"):
            for kpoint in list(spin.iterchildren("set"))[nk:]:
                spin.remove(kpoint)
            for kpoint in spin.iterchildren("set"):
                for band in list(kpoint)[nb:]:
                    kpoint.remove(band)


# Cut a real band run down to nk k-points and nb bands, keeping headers, structures, Fermi energy and all blocks consistent.
def trim_real_file(source: Path, target: Path, nk: int = 2, nb: int = 4) -> None:
    tree = etree.parse(str(source), etree.XMLParser(huge_tree=True))
    root = tree.getroot()
    for nbands in root.xpath("incar/i[@name='NBANDS'] | parameters//i[@name='NBANDS']"):
        nbands.text = f"{nb:6d}"
    for varray in root.xpath("kpoints/varray[@name='kpointlist' or @name='weights']"):
        for row in list(varray)[nk:]:
            varray.remove(row)
    for calculation in root.iterchildren("calculation"):
        for child in list(calculation):
            if child.tag not in ("structure", "eigenvalues", "dos", "projected"):
                calculation.remove(child)
        for dos in calculation.iterchildren("dos"):
            for child in list(dos):
                if child.get("name") != "efermi":
                    dos.remove(child)
        for block in calculation.xpath("eigenvalues | projected"):
            _trim_band_arrays(block, nk, nb)
    tree.write(str(target), xml_declaration=True, encoding="ISO-8859-1")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--trim":
        trim_real_file(Path(sys.argv[2]), OUT_DIR / TRIMMED_NAME)
    else:
        write_all(Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DIR)
