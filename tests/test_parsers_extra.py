"""Tests for the FINDSYM, wannier90.win and INCAR parsers on real files and small inline strings."""

import tempfile
import unittest
import os
from pathlib import Path

from altermagagent.parsers_extra import (diff_incar_dicts, incar_diff, parse_findsym_log, parse_wannier_win,
                                         parse_wannier_win_text, unicode_label, write_findsym_input)

DATA = Path(os.environ.get("ALTERMAG_DATA", Path.home() / "development/osi-hackathon/data"))
CRSE = DATA / "diana_CrSe/CrSe-test/coplanar"
MNTE = DATA / "diana_MnTe/hackaton_exa"

FINDSYM_TEXT = """FINDSYM, Version 4.2, June 2014
Written by Harold T. Stokes, Branton J. Campbell, and Dorian M. Hatch
Brigham Young University

MnTe
Tolerance:    0.00100
Lattice parameters, a,b,c,alpha,beta,gamma:
   4.14800   4.14800   6.71000  90.00000  90.00000 120.00000
Centering: P
------------------------------------------
Space Group 194.268 P6_3'/m'm'c
Number of operations: 24
"""


# Write text to a temporary file and return its path.
def _tmp(text, name):
    d = Path(tempfile.mkdtemp())
    p = d / name
    p.write_text(text)
    return p


class FindsymLog(unittest.TestCase):
    """Check the FINDSYM log parser and the findsym.in writer."""

    # The real CrSe log gives magnetic group 189.225 P-62'm' on space group 189 P-62m.
    @unittest.skipUnless((CRSE / "1-noSOC/findsym.log").exists(), "real CrSe data not on this machine")
    def test_real_crse(self):
        r = parse_findsym_log(CRSE / "1-noSOC/findsym.log")
        self.assertNotIn("error", r)
        self.assertEqual(r["bns_number"], "189.225")
        self.assertEqual(r["bns_label"], "P-62'm'")
        self.assertEqual(r["bns_label_unicode"], "P6̅2′m′")
        self.assertEqual(r["space_group_number"], 189)
        self.assertEqual(r["space_group_symbol"], "P-62m")
        self.assertTrue(r["magnetic"])
        self.assertEqual(r["num_atoms"], 12)
        self.assertEqual(r["lattice"], [6.36355, 6.36355, 5.995, 90.0, 90.0, 120.0])
        self.assertEqual(r["wyckoff"], ["i", "a", "c", "g"])
        self.assertIsInstance(r["notes"], list)

    # An inline log with a screw axis and an operation count is read and converted to Unicode.
    def test_inline(self):
        r = parse_findsym_log(_tmp(FINDSYM_TEXT, "findsym.log"))
        self.assertEqual(r["bns_label"], "P6_3'/m'm'c")
        self.assertEqual(r["bns_label_unicode"], "P6₃′/m′m′c")
        self.assertEqual(r["space_group_number"], 194)
        self.assertEqual(r["num_operations"], 24)
        self.assertEqual(r["tolerance"], 0.001)
        self.assertIsNone(r["lattice"])
        self.assertEqual(r["lattice_input"][0], 4.148)

    # A newer CIF-style output with OG tags and an operation loop is read.
    def test_cif_style(self):
        text = ("_space_group_magn.number_BNS 189.225\n_space_group_magn.name_BNS \"P-62'm'\"\n"
                "_space_group_magn.number_OG 189.5.1455\n_space_group_magn.name_OG \"P-62'm'\"\n"
                "loop_\n_space_group_symop_magn_operation.id\n_space_group_symop_magn_operation.xyz\n"
                "1 x,y,z,+1\n2 -y,x-y,z,+1\n3 -x+y,-x,z,+1\n\n")
        r = parse_findsym_log(_tmp(text, "findsym.cif"))
        self.assertEqual(r["bns_label"], "P-62'm'")
        self.assertEqual(r["og_number"], "189.5.1455")
        self.assertEqual(r["num_operations"], 3)

    # A missing file and an empty file never raise.
    def test_failures(self):
        self.assertIn("error", parse_findsym_log("/nonexistent/findsym.log"))
        r = parse_findsym_log(_tmp("", "findsym.log"))
        self.assertIsNone(r["space_group_number"])
        self.assertTrue(r["notes"])

    # Unicode labels handle bars, primes and subscripts.
    def test_unicode(self):
        self.assertEqual(unicode_label("P-62'm'"), "P6̅2′m′")
        self.assertEqual(unicode_label("P6_3/mmc"), "P6₃/mmc")
        self.assertEqual(unicode_label(None), "")

    # The writer copies the physicist's findsym.in layout and the result reads back in the same order.
    @unittest.skipUnless((CRSE / "1-noSOC/020_scf/POSCAR").exists(), "real CrSe data not on this machine")
    def test_write_input(self):
        moms = [(-1, -1, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0), (0, -1, 0), (-1, 0, 0)] + [(0, 0, 0)] * 6
        out = Path(tempfile.mkdtemp()) / "findsym.in"
        write_findsym_input(CRSE / "1-noSOC/020_scf/POSCAR", moms, out)
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "CrSe")
        self.assertEqual(lines[1:3], ["0.0001", "2"])
        self.assertEqual([float(x) for x in lines[3].split()][3:], [90.0, 90.0, 120.0])
        self.assertEqual(lines[4:7], ["2", "P", "12"])
        self.assertEqual(lines[7].split(), ["Cr"] * 6 + ["Se"] * 6)
        self.assertEqual(lines[8], "magnetic")
        self.assertEqual(len(lines), 21)
        self.assertEqual([float(x) for x in lines[9].split()[3:]], [-1.0, -1.0, 0.0])
        with self.assertRaises(ValueError):
            write_findsym_input(CRSE / "1-noSOC/020_scf/POSCAR", moms[:3], out)


class WannierWin(unittest.TestCase):
    """Check the wannier90.win parser."""

    # The SOC CrSe file gives counts, windows, four projection lines and the flags.
    @unittest.skipUnless((CRSE / "2-SOC/3.1-wannier/wannier90.win").exists(), "real CrSe data not on this machine")
    def test_real_crse_soc(self):
        r = parse_wannier_win(CRSE / "2-SOC/3.1-wannier/wannier90.win")
        self.assertEqual((r["num_wann"], r["num_bands"]), (96, 108))
        self.assertEqual((r["dis_win_min"], r["dis_win_max"]), (-1.0, 10.5))
        self.assertEqual(r["dis_froz_max"], 8.33728721)
        self.assertEqual(r["projections"], ["Cr:d(u)", "Se:p(u)", "Cr:d(d)", "Se:p(d)"])
        self.assertTrue(r["spinors"] and r["write_hr"] and r["bands_plot"] and r["kpoint_path"])
        self.assertTrue(r["summary"].startswith("num_wann 96, bands 108, outer window -1 to 10.5 eV"))

    # The no-SOC CrSe file has two projection lines.
    @unittest.skipUnless((CRSE / "1-noSOC/04_wannier/wannier90.win").exists(), "real CrSe data not on this machine")
    def test_real_crse_nosoc(self):
        r = parse_wannier_win(CRSE / "1-noSOC/04_wannier/wannier90.win")
        self.assertEqual(r["projections"], ["Cr:d", "Se:p"])
        self.assertIn("projections: Cr:d, Se:p", r["summary"])

    # The MnTe file uses .true. values, no band path and windows written first.
    @unittest.skipUnless((MNTE / "9.6-wannier/wannier90.win").exists(), "real MnTe data not on this machine")
    def test_real_mnte(self):
        r = parse_wannier_win(MNTE / "9.6-wannier/wannier90.win")
        self.assertEqual((r["num_wann"], r["num_bands"]), (36, 60))
        self.assertEqual((r["dis_win_max"], r["dis_froz_max"]), (20.0, 8.1))
        self.assertTrue(r["spinors"] and r["write_hr"])
        self.assertFalse(r["bands_plot"] or r["kpoint_path"])
        self.assertEqual(r["projections"][0], "Mn: s;d(u)")

    # Uppercase keys, ":" separators, "!" comments and missing windows are handled.
    def test_inline(self):
        r = parse_wannier_win_text("NUM_WANN : 8\nNum_Bands 12 ! comment\nSPINORS = F\n"
                                   "Begin Projections\nGa:sp3\nEnd Projections\nBANDS_PLOT : T\n")
        self.assertEqual((r["num_wann"], r["num_bands"]), (8, 12))
        self.assertFalse(r["spinors"])
        self.assertTrue(r["bands_plot"])
        self.assertEqual(r["projections"], ["Ga:sp3"])
        self.assertEqual(r["summary"], "num_wann 8, bands 12, outer window: no window set, "
                                       "frozen: no window set, projections: Ga:sp3, flags: bands_plot")
        self.assertIn("error", parse_wannier_win("/nonexistent/wannier90.win"))


class IncarDiff(unittest.TestCase):
    """Check the INCAR difference against a parent run."""

    # The no-SOC Wannier interface run differs from the first SCF in ICHARG, NBANDS and the Wannier tags.
    @unittest.skipUnless((CRSE / "1-noSOC/04_wannier/INCAR").exists(), "real CrSe data not on this machine")
    def test_real_crse(self):
        r = incar_diff(CRSE / "1-noSOC/04_wannier/INCAR", CRSE / "1-noSOC/020_scf/INCAR")
        self.assertEqual(r["changed"]["ICHARG"], (2, 11))
        self.assertEqual(r["added"]["NBANDS"], 120)
        self.assertNotIn("MAGMOM", r["changed"])
        self.assertIn("M_CONSTR changed", r["summary"])
        self.assertIn("ICHARG 2 -> 11", r["summary"])

    # The SOC SCF differs from the no-SOC SCF with 120 bands by LSORBIT and SAXIS.
    @unittest.skipUnless((CRSE / "2-SOC/1-scf/INCAR").exists(), "real CrSe data not on this machine")
    def test_real_soc(self):
        r = incar_diff(CRSE / "2-SOC/1-scf/INCAR", CRSE / "1-noSOC/023_scf_nbands120/INCAR")
        self.assertEqual(r["changed"]["LSORBIT"], (False, True))
        self.assertIn("LSORBIT .FALSE. -> .TRUE.", r["summary"])

    # Logical spellings and MAGMOM rounding within 0.01 count as equal; a real MAGMOM change is named only.
    def test_inline(self):
        parent = {"LSORBIT": "F", "MAGMOM": [3.0, -3.0, 0.0], "ICHARG": 2, "ALGO": "Normal", "NELM": 60}
        this = {"lsorbit": ".TRUE.", "MAGMOM": [3.004, -3.0, 0.0], "ICHARG": 11, "ALGO": "normal", "NBANDS": 120}
        r = diff_incar_dicts(this, parent)
        self.assertEqual(r["summary"], "ICHARG 2 -> 11, LSORBIT .FALSE. -> .TRUE., NBANDS added 120, NELM removed")
        r = diff_incar_dicts({"MAGMOM": [3.0, 3.0], "LWAVE": "T"}, {"MAGMOM": [3.0, -3.0], "LWAVE": True})
        self.assertEqual(r["summary"], "MAGMOM changed")
        self.assertEqual(diff_incar_dicts({"A": 1}, {"A": 1})["summary"], "no changes")
        many = diff_incar_dicts({f"T{i}": i for i in range(8)}, {})
        self.assertTrue(many["summary"].endswith("and 2 more"))
        self.assertIn("error", incar_diff("/nonexistent/INCAR", "/nonexistent/INCAR"))


if __name__ == "__main__":
    unittest.main()
