"""Symmetry labels on the synthetic CrSe demo and the mysymmop.dat parser; every label stays "parsed, pending validation"."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SYMMOP = """  ---------------    1 ---------------
  rotation:
       1  0  0
       0  1  0
       0  0  1
  translation:
       0.00000  0.00000  0.00000
  rotmap:
              1 -> 1

  ---------------    2 ---------------
  rotation:
      -1  0  0
       0 -1  0
       0  0  1
  translation:
       0.00000  0.00000  0.50000
  rotmap:
              1 -> 2
"""


class Symmetry(unittest.TestCase):
    """Magnetic and spin groups of the collinear demo run, and the wannhr_symm operation parser."""

    # Write the demo tree into a temporary folder.
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.base = Path(cls.tmp.name)
        subprocess.run([sys.executable, str(REPO / "demo" / "make_demo.py"), "--out", str(cls.base / "full")],
                       check=True, capture_output=True)

    # Remove the temporary folder.
    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    # The collinear demo run gets a magnetic space group whose operations split into unitary and antiunitary halves.
    def test_magnetic_group_collinear(self):
        from altermagagent import symmetry
        m = symmetry.magnetic_space_group(self.base / "full" / "scf_nosoc")
        self.assertEqual(m["state"], symmetry.PENDING)
        self.assertIn(m["type"], (1, 2, 3, 4))
        self.assertEqual(m["n_unitary"] + m["n_antiunitary"], m["n_operations"])
        self.assertEqual(m["space_group"]["number"], 194)

    # The spin group of the collinear run is either computed with a collinear spin-only group or says "not computed".
    def test_spin_group_collinear(self):
        from altermagagent import symmetry
        s = symmetry.spin_space_group(self.base / "full" / "scf_nosoc")
        self.assertIn(s["state"], (symmetry.PENDING, "not computed"))
        if s["state"] == symmetry.PENDING:
            self.assertEqual(s["spin_only_type"], "collinear")

    # mysymmop.dat blocks parse into rotation, translation and words such as "2-fold screw".
    def test_parse_mysymmop(self):
        from altermagagent import symmetry
        p = self.base / "mysymmop.dat"
        p.write_text(SYMMOP)
        ops = symmetry.parse_mysymmop(p)
        self.assertEqual(len(ops), 2)
        self.assertEqual(ops[1]["translation"], [0.0, 0.0, 0.5])
        self.assertEqual(symmetry.describe(ops[1]["rotation"], ops[1]["translation"]), "2-fold screw")


if __name__ == "__main__":
    unittest.main()
