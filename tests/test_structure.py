"""MAGMOM parsing and structure figures with surplus input values."""

import json
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from altermagagent.structure import plot_structure, read_moments


class StructureMoments(unittest.TestCase):
    """Use only the initial moments belonging to the atoms in the structure."""

    # Keep all synthetic inputs and figure outputs inside a temporary directory.
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.run = self.base / "run"
        self.run.mkdir()

    # Write a synthetic INCAR and read two atoms' initial magnetic moments.
    def moments(self, text):
        (self.run / "INCAR").write_text(text, encoding="utf-8")
        notes = []
        return read_moments(self.run, 2, notes), notes

    # Full and partial surplus vectors must leave the first two vectors unchanged.
    def test_noncollinear_surplus(self):
        for flag in ("LNONCOLLINEAR", "LSORBIT"):
            for extra in ("", " 99", " 99 98", " 99 98 97", " 2*99 4*98"):
                with self.subTest(flag=flag, extra=extra):
                    (vecs, raw, _, _), _ = self.moments(
                        f"{flag} = .TRUE.\nMAGMOM = 0 0 2 0 0 -2{extra}\n")
                    np.testing.assert_allclose(vecs, [[0, 0, 2], [0, 0, -2]])
                    self.assertEqual(raw, [[0, 0, 2], [0, 0, -2]])

    # Collinear calculations consume one scalar per atom regardless of the extra count.
    def test_collinear_surplus(self):
        for extra in ("", " 99", " 3*99"):
            with self.subTest(extra=extra):
                (vecs, raw, _, _), _ = self.moments(f"ISPIN = 2\nMAGMOM = 2 -2{extra}\n")
                np.testing.assert_allclose(vecs, [[0, 0, 2], [0, 0, -2]])
                self.assertEqual(raw, [2, -2])

    # Expand repeats across continued lines before consuming one vector per atom.
    def test_noncollinear_repeats(self):
        (vecs, raw, _, _), _ = self.moments(
            "LNONCOLLINEAR = .TRUE.\nMAGMOM = 2*0 2 \\\n  2*0 -2 6*99 ! ignored tail\n")
        np.testing.assert_allclose(vecs, [[0, 0, 2], [0, 0, -2]])
        self.assertEqual(len(raw), 2)

    # Flat values from INCAR readers must also be truncated before grouping into triples.
    def test_flat_noncollinear_surplus(self):
        (self.run / "INCAR").touch()
        for extra in ([], [99], [99, 98], [99, 98, 97]):
            with self.subTest(extra=extra):
                incar = {"LNONCOLLINEAR": True, "MAGMOM": [0, 0, 2, 0, 0, -2] + extra}
                with patch("pymatgen.io.vasp.Incar.from_file", return_value=incar):
                    vecs, raw, _, _ = read_moments(self.run, 2, [])
                np.testing.assert_allclose(vecs, [[0, 0, 2], [0, 0, -2]])
                self.assertEqual(raw, [[0, 0, 2], [0, 0, -2]])

    # Incomplete inputs still explain why no moment arrows can be drawn.
    def test_too_few_values(self):
        for text in ("ISPIN = 2\nMAGMOM = 2\n",
                     "LNONCOLLINEAR = .TRUE.\nMAGMOM = 0 0 2 0 0\n"):
            with self.subTest(text=text):
                (vecs, _, _, _), notes = self.moments(text)
                self.assertIsNone(vecs)
                self.assertTrue(any("no arrows drawn" in note for note in notes))

    # Surplus values must not affect rotation from the specified spin axis.
    def test_surplus_with_saxis(self):
        (vecs, raw, _, _), _ = self.moments(
            "LSORBIT = .TRUE.\nSAXIS = 1 0 0\nMAGMOM = 0 0 2 0 0 -2 3*99\n")
        np.testing.assert_allclose(vecs, [[2, 0, 0], [-2, 0, 0]], atol=1e-12)
        self.assertEqual(raw, [[0, 0, 2], [0, 0, -2]])

    # Real figure generation must retain arrows and exclude surplus values from the saved moments.
    def test_plot_with_surplus(self):
        (self.run / "POSCAR").write_text(
            "Two Fe atoms\n1.0\n3 0 0\n0 3 0\n0 0 3\nFe\n2\nDirect\n0 0 0\n0.5 0.5 0.5\n",
            encoding="utf-8")
        self.moments("LNONCOLLINEAR = .TRUE.\nMAGMOM = 0 0 2 0 0 -2 3*99\n")
        result = plot_structure(self.run, out_dir=self.base / "figures", supercell=(2, 1, 1))
        self.assertEqual([m["vector"] for m in result["moments"]], [[0, 0, 2], [0, 0, -2]])
        self.assertEqual(result["net_moment"], [0, 0, 0])
        for path in result["paths"].values():
            self.assertGreater(Path(path).stat().st_size, 0)
        html = Path(result["paths"]["html"]).read_text(encoding="utf-8")
        data = json.loads(re.search(r"const data = (.*);", html).group(1))
        self.assertEqual(len(data["arrows"]), 4)
        saved = json.loads(Path(result["paths"]["png"]).with_suffix(".json").read_text(encoding="utf-8"))
        self.assertEqual(saved["moments"], result["moments"])


if __name__ == "__main__":
    unittest.main()
