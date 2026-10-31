"""Stage 2 of the orbital plotter: k-path distances, segment breaks, label sources, and a check against VASPKIT output."""

import dataclasses
import math
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from vasp_orbital_plotter.kpath import (KPathError, build_kpath, corner_to_index, cumulative_distance,
                                        labels_from_vasprun, parse_label_option, read_kpoints_line_mode)
from vasp_orbital_plotter.parser import parse_vasprun

from _fixtures import tmp_folder

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "orbital"
REAL_FOLDER = Path("/home/clawc/development/osi-hackathon/data/diana_MnTe/hackaton_exa/9.2-band")
# Hand value: 2 pi * 0.25 * |b1| with |b1| = |b2| = 0.32075015 1/Angstrom of the fixture's hexagonal rec_basis.
QUARTER_STEP = 0.5038331560045315
TOL = 1e-6

KPOINTS_MATCHING = """path matching the fixture
2
Line-mode
Reciprocal
 0.0 0.0 0.0 ! GAMMA
 0.25 0.0 0.0 ! A

 0.5 0.0 0.0 ! M
 0.5 0.25 0.0 ! K
"""


# Copy one fixture into a fresh temporary folder as vasprun.xml, with optional sibling files; return the vasprun path.
def copy_with_siblings(test: unittest.TestCase, fixture: str, siblings: dict[str, str]) -> Path:
    folder = tmp_folder(test, "orbital_kpath_")
    shutil.copy(FIXTURES / fixture, folder / "vasprun.xml")
    for name, text in siblings.items():
        (folder / name).write_text(text)
    return folder / "vasprun.xml"


class Distances(unittest.TestCase):
    """Distances use 2 pi times VASP's rec_basis, in 1/Angstrom."""

    # Parse the nonmagnetic fixture once.
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = parse_vasprun(FIXTURES / "nonmagnetic_vasprun.xml", progress=False)

    # Without breaks, the four fixture k-points are three quarter steps of b1, b1 and b2 apart.
    def test_continuous_distance(self) -> None:
        x = cumulative_distance(self.data.kpoints, self.data.reciprocal_lattice)
        np.testing.assert_allclose(x, [0.0, QUARTER_STEP, 2 * QUARTER_STEP, 3 * QUARTER_STEP], atol=TOL)

    # A break before index 2 adds no distance between k-points 2 and 3 (1-based).
    def test_break_adds_nothing(self) -> None:
        x = cumulative_distance(self.data.kpoints, self.data.reciprocal_lattice, breaks=[2])
        np.testing.assert_allclose(x, [0.0, QUARTER_STEP, QUARTER_STEP, 2 * QUARTER_STEP], atol=TOL)

    # Break index 0 (and any index below 1) is ignored instead of zeroing the last step.
    def test_break_zero_ignored(self) -> None:
        full = cumulative_distance(self.data.kpoints, self.data.reciprocal_lattice)
        for breaks in ([0], [-1, 0]):
            with self.subTest(breaks=breaks):
                np.testing.assert_allclose(cumulative_distance(self.data.kpoints, self.data.reciprocal_lattice, breaks), full)

    # The 2 pi convention: a full reciprocal vector b1 has length 2 pi |rec_basis row 1|.
    def test_two_pi_convention(self) -> None:
        x = cumulative_distance(np.array([[0, 0, 0], [1, 0, 0]]), self.data.reciprocal_lattice)
        self.assertAlmostEqual(x[1], 2 * math.pi * float(np.linalg.norm(self.data.reciprocal_lattice[0])), places=9)
        self.assertAlmostEqual(x[1], 4 * QUARTER_STEP, places=6)

    # The fixture's listgenerated divisions (2) split it into two segments; the jump from 0.25 to 0.5 is a break.
    def test_vasprun_divisions(self) -> None:
        kpath = build_kpath(self.data, folder=Path(tempfile.gettempdir()) / "no_such_orbital_folder")
        self.assertEqual(kpath.segment_source, "vasprun line-mode divisions")
        self.assertEqual((kpath.segment_starts, kpath.branches), ([2], [(0, 1), (2, 3)]))
        np.testing.assert_allclose(kpath.x, [0.0, QUARTER_STEP, QUARTER_STEP, 2 * QUARTER_STEP], atol=TOL)

    # Without divisions, a repeated k-point ends a segment but is no break, since both ends are the same point.
    def test_repeated_point(self) -> None:
        k = np.array([[0, 0, 0], [0.25, 0, 0], [0.25, 0, 0], [0.5, 0, 0]], dtype=float)
        data = dataclasses.replace(self.data, kpoints=k, kpoint_generation=None, kpoint_divisions=None)
        kpath = build_kpath(data, folder=Path(tempfile.gettempdir()) / "no_such_orbital_folder")
        self.assertEqual((kpath.segment_starts, kpath.segment_source), ([2], "repeated k-points"))
        np.testing.assert_allclose(kpath.x, [0.0, QUARTER_STEP, QUARTER_STEP, 2 * QUARTER_STEP], atol=TOL)

    # Without any segment structure, even a large jump stays one continuous line, and a note says so.
    def test_no_structure_is_one_line(self) -> None:
        k = np.array([[0, 0, 0], [0.01, 0, 0], [0.02, 0, 0], [0.03, 0, 0], [0.5, 0.25, 0], [0.5, 0.26, 0]], dtype=float)
        data = dataclasses.replace(self.data, kpoints=k, kpoint_generation=None, kpoint_divisions=None)
        kpath = build_kpath(data, folder=Path(tempfile.gettempdir()) / "no_such_orbital_folder")
        self.assertEqual((kpath.segment_starts, kpath.branches), ([], [(0, 5)]))
        self.assertEqual(kpath.segment_source, "none (one continuous segment)")
        self.assertGreater(kpath.x[4], kpath.x[3])
        self.assertTrue(any("drawn as one line" in note for note in kpath.notes))

    # A matching line-mode KPOINTS file is the first segment source and gives labels when nothing else does.
    def test_kpoints_file(self) -> None:
        vasprun = copy_with_siblings(self, "nonmagnetic_vasprun.xml", {"KPOINTS": KPOINTS_MATCHING})
        kpath = build_kpath(parse_vasprun(vasprun, progress=False))
        self.assertEqual((kpath.segment_source, kpath.label_source), ("KPOINTS", "KPOINTS comments"))
        self.assertEqual([label for _, label in kpath.ticks], ["Γ", "A|M", "K"])

    # A KPOINTS file whose corners do not match the k-points is ignored with a note.
    def test_kpoints_mismatch(self) -> None:
        vasprun = copy_with_siblings(self, "nonmagnetic_vasprun.xml", {"KPOINTS": KPOINTS_MATCHING.replace("0.5 0.25", "0.5 0.4")})
        kpath = build_kpath(parse_vasprun(vasprun, progress=False))
        self.assertEqual(kpath.segment_source, "vasprun line-mode divisions")
        self.assertTrue(any("do not match" in note for note in kpath.notes))

    # The KPOINTS reader accepts labels after "!" or as a bare fourth column, and refuses automatic meshes.
    def test_kpoints_reader(self) -> None:
        folder = tmp_folder(self, "orbital_kpts_")
        (folder / "KPOINTS").write_text("x\n10\nLine\nRec\n0 0 0 GAMMA\n0.5 0 0 X\n")
        parsed, note = read_kpoints_line_mode(folder / "KPOINTS")
        self.assertIsNone(note)
        self.assertEqual((parsed.divisions, parsed.nsegments, [c[1] for c in parsed.corners]), (10, 1, ["GAMMA", "X"]))
        (folder / "KPOINTS").write_text("x\n0\nGamma\n4 4 4\n")
        parsed, note = read_kpoints_line_mode(folder / "KPOINTS")
        self.assertIsNone(parsed)
        self.assertIn("not in line mode", note)


class Labels(unittest.TestCase):
    """Label sources in order: --high-symmetry, KLABELS, <kpoints_labels>, KPOINTS comments, none."""

    # The option uses 1-based k-point numbers and converts them to 0-based indices.
    def test_label_option(self) -> None:
        self.assertEqual(parse_label_option("1:G, 41:M,81:K,121:G", 121), [(0, "G"), (40, "M"), (80, "K"), (120, "G")])
        for bad in ("0:G", "122:G", "G", "1-G", ""):
            with self.subTest(bad=bad), self.assertRaises(KPathError):
                parse_label_option(bad, 121)

    # Corner line n of a line-mode path with N points per segment is k-point (n // 2) * N or that plus N - 1.
    def test_corner_to_index(self) -> None:
        self.assertEqual([corner_to_index(n, 20) for n in range(1, 7)], [0, 19, 20, 39, 40, 59])

    # The SOC fixture's <kpoints_labels> (corners 1..4, 2 points per segment) map onto k-points 1..4.
    def test_vasprun_labels(self) -> None:
        data = parse_vasprun(FIXTURES / "soc_vasprun.xml", progress=False)
        self.assertEqual(labels_from_vasprun(data), [(0, "GAMMA"), (1, "M"), (2, "M"), (3, "K")])
        kpath = build_kpath(data, folder=Path(tempfile.gettempdir()) / "no_such_orbital_folder")
        self.assertEqual(kpath.label_source, "vasprun <kpoints_labels>")
        self.assertEqual([label for _, label in kpath.ticks], ["Γ", "M", "K"])

    # A <kpoints_labels> block that lists only the labelled corners (here without the last one) is still used.
    def test_vasprun_labels_missing_last_corner(self) -> None:
        data = parse_vasprun(FIXTURES / "soc_vasprun.xml", progress=False)
        partial = dataclasses.replace(data, kpoint_label_entries=data.kpoint_label_entries[:-1])
        self.assertEqual(labels_from_vasprun(partial), [(0, "GAMMA"), (1, "M"), (2, "M")])
        self.assertIsNone(labels_from_vasprun(dataclasses.replace(data, kpoint_label_entries=[(1, "G"), (5, "X")])))
        self.assertIsNone(labels_from_vasprun(dataclasses.replace(data, kpoint_divisions=3)))

    # --high-symmetry wins over KLABELS, KLABELS over <kpoints_labels>, and <kpoints_labels> over KPOINTS comments.
    def test_source_order(self) -> None:
        klabels = "K-Label  Coordinate\n  GAMMA 0.000\n  Y 0.504\n  Z 1.008\n \n* Tip\n"
        vasprun = copy_with_siblings(self, "soc_vasprun.xml", {"KPOINTS": KPOINTS_MATCHING, "KLABELS": klabels})
        data = parse_vasprun(vasprun, progress=False)
        self.assertEqual(build_kpath(data, labels="1:G,4:X").label_source, "--high-symmetry option")
        self.assertEqual([lab for _, lab in build_kpath(data, labels="1:G,4:X").ticks], ["Γ", "X"])
        kpath = build_kpath(data)
        self.assertEqual(kpath.label_source, "KLABELS")
        self.assertEqual([lab for _, lab in kpath.ticks], ["Γ", "Y", "Z"])
        (vasprun.parent / "KLABELS").unlink()
        self.assertEqual(build_kpath(parse_vasprun(vasprun, progress=False)).label_source, "vasprun <kpoints_labels>")

    # KLABELS written without the factor 2 pi is rescaled to the plotted length, with a note.
    def test_klabels_rescaled(self) -> None:
        length = 2 * QUARTER_STEP
        klabels = f"K-Label  Coordinate\n  GAMMA 0.000\n  M {QUARTER_STEP / (2 * math.pi):.5f}\n  K {length / (2 * math.pi):.5f}\n"
        vasprun = copy_with_siblings(self, "nonmagnetic_vasprun.xml", {"KLABELS": klabels})
        kpath = build_kpath(parse_vasprun(vasprun, progress=False))
        self.assertEqual(kpath.label_source, "KLABELS")
        np.testing.assert_allclose([p for p, _ in kpath.ticks], [0.0, QUARTER_STEP, length], atol=1e-9)
        self.assertTrue(any("rescaled" in note for note in kpath.notes))

    # With no label source, ticks sit at the segment ends without names.
    def test_no_labels(self) -> None:
        data = parse_vasprun(FIXTURES / "nonmagnetic_vasprun.xml", progress=False)
        kpath = build_kpath(data, folder=Path(tempfile.gettempdir()) / "no_such_orbital_folder")
        self.assertEqual(kpath.label_source, "none")
        self.assertTrue(all(label == "" for _, label in kpath.ticks))


# Read VASPKIT's BAND.dat into an array (band, kpoint, [distance, energy]); odd bands run backwards and are reversed.
def read_band_dat(path: Path) -> np.ndarray:
    blocks: list[list[list[float]]] = []
    for line in path.read_text().splitlines():
        if line.startswith("# Band-Index"):
            blocks.append([])
        elif blocks and line.strip() and not line.startswith("#"):
            blocks[-1].append([float(v) for v in line.split()])
    array = np.array(blocks)
    array[1::2] = array[1::2, ::-1]
    return array


@unittest.skipUnless((REAL_FOLDER / "vasprun.xml").is_file(), "Diana's MnTe 9.2-band run is not on this machine")
class RealAgainstVaspkit(unittest.TestCase):
    """MnTe 9.2-band, read-only: energies and distances against VASPKIT's BAND.dat, REFORMATTED_BAND.dat and KLABELS."""

    # Parse the real run once and build its k-path from the files next to it.
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = parse_vasprun(REAL_FOLDER / "vasprun.xml", progress=False)
        cls.kpath = build_kpath(cls.data)

    # Segments come from the KPOINTS file, labels from KLABELS; the path is continuous with five segments of 20 points.
    def test_sources(self) -> None:
        self.assertEqual((self.kpath.segment_source, self.kpath.label_source), ("KPOINTS", "KLABELS"))
        self.assertEqual(self.kpath.segment_starts, [20, 40, 60, 80])
        self.assertEqual([lab for _, lab in self.kpath.ticks], ["Γ", "Y", "C₀", "Σ₀", "Γ", "S"])

    # BAND.dat energies (E - E_F, 6 decimals from EIGENVAL) agree within the 4-decimal precision of vasprun.xml.
    def test_band_dat_energies_and_distances(self) -> None:
        reference = read_band_dat(REAL_FOLDER / "BAND.dat")
        self.assertEqual(reference.shape, (60, 100, 2))
        energies = self.data.eigenvalues[0] - self.data.efermi
        self.assertLess(float(np.abs(reference[:, :, 1].T - energies).max()), 5.1e-5 + 2e-7)
        # BAND.dat distances (5 decimals) equal ours, so VASPKIT also uses 1/Angstrom with the factor 2 pi.
        self.assertLess(float(np.abs(reference[:, :, 0].T - self.kpath.x[:, None]).max()), 5.1e-6)

    # REFORMATTED_BAND.dat has 4 decimals in both columns; agreement is within its rounding plus vasprun's.
    def test_reformatted_band(self) -> None:
        reference = np.loadtxt(REAL_FOLDER / "REFORMATTED_BAND.dat")
        energies = self.data.eigenvalues[0] - self.data.efermi
        self.assertLess(float(np.abs(reference[:, 1:] - energies).max()), 1.0e-4 + 2e-7)
        self.assertLess(float(np.abs(reference[:, 0] - self.kpath.x).max()), 5.1e-5)

    # KLABELS distances (3 decimals) sit at our segment ends.
    def test_klabels_positions(self) -> None:
        rows = [line.split() for line in (REAL_FOLDER / "KLABELS").read_text().splitlines()[1:7]]
        expected = [float(row[-1]) for row in rows]
        ends = [self.kpath.x[i] for i in (0, 19, 39, 59, 79, 99)]
        np.testing.assert_allclose(ends, expected, atol=5.1e-4)


if __name__ == "__main__":
    unittest.main()
