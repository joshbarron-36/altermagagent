"""Stage 1 of the orbital plotter: the vasprun.xml parser on synthetic fixtures and, when present, on Diana's MnTe band run."""

import contextlib
import gzip
import io
import re
import shutil
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "fixtures"))

import orbital_vasprun_builder as fx  # noqa: E402
from vasp_orbital_plotter import cli  # noqa: E402
from vasp_orbital_plotter.parser import (LAYOUT_COLLINEAR, LAYOUT_NONCOLLINEAR, LAYOUT_NONMAGNETIC,  # noqa: E402
                                         ProjectionsMissingError, VasprunLayoutError, VasprunParseError, parse_vasprun)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "orbital"
REAL_BAND_RUN = Path("/home/clawc/development/osi-hackathon/data/diana_MnTe/hackaton_exa/9.2-band/vasprun.xml")
TOL = 1e-6
MNTE_EFERMI = 5.59960418


# Expected projection array built from the fixture formula, shape (component, kpoint, band, ion, orbital).
def expected_weights(components: int, norb: int, signed: bool) -> np.ndarray:
    nk, nb, nion = len(fx.KPOINTS), fx.NBANDS, len(fx.SPECIES)
    out = np.zeros((components, nk, nb, nion, norb))
    for index in np.ndindex(out.shape):
        out[index] = fx.weight(*index, signed=signed)
    return out


# Expected eigenvalues from the fixture formula, shape (spin, kpoint, band).
def expected_eigenvalues(spins: int) -> np.ndarray:
    return np.array([[[fx.eigenvalue(s, k, b) for b in range(fx.NBANDS)] for k in range(len(fx.KPOINTS))] for s in range(spins)])


class FixtureFilesTest(unittest.TestCase):
    """The stored fixture files are exactly what the builder writes, so formulas and files cannot drift apart."""

    # Each stored fixture equals the builder's current output.
    def test_files_match_builder(self):
        for name, text in fx.all_fixtures().items():
            with self.subTest(name=name):
                self.assertEqual((FIXTURES / name).read_text(encoding="latin-1"), text)


class NonmagneticTest(unittest.TestCase):
    """Nonmagnetic layout: one eigenvalue spin block and one projection component."""

    # Parse the nonmagnetic fixture once for all tests.
    @classmethod
    def setUpClass(cls):
        cls.data = parse_vasprun(FIXTURES / "nonmagnetic_vasprun.xml")

    # Layout, component names and array shapes follow the documented convention.
    def test_layout_and_shapes(self):
        d = self.data
        self.assertEqual(d.layout, LAYOUT_NONMAGNETIC)
        self.assertEqual(d.projection_components, ("total",))
        self.assertEqual(d.eigenvalues.shape, (1, 4, 3))
        self.assertEqual(d.occupations.shape, (1, 4, 3))
        self.assertEqual(d.projections.shape, (1, 4, 3, 3, 9))

    # Eigenvalues, occupations, Fermi energy, k-points and weights are read exactly.
    def test_values(self):
        d = self.data
        np.testing.assert_allclose(d.eigenvalues, expected_eigenvalues(1), atol=TOL)
        np.testing.assert_allclose(d.occupations[0, :, 0], 1.0)
        np.testing.assert_allclose(d.occupations[0, :, 1:], 0.0)
        self.assertAlmostEqual(d.efermi, fx.EFERMI, places=8)
        np.testing.assert_allclose(d.kpoints, np.array(fx.KPOINTS), atol=TOL)
        np.testing.assert_allclose(d.kpoint_weights, 0.25)
        np.testing.assert_allclose(d.projections, expected_weights(1, 9, signed=False), atol=TOL)

    # Structure: species per ion, lattice, reciprocal lattice and positions, taken from finalpos.
    def test_structure(self):
        d = self.data
        self.assertEqual(d.species, fx.SPECIES)
        self.assertEqual(d.structure_source, "finalpos")
        np.testing.assert_allclose(d.lattice, np.array(fx.LATTICE), atol=TOL)
        np.testing.assert_allclose(d.reciprocal_lattice, np.array(fx.RECIPROCAL), atol=TOL)
        np.testing.assert_allclose(d.positions, np.array(fx.POSITIONS), atol=TOL)

    # VASP ion n (1-based) is array index n - 1; out-of-range and non-integer numbers are refused.
    def test_atom_indexing(self):
        d = self.data
        self.assertEqual(d.ion_index(1), 0)
        self.assertEqual(d.ion_index(3), 2)
        for bad in (0, 4, -1):
            with self.assertRaises(IndexError):
                d.ion_index(bad)
        with self.assertRaises(TypeError):
            d.ion_index(1.0)
        with self.assertRaises(TypeError):
            d.ion_index(True)
        self.assertEqual(d.ions_of_element("Cr"), [1, 3])
        self.assertEqual(d.ions_of_element("Se"), [2])
        with self.assertRaises(KeyError):
            d.ions_of_element("Te")
        # Ion 3 (Cr), orbital dxy, k-point 2, band 3 decodes to the formula value for indices (0, 1, 2, 2, 4).
        value = d.projections[0, 1, 2, d.ion_index(3), d.orbital_index("dxy")]
        self.assertAlmostEqual(value, fx.weight(0, 1, 2, 2, 4), places=8)

    # Orbital labels come from the <field> entries in file order, with padding stripped.
    def test_orbital_indexing(self):
        d = self.data
        self.assertEqual(d.orbital_labels, fx.ORBITALS_9)
        self.assertEqual(d.orbital_index("s"), 0)
        self.assertEqual(d.orbital_index("px"), 3)
        self.assertEqual(d.orbital_index("x2-y2"), 8)
        with self.assertRaises(KeyError):
            d.orbital_index("dx2-y2")

    # Charge weights are nonnegative and at most one, and the partial DOS is recorded but not used.
    def test_weight_range(self):
        self.assertGreaterEqual(self.data.projections.min(), 0.0)
        self.assertLessEqual(self.data.projections.max(), 1.0)
        self.assertTrue(self.data.has_partial_dos)
        self.assertLess(self.data.projections.max(), 0.9)


class CollinearTest(unittest.TestCase):
    """ISPIN = 2 layout: two eigenvalue spin blocks and two projection components, with f orbitals in the fields."""

    # Parse the collinear fixture once for all tests.
    @classmethod
    def setUpClass(cls):
        cls.data = parse_vasprun(FIXTURES / "collinear_f_vasprun.xml")

    # Two spin blocks, components up and down, 16 orbital columns.
    def test_layout_and_shapes(self):
        d = self.data
        self.assertEqual(d.layout, LAYOUT_COLLINEAR)
        self.assertEqual(d.projection_components, ("up", "down"))
        self.assertEqual(d.eigenvalues.shape, (2, 4, 3))
        self.assertEqual(d.projections.shape, (2, 4, 3, 3, 16))
        self.assertEqual(d.orbital_labels, fx.ORBITALS_16)

    # A file whose <dimension> entries put band outside kpoint gives the same arrays after the transpose.
    def test_permuted_dimensions(self):
        permuted = parse_vasprun(FIXTURES / "permuted_dims_vasprun.xml")
        np.testing.assert_array_equal(permuted.projections, self.data.projections)
        np.testing.assert_array_equal(permuted.eigenvalues, self.data.eigenvalues)

    # Spin-up and spin-down eigenvalues and weights land in their own blocks.
    def test_values(self):
        d = self.data
        np.testing.assert_allclose(d.eigenvalues, expected_eigenvalues(2), atol=TOL)
        np.testing.assert_allclose(d.projections, expected_weights(2, 16, signed=False), atol=TOL)
        self.assertAlmostEqual(d.projections[1, 0, 0, 0, 0] - d.projections[0, 0, 0, 0, 0], 0.3, places=8)
        self.assertEqual(d.orbital_index("fx3"), 15)


class NoncollinearTest(unittest.TestCase):
    """SOC layout: one eigenvalue block, four projection blocks (total, mx, my, mz), signed magnetisation."""

    # Parse the SOC fixture once for all tests.
    @classmethod
    def setUpClass(cls):
        cls.data = parse_vasprun(FIXTURES / "soc_vasprun.xml")

    # One eigenvalue block and four named components.
    def test_layout_and_shapes(self):
        d = self.data
        self.assertEqual(d.layout, LAYOUT_NONCOLLINEAR)
        self.assertEqual(d.projection_components, ("total", "mx", "my", "mz"))
        self.assertEqual(d.eigenvalues.shape, (1, 4, 3))
        self.assertEqual(d.projections.shape, (4, 4, 3, 3, 9))

    # Total charge is nonnegative, magnetisation components keep their sign, and the INCAR override is noted.
    def test_values_and_signs(self):
        d = self.data
        np.testing.assert_allclose(d.projections, expected_weights(4, 9, signed=True), atol=TOL)
        self.assertGreaterEqual(d.projections[0].min(), 0.0)
        self.assertLess(d.projections[1].max(), 0.0)
        self.assertGreater(d.projections[2].min(), 0.0)
        self.assertLess(d.projections[3].max(), 0.0)
        self.assertEqual(d.incar["ISPIN"], 2)
        self.assertEqual(d.parameters["ISPIN"], 1)
        self.assertTrue(any("INCAR ISPIN = 2" in note for note in d.notes))

    # The duplicate eigenvalues inside <projected> are read and agree with the main block.
    def test_projected_eigenvalues(self):
        np.testing.assert_allclose(self.data.projected_eigenvalues, self.data.eigenvalues, atol=0)

    # The raw <kpoints_labels> entries and the line-mode divisions are kept as written.
    def test_kpoint_labels(self):
        d = self.data
        self.assertEqual(d.kpoint_label_entries, [(1, "GAMMA"), (2, "M"), (3, "M"), (4, "K")])
        self.assertEqual((d.kpoint_generation, d.kpoint_divisions), ("listgenerated", 2))


class TrimmedRealFileTest(unittest.TestCase):
    """A fixture cut from Diana's MnTe 9.2-band run: 2 k-points, 4 bands, every header block as VASP wrote it."""

    # Parse the trimmed fixture once for all tests.
    @classmethod
    def setUpClass(cls):
        cls.data = parse_vasprun(FIXTURES / fx.TRIMMED_NAME)

    # Shape (component, kpoint, band, ion, orbital) and the full file's Fermi energy.
    def test_shape_and_fermi_energy(self):
        self.assertEqual(self.data.projections.shape, (4, 2, 4, 4, 9))
        self.assertEqual(self.data.layout, LAYOUT_NONCOLLINEAR)
        self.assertEqual(self.data.species, ["Mn", "Mn", "Te", "Te"])
        self.assertEqual(self.data.efermi, MNTE_EFERMI)
        self.assertEqual(self.data.projection_components, ("total", "mx", "my", "mz"))
        self.assertEqual(self.data.eigenvalues.shape[0], 1)

    # The trimmed arrays equal the matching slice of the full file.
    @unittest.skipUnless(REAL_BAND_RUN.is_file(), "Diana's MnTe band run is not on this machine")
    def test_matches_full_file(self):
        full = parse_vasprun(REAL_BAND_RUN, progress=False)
        self.assertEqual(self.data.efermi, full.efermi)
        np.testing.assert_array_equal(self.data.projections, full.projections[:, :2, :4])
        np.testing.assert_array_equal(self.data.eigenvalues, full.eigenvalues[:, :2, :4])
        np.testing.assert_array_equal(self.data.kpoints, full.kpoints[:2])


class ErrorTest(unittest.TestCase):
    """Missing projections, wrong dimensions, extra blocks, truncated files and compressed input."""

    # No <projected> block: informative error naming LORBIT = 11, unless projections are optional.
    def test_missing_projections(self):
        path = FIXTURES / "no_projected_vasprun.xml"
        with self.assertRaises(ProjectionsMissingError) as caught:
            parse_vasprun(path)
        self.assertIn("LORBIT = 11", str(caught.exception))
        self.assertIn("<dos><partial>", str(caught.exception))
        data = parse_vasprun(path, require_projections=False)
        self.assertIsNone(data.projections)
        self.assertEqual(data.eigenvalues.shape, (1, 4, 3))

    # Without projections: the numbers are not stored, but labels, components and the block's presence are.
    def test_without_projections(self):
        full = parse_vasprun(FIXTURES / "soc_vasprun.xml")
        light = parse_vasprun(FIXTURES / "soc_vasprun.xml", with_projections=False)
        self.assertIsNone(light.projections)
        self.assertIsNone(light.projected_eigenvalues)
        self.assertTrue(light.has_projected_block and full.has_projected_block)
        self.assertEqual((light.orbital_labels, light.projection_components), (full.orbital_labels, full.projection_components))
        np.testing.assert_array_equal(light.eigenvalues, full.eigenvalues)
        self.assertIn("present, not read; shape (4, 4, 3, 3, 9)", "\n".join(light.summary_lines()))
        missing = parse_vasprun(FIXTURES / "no_projected_vasprun.xml", require_projections=False, with_projections=False)
        self.assertFalse(missing.has_projected_block)
        with self.assertRaises(ProjectionsMissingError):
            parse_vasprun(FIXTURES / "no_projected_vasprun.xml", with_projections=False)

    # An empty <field></field> orbital label is refused as a parse error, with and without the projection numbers.
    def test_empty_orbital_field(self):
        head, tail = (FIXTURES / "nonmagnetic_vasprun.xml").read_text(encoding="latin-1").split("<projected>")
        self.assertIn("<field>   py</field>", tail)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vasprun.xml"
            path.write_text(head + "<projected>" + tail.replace("<field>   py</field>", "<field></field>", 1), encoding="latin-1")
            for full in (True, False):
                with self.subTest(with_projections=full):
                    with self.assertRaises(VasprunParseError) as caught:
                        parse_vasprun(path, with_projections=full)
                    self.assertIn("empty <field> orbital label (column 2)", str(caught.exception))

    # A k-point block with one band too few is refused with its position.
    def test_missing_band(self):
        with self.assertRaises(VasprunLayoutError) as caught:
            parse_vasprun(FIXTURES / "missing_band_vasprun.xml")
        self.assertIn("k-point 3", str(caught.exception))

    # More projection blocks than the layout allows are refused instead of guessed.
    def test_extra_block(self):
        with self.assertRaises(VasprunLayoutError) as caught:
            parse_vasprun(FIXTURES / "extra_block_vasprun.xml")
        self.assertIn("more than 2 spin blocks", str(caught.exception))

    # A <projected><array> whose dimension names are not spin, kpoint, band and ion is refused.
    def test_unknown_dimension(self):
        text = (FIXTURES / "nonmagnetic_vasprun.xml").read_text(encoding="latin-1")
        text = text.replace('<dimension dim="1">ion</dimension>\n    <dimension dim="2">band', '<dimension dim="1">atom</dimension>\n    <dimension dim="2">band', 1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vasprun.xml"
            path.write_text(text, encoding="latin-1")
            with self.assertRaises(VasprunLayoutError) as caught:
                parse_vasprun(path)
        self.assertIn("in any order", str(caught.exception))

    # A file cut in the middle gives a parse error that says so.
    def test_truncated(self):
        with self.assertRaises(VasprunParseError) as caught:
            parse_vasprun(FIXTURES / "truncated_vasprun.xml")
        self.assertIn("not complete", str(caught.exception))

    # A k-point list longer than the band data is refused.
    def test_too_few_kpoint_blocks(self):
        text = (FIXTURES / "nonmagnetic_vasprun.xml").read_text(encoding="latin-1")
        extra = '   <v>       0.10000000       0.10000000       0.10000000 </v>\n  </varray>\n  <varray name="weights" >'
        text = text.replace('  </varray>\n  <varray name="weights" >', extra, 1)
        text = text.replace("   <v>       0.25000000 </v>\n  </varray>", "   <v>       0.25000000 </v>\n   <v>       0.25000000 </v>\n  </varray>", 1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vasprun.xml"
            path.write_text(text, encoding="latin-1")
            with self.assertRaises(VasprunLayoutError) as caught:
                parse_vasprun(path)
        self.assertIn("5 k-points are listed", str(caught.exception))

    # A gzip-compressed file gives the same arrays as the plain file.
    def test_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "vasprun.xml.gz"
            with open(FIXTURES / "soc_vasprun.xml", "rb") as src, gzip.open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
            packed = parse_vasprun(target)
        plain = parse_vasprun(FIXTURES / "soc_vasprun.xml")
        np.testing.assert_array_equal(packed.projections, plain.projections)

    # Progress on a .xml.gz file is reported on compressed bytes read and reaches 100 percent.
    def test_gzip_progress(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "vasprun.xml.gz"
            with open(FIXTURES / "soc_vasprun.xml", "rb") as src, gzip.open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                parse_vasprun(target, progress=True)
        self.assertIn("reading vasprun.xml.gz (compressed bytes): 100%", err.getvalue())

    # An empty <v></v> row in the k-point list is a parse error, not a traceback.
    def test_empty_v_row(self):
        text = (FIXTURES / "nonmagnetic_vasprun.xml").read_text(encoding="latin-1")
        text = re.sub(r'(<varray name="kpointlist" >\n)   <v>[^<]*</v>', r"\1   <v></v>", text, count=1)
        self.assertIn("<v></v>", text)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vasprun.xml"
            path.write_text(text, encoding="latin-1")
            with self.assertRaises(VasprunParseError) as caught:
                parse_vasprun(path)
            self.assertIn("kpointlist", str(caught.exception))
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                code = cli.main(["inspect", "--vasprun", str(path)])
        self.assertEqual(code, 2)
        self.assertTrue(err.getvalue().startswith("error: "))

    # A missing file raises FileNotFoundError.
    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            parse_vasprun(FIXTURES / "does_not_exist.xml")


# Text of the nonmagnetic fixture with one earlier <calculation> before the real one; the earlier one has every <r> value replaced.
def two_calculation_text(earlier_projected: bool, last_projected: bool) -> str:
    text = (FIXTURES / "nonmagnetic_vasprun.xml").read_text(encoding="latin-1")
    start, end = text.index(" <calculation>"), text.index(" </calculation>\n") + len(" </calculation>\n")
    block = text[start:end]
    without = re.sub(r"  <projected>.*?</projected>\n", "", block, flags=re.S)
    earlier = re.sub(r"<r>([^<]*)</r>", lambda m: "<r>" + " 0.7777" * len(m.group(1).split()) + " </r>",
                     block if earlier_projected else without)
    return text[:start] + earlier + (block if last_projected else without) + text[end:]


class SeveralCalculationsTest(unittest.TestCase):
    """A file with two <calculation> blocks gives the band data of the last one only."""

    # Write the text to a temporary vasprun.xml and parse it.
    def parse_text(self, text: str):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vasprun.xml"
            path.write_text(text, encoding="latin-1")
            return parse_vasprun(path)

    # The earlier calculation's values are discarded and a note says the last one was used.
    def test_last_calculation_used(self):
        for earlier_projected in (True, False):
            with self.subTest(earlier_projected=earlier_projected):
                data = self.parse_text(two_calculation_text(earlier_projected, last_projected=True))
                plain = parse_vasprun(FIXTURES / "nonmagnetic_vasprun.xml")
                np.testing.assert_array_equal(data.eigenvalues, plain.eigenvalues)
                np.testing.assert_array_equal(data.occupations, plain.occupations)
                np.testing.assert_array_equal(data.projections, plain.projections)
                self.assertIn("2 calculation blocks; last one used", data.notes)

    # Projections only in an earlier calculation are refused rather than paired with the last eigenvalues.
    def test_projections_only_in_earlier_calculation(self):
        with self.assertRaises(VasprunParseError) as caught:
            self.parse_text(two_calculation_text(earlier_projected=True, last_projected=False))
        self.assertIn("earlier one has a <projected> block", str(caught.exception))


class CliTest(unittest.TestCase):
    """The stage 1 "inspect" command prints the summary or a one-line error."""

    # Inspect prints the layout, the shape convention and the orbital labels.
    def test_inspect(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = cli.main(["inspect", "--vasprun", str(FIXTURES / "soc_vasprun.xml")])
        self.assertEqual(code, 0)
        text = out.getvalue()
        self.assertIn("layout: noncollinear", text)
        self.assertIn("(component, kpoint, band, ion, orbital)", text)
        self.assertIn("orbitals: s, py, pz, px, dxy, dyz, dz2, dxz, x2-y2", text)

    # Inspect on a file without projections exits with code 2 and an error on stderr.
    def test_inspect_missing(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code = cli.main(["inspect", "--vasprun", str(FIXTURES / "no_projected_vasprun.xml")])
        self.assertEqual(code, 2)
        self.assertIn("LORBIT = 11", err.getvalue())


# Parse a file with pymatgen's Vasprun, silencing its warnings about missing POTCAR files.
def pymatgen_vasprun(path: Path):
    from pymatgen.io.vasp.outputs import Vasprun
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Vasprun(str(path), parse_projected_eigen=True, parse_potcar_file=False)


# Assert that the parser and pymatgen agree on eigenvalues, Fermi energy, k-points and weights of one file.
def assert_agrees_with_pymatgen(case: unittest.TestCase, path: Path) -> None:
    from pymatgen.electronic_structure.core import Spin
    tool = parse_vasprun(path, progress=False)
    ref = pymatgen_vasprun(path)
    spins = [Spin.up, Spin.down][: tool.eigenvalues.shape[0]]
    for s, spin in enumerate(spins):
        np.testing.assert_allclose(tool.eigenvalues[s], ref.eigenvalues[spin][..., 0], atol=TOL)
        np.testing.assert_allclose(tool.occupations[s], ref.eigenvalues[spin][..., 1], atol=TOL)
    case.assertAlmostEqual(tool.efermi, ref.efermi, places=8)
    np.testing.assert_allclose(tool.kpoints, np.array(ref.actual_kpoints), atol=TOL)
    np.testing.assert_allclose(tool.kpoint_weights, np.array(ref.actual_kpoints_weights), atol=TOL)
    case.assertEqual(tool.species, list(ref.atomic_symbols))
    if tool.layout == LAYOUT_NONCOLLINEAR:
        np.testing.assert_allclose(tool.projections[0], ref.projected_eigenvalues[Spin.up], atol=TOL)
        np.testing.assert_allclose(np.moveaxis(tool.projections[1:4], 0, -1), ref.projected_magnetization, atol=TOL)
    else:
        for c, spin in enumerate(spins):
            np.testing.assert_allclose(tool.projections[c], ref.projected_eigenvalues[spin], atol=TOL)


class PymatgenCrossCheckTest(unittest.TestCase):
    """Independent reference: pymatgen's Vasprun on the three synthetic layouts and on one real band run."""

    # The three synthetic layouts agree with pymatgen to 1e-6.
    def test_fixtures(self):
        for name in ("nonmagnetic_vasprun.xml", "collinear_f_vasprun.xml", "soc_vasprun.xml"):
            with self.subTest(name=name):
                assert_agrees_with_pymatgen(self, FIXTURES / name)

    # Diana's MnTe SOC band run (VASP 6.3.2, 100 k-points, 60 bands, 4 ions) agrees with pymatgen to 1e-6.
    @unittest.skipUnless(REAL_BAND_RUN.is_file(), "Diana's MnTe band run is not on this machine")
    def test_real_mnte_band_run(self):
        assert_agrees_with_pymatgen(self, REAL_BAND_RUN)
        data = parse_vasprun(REAL_BAND_RUN, progress=False)
        self.assertEqual(data.projections.shape, (4, 100, 60, 4, 9))
        self.assertEqual(data.species, ["Mn", "Mn", "Te", "Te"])
        self.assertAlmostEqual(data.efermi, MNTE_EFERMI, places=8)
        self.assertEqual(data.projection_components, ("total", "mx", "my", "mz"))
        self.assertEqual(data.eigenvalues.shape[0], 1)


if __name__ == "__main__":
    unittest.main()
