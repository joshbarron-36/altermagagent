"""Stage 2 of the orbital plotter: atom and orbital selection, groups, name matching, weight sums and normalisation."""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "fixtures"))

import orbital_vasprun_builder as fx  # noqa: E402
from vasp_orbital_plotter.parser import parse_vasprun  # noqa: E402
from vasp_orbital_plotter.projections import (SelectionError, available_groups, build_selection,  # noqa: E402
                                              orbital_column, parse_atoms, resolve_orbital, resolve_orbitals,
                                              select_component, select_weights, shell_columns)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "orbital"
TOL = 1e-9
VASP_LABELS = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"]


# Weights of a selection from the block that select_component picks.
def weights_of(data, selection, normalize: bool = False) -> np.ndarray:
    return select_weights(data, selection, select_component(data)[0], normalize=normalize)


# Hand sum of the fixture formula over 0-based ions and orbital columns for component c, shape (kpoint, band).
def hand_sum(ions0: list[int], columns: list[int], c: int = 0, signed: bool = False) -> np.ndarray:
    out = np.zeros((len(fx.KPOINTS), fx.NBANDS))
    for k in range(len(fx.KPOINTS)):
        for b in range(fx.NBANDS):
            out[k, b] = sum(fx.weight(c, k, b, i, o, signed=signed) for i in ions0 for o in columns)
    return out


class AtomSelection(unittest.TestCase):
    """Ions are named in VASP's 1-based numbering and converted to 0-based indices explicitly."""

    # Parse the nonmagnetic fixture once: ions 1:Cr, 2:Se, 3:Cr.
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = parse_vasprun(FIXTURES / "nonmagnetic_vasprun.xml", progress=False)

    # Numbers, ranges and element symbols give 1-based ion numbers without duplicates.
    def test_numbers_ranges_and_elements(self) -> None:
        self.assertEqual(parse_atoms(["1", "3"], self.data), [1, 3])
        self.assertEqual(parse_atoms(["1-2"], self.data), [1, 2])
        self.assertEqual(parse_atoms(["Cr"], self.data), [1, 3])
        self.assertEqual(parse_atoms(["Cr", "1", "2-3"], self.data), [1, 3, 2])

    # Ion 0, an ion beyond the last, a backwards range and an unknown element are refused with a SelectionError.
    def test_bad_atoms(self) -> None:
        for items in (["0"], ["4"], ["3-1"], ["Fe"]):
            with self.subTest(items=items), self.assertRaises(SelectionError):
                parse_atoms(items, self.data)

    # VASP ion 1 reads array index 0 and ion 3 reads index 2: the weight of ion n equals the formula at ion n - 1.
    def test_one_based_to_zero_based(self) -> None:
        for ion in (1, 2, 3):
            selection = build_selection(self.data, atoms=str(ion), orbitals="s")
            np.testing.assert_allclose(weights_of(self.data, selection), hand_sum([ion - 1], [0]), atol=TOL)

    # The selection label names the element when the ions are exactly that element's ions.
    def test_label(self) -> None:
        self.assertEqual(build_selection(self.data, elements="Cr", orbitals="d").label, "Cr d")
        self.assertEqual(build_selection(self.data, atoms="1", orbitals="pz").label, "Cr1 pz")
        self.assertEqual(build_selection(self.data).label, "all atoms all orbitals")


class OrbitalNames(unittest.TestCase):
    """The file's <field> labels are the authority; one normalising rule maps user names and groups onto them."""

    # Case, spaces, "_", "^", "-" and the shell letter are optional, on the user's side and on the file's side.
    def test_name_matching(self) -> None:
        for name in ("x2-y2", "dx2-y2", "DX2-Y2", "dx2y2", "d_x^2-y^2"):
            self.assertEqual(orbital_column(name, VASP_LABELS), 8, name)
        for name in ("dz2", "z2", "Dz2", "d_z^2"):
            self.assertEqual(orbital_column(name, VASP_LABELS), 6, name)
        self.assertEqual(orbital_column("xy", VASP_LABELS), 4)
        self.assertIsNone(orbital_column("dx", VASP_LABELS))
        self.assertEqual(resolve_orbital("dx2-y2", VASP_LABELS), [8])

    # A file that spells the last d label with its shell letter is matched by VASP's spelling too.
    def test_file_with_other_spelling(self) -> None:
        labels = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2-y2"]
        self.assertEqual(resolve_orbital("x2-y2", labels), [8])
        self.assertEqual(resolve_orbital("eg", labels), [6, 8])

    # Shells come from the file labels: x2-y2 has no shell letter and continues the d shell before it.
    def test_shell_columns(self) -> None:
        self.assertEqual(shell_columns(VASP_LABELS), {"s": [0], "p": [1, 2, 3], "d": [4, 5, 6, 7, 8]})
        self.assertEqual(shell_columns(fx.ORBITALS_16)["f"], list(range(9, 16)))
        self.assertEqual(shell_columns(["", "s", " ", "py"]), {"s": [1, 2], "p": [3]})

    # Shell groups resolve to their columns in file order, t2g and eg in the order of their definition.
    def test_groups(self) -> None:
        self.assertEqual(resolve_orbital("p", VASP_LABELS), [1, 2, 3])
        self.assertEqual(resolve_orbital("d", VASP_LABELS), [4, 5, 6, 7, 8])
        self.assertEqual(resolve_orbital("t2g", VASP_LABELS), [4, 7, 5])
        self.assertEqual(resolve_orbital("eg", VASP_LABELS), [6, 8])
        self.assertEqual(resolve_orbital("all", VASP_LABELS), list(range(9)))
        self.assertEqual(resolve_orbitals(["dxy", "t2g", "s"], VASP_LABELS), [4, 7, 5, 0])
        self.assertEqual(available_groups(VASP_LABELS), ["s", "p", "d", "t2g", "eg"])

    # A group with a missing member is refused, while groups whose members all exist still work.
    def test_group_needs_all_members(self) -> None:
        labels = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz"]
        for group, missing in (("eg", "missing: x2-y2"), ("d", "missing: 1 of the 5 d components")):
            with self.subTest(group=group), self.assertRaises(SelectionError) as caught:
                resolve_orbital(group, labels)
            self.assertIn(missing, str(caught.exception))
        self.assertEqual(resolve_orbital("t2g", labels), [4, 7, 5])
        self.assertNotIn("eg", available_groups(labels))
        self.assertNotIn("d", available_groups(labels))

    # An unknown orbital is refused with a message that lists the file's orbitals and groups.
    def test_unknown_orbital_lists_available(self) -> None:
        with self.assertRaises(SelectionError) as caught:
            resolve_orbital("dq", VASP_LABELS)
        message = str(caught.exception)
        self.assertIn("available orbitals: s, py, pz, px, dxy, dyz, dz2, dxz, x2-y2", message)
        self.assertIn("groups: s, p, d, t2g, eg, all", message)

    # The 16-column file of the collinear fixture offers the f group.
    def test_f_group(self) -> None:
        self.assertEqual(resolve_orbital("f", fx.ORBITALS_16), list(range(9, 16)))


class WeightSums(unittest.TestCase):
    """The weight of a state is the plain sum of stored values over the selected ions and orbitals."""

    # Parse the three layouts once.
    @classmethod
    def setUpClass(cls) -> None:
        cls.nonmag = parse_vasprun(FIXTURES / "nonmagnetic_vasprun.xml", progress=False)
        cls.collinear = parse_vasprun(FIXTURES / "collinear_f_vasprun.xml", progress=False)
        cls.soc = parse_vasprun(FIXTURES / "soc_vasprun.xml", progress=False)

    # Ions 1 and 3 with dxy and pz: four (ion, orbital) pairs summed per state.
    def test_multiple_atoms_and_orbitals(self) -> None:
        selection = build_selection(self.nonmag, atoms="1,3", orbitals="dxy,pz")
        np.testing.assert_allclose(weights_of(self.nonmag, selection), hand_sum([0, 2], [4, 2]), atol=TOL)

    # An element and a group: Cr t2g equals the hand sum over ions 0 and 2 and the three t2g columns.
    def test_element_and_group(self) -> None:
        selection = build_selection(self.nonmag, elements="Cr", orbitals="t2g")
        np.testing.assert_allclose(weights_of(self.nonmag, selection), hand_sum([0, 2], [4, 7, 5]), atol=TOL)

    # Known value: Cr d at k-point 1, band 3 (VASP numbering) is a literal number from the fixture formula.
    def test_literal_value(self) -> None:
        selection = build_selection(self.nonmag, elements="Cr", orbitals="d")
        # 10 pairs: 0.2 * 10 + 0.01 * (1 + 3) * 5 + 0.001 * (5 + 6 + 7 + 8 + 9) * 2 = 2.0 + 0.2 + 0.07
        self.assertAlmostEqual(float(weights_of(self.nonmag, selection)[0, 2]), 2.27, places=9)

    # Repeated --select terms are summed over the union of pairs, so a pair named twice counts once.
    def test_select_terms_union(self) -> None:
        selection = build_selection(self.nonmag, select=["1:s", "Cr:s,pz", "2:eg"])
        self.assertEqual(sorted(selection.pairs), [(1, 0), (1, 2), (2, 6), (2, 8), (3, 0), (3, 2)])
        expected = hand_sum([0, 2], [0, 2]) + hand_sum([1], [6, 8])
        np.testing.assert_allclose(weights_of(self.nonmag, selection), expected, atol=TOL)

    # --select cannot be mixed with the plain options, and a term without a colon is refused.
    def test_select_errors(self) -> None:
        with self.assertRaises(SelectionError):
            build_selection(self.nonmag, atoms="1", select=["Cr:d"])
        with self.assertRaises(SelectionError):
            build_selection(self.nonmag, select=["Cr-d"])

    # --normalize divides by the sum over all ions and orbitals of the same state.
    def test_normalisation(self) -> None:
        selection = build_selection(self.nonmag, atoms="2", orbitals="p")
        total = hand_sum([0, 1, 2], list(range(9)))
        np.testing.assert_allclose(weights_of(self.nonmag, selection, normalize=True),
                                   hand_sum([1], [3, 1, 2]) / total, atol=TOL)
        whole = build_selection(self.nonmag)
        np.testing.assert_allclose(weights_of(self.nonmag, whole, normalize=True), 1.0, atol=TOL)

    # Raw weights are the default: they are not scaled to one.
    def test_raw_by_default(self) -> None:
        whole = build_selection(self.nonmag)
        np.testing.assert_allclose(weights_of(self.nonmag, whole), hand_sum([0, 1, 2], list(range(9))), atol=TOL)

    # Noncollinear files use the "total" block, never mx, my or mz.
    def test_noncollinear_uses_total(self) -> None:
        index, label, note = select_component(self.soc)
        self.assertEqual((index, label, note), (0, "total", None))
        selection = build_selection(self.soc, atoms="2", orbitals="d")
        np.testing.assert_allclose(weights_of(self.soc, selection), hand_sum([1], [4, 5, 6, 7, 8], signed=True), atol=TOL)

    # Collinear files use the spin-up block in stage 2 and say so in a note.
    def test_collinear_uses_spin_up_with_note(self) -> None:
        index, label, note = select_component(self.collinear)
        self.assertEqual((index, label), (0, "up"))
        self.assertIn("spin-up", note)
        selection = build_selection(self.collinear, atoms="3", orbitals="f")
        np.testing.assert_allclose(weights_of(self.collinear, selection), hand_sum([2], list(range(9, 16))), atol=TOL)


if __name__ == "__main__":
    unittest.main()
