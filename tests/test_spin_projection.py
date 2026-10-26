"""Spin projection: Pauli expectations on a synthetic spinor model, spin-order detection, degenerate rotation, errors."""

import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from _fixtures import collinear_wannier_project, demo_tree, temp_index
from altermagagent import importer
from altermagagent.spin_projection import (NO_WIN_NOTE, NOT_CHECKED, apply_operator, degenerate_groups, expectations, group_rotation, project_path,
                                           refresh_spin_step, resolve_order, spin_line, spin_operators,
                                           spin_order_from_hamiltonian, spin_projection, summary_lines)
from altermagagent.wannier_compare import hr_parsed, read_hr

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
# Rows (up, down) of orbitals 0 and 1 in the Wannier basis for each spin order of the two-orbital model.
ROWS = {"udud": ((0, 1), (2, 3)), "uudd": ((0, 2), (1, 3))}


# Text of a two-orbital spinor chain: orbital 0 at -1 eV with 0.2 (cos, sin, 0).sigma, orbital 1 at +1 eV with 0.2 sigma_z.
def spinor_hr(order, header=" synthetic spinor chain"):
    blocks = {(0, 0, 0): [-np.eye(2), np.eye(2) + 0.2 * SZ],
              (1, 0, 0): [0.1 * SX - 0.1j * SY, np.zeros((2, 2))],
              (-1, 0, 0): [0.1 * SX + 0.1j * SY, np.zeros((2, 2))]}
    lines = [header, "4", "3", "1 1 1"]
    for r, orbitals in blocks.items():
        full = np.zeros((4, 4), complex)
        for i, block in enumerate(orbitals):
            idx = ROWS[order][i]
            full[np.ix_(idx, idx)] = block
        for n in range(4):
            for m in range(4):
                lines.append(f"{r[0]:5d}{r[1]:5d}{r[2]:5d}{m + 1:5d}{n + 1:5d}{full[m, n].real:12.6f}{full[m, n].imag:12.6f}")
    return "\n".join(lines) + "\n"


class SpinProjection(unittest.TestCase):
    """Check the expectations on known spin textures and the guards around them."""

    # Parse the two-orbital spinor Hamiltonian once per spin order.
    @classmethod
    def setUpClass(cls):
        cls.hr = {}
        with tempfile.TemporaryDirectory() as tmp:
            for order in ROWS:
                path = Path(tmp) / f"wannier90_hr_{order}.dat"
                path.write_text(spinor_hr(order))
                cls.hr[order] = read_hr(path)

    # Both orders give the known texture: orbital 0 spins follow (cos, sin, 0) of 2 pi k, orbital 1 is pure sigma_z.
    def test_known_texture_both_orders(self):
        k = np.array([[0.0, 0, 0], [0.125, 0, 0], [0.3, 0, 0]])
        theta = 2 * np.pi * k[:, 0]
        expected = {"x": [-np.cos(theta), np.cos(theta), 0 * theta, 0 * theta],
                    "y": [-np.sin(theta), np.sin(theta), 0 * theta, 0 * theta],
                    "z": [0 * theta, 0 * theta, -1 + 0 * theta, 1 + 0 * theta]}
        for order, hr in self.hr.items():
            energies, spins, residual = project_path(hr, k, order)
            np.testing.assert_allclose(energies, [[-1.2] * 3, [-0.8] * 3, [0.8] * 3, [1.2] * 3], atol=1e-9)
            for c in "xyz":
                np.testing.assert_allclose(spins[c], expected[c], atol=1e-9, err_msg=f"{order} sigma_{c}")
            self.assertLess(residual, 1e-12)
            other = "uudd" if order == "udud" else "udud"
            _, wrong, _ = project_path(hr, k, other)
            self.assertFalse(np.allclose(wrong["z"], expected["z"], atol=1e-6))

    # The operators are Hermitian Pauli matrices on the spin rows of each orbital.
    def test_operators(self):
        for order in ("udud", "uudd"):
            ops = spin_operators(4, order)
            for c, pauli in zip("xyz", (SX, SY, SZ)):
                np.testing.assert_allclose(ops[c], ops[c].conj().T)
                for i in range(2):
                    idx = ROWS[order][i]
                    np.testing.assert_allclose(ops[c][np.ix_(idx, idx)], pauli)
        with self.assertRaises(ValueError):
            spin_operators(3, "udud")

    # The order comes from "->udud" in the header, then _udud or _uudd in the file name, else None.
    def test_order_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            cases = [("wannier90_hr_symmetrized.dat", "Symmed by someone on 2026-09-15; permuted uudd->udud", "udud"),
                     ("wannier90_hr_uudd.dat", " written on  3Sep2026 at 17:13:34", "uudd"),
                     ("wannier90_hr_symmetrized_udud.dat", " written on  3Sep2026", "udud"),
                     ("wannier90_hr.dat", " written on  3Sep2026 at 17:13:34", None)]
            for name, header, want in cases:
                path = Path(tmp) / name
                path.write_text(header + "\n4\n")
                self.assertEqual(spin_order_from_hamiltonian(path), want, name)

    # A spin-degenerate pair in a mixed basis gives 0 before the rotation and +1 and -1 after it.
    def test_degenerate_rotation(self):
        vectors = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
        self.assertAlmostEqual(float(np.real(vectors[:, 0].conj() @ SZ @ vectors[:, 0])), 0.0)
        diagonal, rotation = group_rotation(vectors, ([0], [1]), "z")
        rotated = vectors @ rotation
        np.testing.assert_allclose(sorted(np.real(np.sum(rotated.conj() * (SZ @ rotated), axis=0))), [-1, 1], atol=1e-12)
        np.testing.assert_allclose(diagonal, [-1, 1], atol=1e-12)
        _, values = expectations(np.zeros((2, 2), complex), ([0], [1]))
        np.testing.assert_allclose(sorted(values["z"]), [-1, 1], atol=1e-12)

    # Each component is diagonalised with its own operator: a degenerate pair gives +1 and -1 for x, y and z alike.
    def test_degenerate_rotation_per_component(self):
        up_down = np.eye(2, dtype=complex)
        self.assertAlmostEqual(float(np.real(up_down[:, 0].conj() @ SX @ up_down[:, 0])), 0.0)
        np.testing.assert_allclose(group_rotation(up_down, ([0], [1]), "x")[0], [-1, 1], atol=1e-12)
        _, values = expectations(np.zeros((2, 2), complex), ([0], [1]))
        for c in "xyz":
            np.testing.assert_allclose(sorted(values[c]), [-1, 1], atol=1e-12, err_msg=f"sigma_{c}")

    # apply_operator equals the dense Pauli matrix times the vectors for every component and both orders.
    def test_apply_operator_matches_dense(self):
        vectors = np.random.default_rng(1).normal(size=(4, 3)) + 1j * np.random.default_rng(2).normal(size=(4, 3))
        for order in ROWS:
            ops = spin_operators(4, order)
            rows = (np.array([ROWS[order][0][0], ROWS[order][1][0]]), np.array([ROWS[order][0][1], ROWS[order][1][1]]))
            for c in "xyz":
                np.testing.assert_allclose(apply_operator(vectors, rows, c), ops[c] @ vectors, atol=1e-14)

    # Groups are bounded by their first state: 0, 0.8e-5 and 1.6e-5 eV give (0, 1) and (2), not one chained group.
    def test_bounded_spread_groups(self):
        self.assertEqual(degenerate_groups(np.array([0, 0.8e-5, 1.6e-5]), 1e-5), [(0, 2), (2, 3)])

    # The step (e) suggestion appears only when the band comparison parsed the Hamiltonian.
    def test_suggestion_gate(self):
        self.assertIn("Next: spin projection of w/wannier90_hr.dat", spin_line("w", {"status": "ok"}))
        self.assertIn("confirm the Wannier spin order first", spin_line("w", {"status": "no_mapping"}))
        for res in (None, {"status": "no_hr"}, {"status": "error"}, {"status": "no_reference"}):
            self.assertIsNone(spin_line("w", res))

    # hr_parsed is true exactly for the band-match statuses that read wannier90_hr.dat.
    def test_hr_parsed(self):
        self.assertTrue(hr_parsed({"status": "ok"}))
        self.assertTrue(hr_parsed({"status": "no_mapping"}))
        for res in (None, {}, {"status": "no_hr"}, {"status": "error"}, {"status": "no_reference"}):
            self.assertFalse(hr_parsed(res))

    # An explicit Fermi energy is printed as given, not turned into a path below the working folder.
    def test_explicit_fermi_source_line(self):
        result = {"caption": "c", "spin_order": "udud", "spin_order_source": "given", "emin": -1, "emax": 1,
                  "states_in_window": 0, "window_stats": {c: None for c in "xyz"}, "fermi_eV": 6.0,
                  "fermi_source": "explicit --fermi-eV", "n_kpoints": 3, "reference_run": "bands", "num_wann": 4,
                  "max_hermiticity_residual_eV": 0.0, "notes": []}
        line = summary_lines(result)[2]
        self.assertIn("from explicit --fermi-eV;", line)
        self.assertNotIn("~", line)
        self.assertNotIn("/", line)

    # A step stored before the spin entry existed keeps the literal line until refreshed, then gets the suggestion.
    def test_refresh_spin_step_without_stored_entry(self):
        step = {"evidence": ["band match: ok", NOT_CHECKED, "wannier90.win: num_wann 4"]}
        refresh_spin_step(step, "w", {"status": "ok"})
        self.assertTrue(step["spin"]["available"])
        self.assertTrue(step["evidence"][1].startswith("Next: spin projection of w/wannier90_hr.dat"))
        self.assertEqual(step["evidence"][2], "wannier90.win: num_wann 4")
        refresh_spin_step(step, "w", {"status": "error"})
        self.assertEqual(step["evidence"][1], NOT_CHECKED)
        self.assertFalse(step["spin"]["available"])

    # A given order beats the file, a stated order is detected, and udud is assumed when nothing states it.
    def test_resolve_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            stated, plain = Path(tmp) / "wannier90_hr_uudd.dat", Path(tmp) / "wannier90_hr.dat"
            for path in (stated, plain):
                path.write_text(" written on  3Sep2026\n4\n")
            self.assertEqual(resolve_order(stated, "udud"), ("udud", "given; the file states uudd"))
            self.assertEqual(resolve_order(stated, "uudd"), ("uudd", "given"))
            self.assertEqual(resolve_order(stated), ("uudd", "detected"))
            self.assertEqual(resolve_order(plain), ("udud", "assumed"))
            self.assertEqual(resolve_order(plain, "uudd"), ("uudd", "given"))


class SpinProjectionProject(unittest.TestCase):
    """Run spin_projection against a demo tree in a temporary index."""

    # Write the demo tree once into a temporary folder.
    @classmethod
    def setUpClass(cls):
        cls.tmp, cls.demo = demo_tree()

    # Remove the temporary folder.
    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    # Point the index at a fresh temporary home for each test.
    def setUp(self):
        self.out = temp_index(self, self.tmp.name) / "figures"

    # A project without any wannier90_hr.dat gives a clear ValueError and writes no figure.
    def test_no_hamiltonian(self):
        copy = Path(self.tmp.name) / "no_hr"
        shutil.copytree(self.demo, copy)
        for path in copy.rglob("wannier90_hr.dat"):
            path.unlink()
        importer.import_project(copy, "no_hr")
        with self.assertRaisesRegex(ValueError, "no unique indexed wannier90_hr.dat"):
            spin_projection("no_hr", out_dir=self.out)
        with self.assertRaisesRegex(ValueError, "not indexed"):
            spin_projection("missing", out_dir=self.out)
        self.assertFalse(self.out.exists())

    # Two calls write into two fresh folders, so concurrent clicks never share a PNG or JSON path.
    def test_unique_output_folders(self):
        collinear_wannier_project(self.demo, self.tmp.name, "spin_twice", spinor_hr("udud"), "num_wann = 4\nspinors = true\n")
        first = spin_projection("spin_twice", run_name="wannier", component="all", out_dir=self.out)
        second = spin_projection("spin_twice", run_name="wannier", component="all", out_dir=self.out)
        self.assertNotEqual(first["paths"]["png"], second["paths"]["png"])
        self.assertNotEqual(Path(first["paths"]["png"]).parent, Path(second["paths"]["png"]).parent)
        for res in (first, second):
            self.assertTrue(all(Path(path).is_file() for path in res["paths"].values()))
            self.assertIn("each component diagonalised inside degenerate groups", res["caption"])
            self.assertEqual(res["degenerate_rotation"], "each component diagonalised inside degenerate groups")

    # spinors = false in wannier90.win is refused; without any .win the spinor basis is assumed and noted.
    def test_spinor_check(self):
        collinear_wannier_project(self.demo, self.tmp.name, "collinear_win", spinor_hr("udud"), "num_wann = 4\nspinors = false\n")
        with self.assertRaisesRegex(ValueError, "spinors = false"):
            spin_projection("collinear_win", run_name="wannier", out_dir=self.out)
        collinear_wannier_project(self.demo, self.tmp.name, "no_win", spinor_hr("udud"))
        res = spin_projection("no_win", run_name="wannier", out_dir=self.out)
        self.assertIn(NO_WIN_NOTE, res["notes"])
        self.assertIn(NO_WIN_NOTE, res["caption"])


if __name__ == "__main__":
    unittest.main()
