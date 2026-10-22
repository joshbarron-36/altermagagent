"""Wannier comparison: read_hr and error_metrics on a tiny synthetic Hamiltonian, and the demo tree Hamiltonian."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from altermagagent.wannier_compare import (WS_FALLBACK_NOTE, auto_band_start, band_match_text, error_box_text,
                                           error_metrics, interpolate, interpolation_terms, read_hr)

DEMO_HR = Path(__file__).resolve().parents[1] / "demo" / "CrSe_demo"


# Text of a two-orbital chain Hamiltonian: on-site e0 and e1, hopping t between neighbour cells along a.
def chain_hr(e0=0.5, e1=-0.5, t=0.25):
    vectors = [(-1, 0, 0), (0, 0, 0), (1, 0, 0)]
    lines = [" synthetic chain", "2", "3", "1 1 1"]
    for r in vectors:
        for n in (1, 2):
            for m in (1, 2):
                val = (e0 if m == 1 else e1) if r == (0, 0, 0) and m == n else (t if m == n and r != (0, 0, 0) else 0.0)
                lines.append(f"{r[0]:5d}{r[1]:5d}{r[2]:5d}{m:5d}{n:5d}{val:12.6f}{0.0:12.6f}")
    return "\n".join(lines) + "\n"


class WannierCompare(unittest.TestCase):
    """Check parsing, interpolation and the error numbers on known values."""

    # read_hr returns the dimensions, R vectors and matrices, and the interpolated bands are e + 2 t cos(2 pi k).
    def test_read_hr_and_interpolate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wannier90_hr.dat"
            path.write_text(chain_hr())
            hr = read_hr(path)
        self.assertEqual(hr["num_wann"], 2)
        self.assertEqual(hr["matrices"].shape, (3, 2, 2))
        self.assertEqual(sorted(map(tuple, hr["vectors"])), [(-1, 0, 0), (0, 0, 0), (1, 0, 0)])
        k = np.array([[0.0, 0, 0], [0.25, 0, 0], [0.5, 0, 0]])
        bands, residual = interpolate(hr, k)
        cos = np.cos(2 * np.pi * k[:, 0])
        np.testing.assert_allclose(bands, [-0.5 + 0.5 * cos, 0.5 + 0.5 * cos], atol=1e-12)
        self.assertLess(residual, 1e-12)

    # A truncated file is rejected rather than read as a smaller Hamiltonian.
    def test_read_hr_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wannier90_hr.dat"
            path.write_text("\n".join(chain_hr().splitlines()[:-2]) + "\n")
            with self.assertRaises(ValueError):
                read_hr(path)

    # MAE, RMSE and max in meV over the samples inside the window, with the band offset applied.
    def test_error_metrics(self):
        reference = np.array([[-5.0, -5.0, -5.0], [-0.5, 0.0, 0.5], [0.2, 0.4, 2.0]])
        wannier = np.array([[-0.499, 0.002, 0.497], [0.204, 0.4, 2.5]])
        m = error_metrics(reference, wannier, 2, -1.0, 1.0)
        errors = np.array([1, 2, 3, 4, 0]) * 1.0
        self.assertEqual(m["samples"], 5)
        self.assertAlmostEqual(m["mae_meV"], errors.mean(), places=6)
        self.assertAlmostEqual(m["rmse_meV"], np.sqrt((errors ** 2).mean()), places=6)
        self.assertAlmostEqual(m["max_meV"], 4.0, places=6)
        self.assertIsNone(error_metrics(reference, wannier, None, -1.0, 1.0))
        with self.assertRaises(ValueError):
            error_metrics(reference, wannier, 3, -1.0, 1.0)

    # The automatic offset is the smallest-error one and says it is not verified.
    def test_auto_offset_note(self):
        reference = np.array([[-5.0, -5.0, -5.0], [-0.5, 0.0, 0.5], [0.2, 0.4, 0.6]])
        wannier = np.array([[-0.499, 0.002, 0.497], [0.204, 0.4, 0.6]])
        with tempfile.TemporaryDirectory() as tmp:
            start, note = auto_band_start(reference, wannier, Path(tmp) / "wannier90_hr.dat", -1.0, 1.0)
        self.assertEqual(start, 2)
        self.assertEqual(note, "band offset 2 chosen by smallest error, not verified")

    # The error box and the step (e) line carry the 95th percentile, the count above 20 meV and the offset note.
    def test_box_and_step_e_line(self):
        reference = np.array([[-0.5, 0.0, 0.5]])
        m = error_metrics(reference, np.array([[-0.47, 0.0, 0.5]]), 1, -1.0, 1.0)
        self.assertEqual(m["samples_over_20_meV"], 1)
        m["mapping_source"] = "band offset 1 chosen by smallest error, not verified"
        box = error_box_text({"metrics": m, "metric_emin": -1, "metric_emax": 1})
        for text in ("95th percentile", "1 of 3 samples above 20 meV", "chosen by smallest error, not verified"):
            self.assertIn(text, box)
        line = band_match_text({"status": "ok", "metrics": m, "window": [-1, 1], "reference_run": "bands",
                                "hr_kind": "raw", "png": "x.png",
                                "evidence_state": "preliminary: structure matched; required inputs differ or are missing"})
        self.assertTrue(line.startswith("preliminary: band match:"))
        self.assertIn("not verified", line)

    # Without wsvec and without a use_ws_distance=false declaration the direct-HR fallback is named; a declaration silences it.
    def test_ws_fallback_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wannier90_hr.dat"
            path.write_text(chain_hr())
            with self.assertRaises(ValueError):
                interpolation_terms(path, read_hr(path))
            _, convention = interpolation_terms(path, read_hr(path), use_ws_distance=False)
            self.assertEqual(convention["convention_note"], WS_FALLBACK_NOTE)
            (Path(tmp) / "wannier90.win").write_text("num_wann = 2\nuse_ws_distance = false\n")
            _, convention = interpolation_terms(path, read_hr(path))
            self.assertIsNone(convention["convention_note"])

    # The demo tree's wannier90_hr.dat is a header excerpt only, so the real comparison is skipped with the reason.
    def test_demo_tree(self):
        files = sorted(DEMO_HR.rglob("wannier90_hr.dat")) if DEMO_HR.exists() else []
        if not files:
            self.skipTest("demo tree not generated (python demo/make_demo.py) or it has no wannier90_hr.dat")
        try:
            read_hr(files[0])
        except ValueError as exc:
            self.skipTest(f"demo wannier90_hr.dat is a header excerpt, not a full matrix ({exc})")
        self.fail("demo Hamiltonian unexpectedly complete; add a real comparison test")


if __name__ == "__main__":
    unittest.main()
