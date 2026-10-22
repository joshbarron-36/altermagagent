"""Band tick label normalisation: Greek names, subscripts, path breaks and stray weight digits."""

import unittest

from altermagagent.bands import build_ticks, pretty_label


class BandLabels(unittest.TestCase):
    """Check the printed form of KPOINTS and KLABELS names."""

    # Greek names and subscripts map to their symbols.
    def test_names(self):
        cases = {"GAMMA": "Γ", "G": "Γ", "Gamma": "Γ", "\\Gamma": "Γ", "SIGMA_0": "Σ₀", "C_0": "C₀", "Y_1": "Y₁",
                 "Y1": "Y₁", "DELTA": "Δ", "LAMBDA_2": "Λ₂", "1    GAMMA": "Γ", "1SIGMA_": "Σ", "1Y|1GAMMA": "Y|Γ",
                 "1C_0|1SIGMA_0": "C₀|Σ₀", "Γ|Γ": "Γ", "1": ""}
        for raw, want in cases.items():
            self.assertEqual(pretty_label(raw), want, raw)
        self.assertEqual(pretty_label("1", keep_numbers=True), "1")

    # Labels at one position join with "|", close ticks merge, repeats collapse and empty labels drop.
    def test_ticks(self):
        ticks = build_ticks([(0.0, "Γ"), (1.0, "Y"), (1.0, "Γ"), (1.0, "Γ"), (5.0, "Z"), (5.05, "T"), (10.0, "")], 10.0)
        self.assertEqual(ticks, [(0.0, "Γ"), (1.0, "Y|Γ"), (5.0, "Z|T")])


if __name__ == "__main__":
    unittest.main()
