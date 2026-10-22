"""Diana's 22 Oct answers: a standalone Wannier90 run counts as completed, and SOC suggestions carry ISPIN = 2."""

import tempfile
import unittest
from pathlib import Path

from altermagagent.branches import MISSING, next_calculation
from altermagagent.importer import find_completion


class DianaAnswersTest(unittest.TestCase):
    # A folder with no OUTCAR but a wannier90.wout ending in "All done" is a completed run.
    def test_wannier90_wout_is_completion_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "wannier90.wout").write_text("some output\n All done: wannier90 exiting\n")
            self.assertEqual(find_completion(folder, {"wannier90.wout", "wannier90.win"}), "wannier90.wout")
            (folder / "wannier90.wout").write_text("Wannierisation convergence criteria satisfied\n")
            self.assertIsNone(find_completion(folder, {"wannier90.wout"}))

    # The SCF suggestion for a SOC branch names ISPIN = 2 next to LNONCOLLINEAR and LSORBIT.
    def test_soc_scf_suggestion_has_ispin_2(self):
        branch = {"name": "collinear-SOC", "realization": "collinear"}
        steps = [{"step": "c", "state": MISSING, "action": ""}]
        text = next_calculation(branch, steps, [], True)
        self.assertIn("ISPIN = 2, LNONCOLLINEAR = .TRUE., LSORBIT = .TRUE.", text)
        text = next_calculation(branch, steps, [], False)
        self.assertIn("ISPIN = 2, LSORBIT = .FALSE.", text)


if __name__ == "__main__":
    unittest.main()
