"""Unittest wrapper around tests/accuracy.py: runs the 20-tree measurement and checks it completes with sane numbers."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import accuracy  # noqa: E402


class Accuracy(unittest.TestCase):
    """Run the measurement once and check every run was found and scored."""

    # Generate and import the 20 trees, print the report.
    @classmethod
    def setUpClass(cls):
        cls.scores, cls.mismatches = accuracy.measure(20, 1)
        accuracy.report(cls.scores, cls.mismatches)

    # Every truth run is found by the importer, so no mismatch reads "not found".
    def test_all_runs_found(self):
        self.assertFalse([m for m in self.mismatches if str(m[4]).startswith("not found")])

    # Phase and status stay at or above the measured level on this branch; a drop means a regression.
    def test_floor(self):
        for field, floor in (("phase", 0.95), ("status", 0.95), ("parent", 0.95)):
            ok, tot = self.scores[field]
            self.assertGreaterEqual(ok / tot, floor, field)


if __name__ == "__main__":
    unittest.main()
