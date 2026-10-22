"""Physicist validations upgrade a matching tag and flag a differing value."""

import unittest

from altermagagent.summary import PENDING_LABEL, find_validation, validation_tag


class ValidationTagTest(unittest.TestCase):
    # A matching computed label becomes "validated by physicist (Diana, 22 Oct)".
    def test_match(self):
        v = find_validation("CrSe_real", "magnetic space group", "collinear-SOC", "branch")
        self.assertEqual(validation_tag(v, "P6_3'/m'm'c (BNS 194.268)", PENDING_LABEL), "validated by physicist (Diana, 22 Oct)")

    # A different computed label keeps the pending tag and says it differs from the validated value.
    def test_differs(self):
        v = find_validation("MnTe", "magnetic space group", "realization-1-SOC", "branch")
        tag = validation_tag(v, "Cmcm (BNS 63.457)", PENDING_LABEL)
        self.assertTrue(tag.startswith(PENDING_LABEL) and "differs from validated value Cm'c'm (BNS 63.462)" in tag)

    # No entry leaves the fallback tag unchanged.
    def test_none(self):
        self.assertEqual(validation_tag(find_validation("MnTe", "scf convergence", "9-scf", "run"), "reached", "parsed"), "parsed")


if __name__ == "__main__":
    unittest.main()
