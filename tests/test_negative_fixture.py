"""Negative fixture test: without the bands_nosoc run the per-branch suggestions mark step (d) missing and the splitting row is unknown."""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class NegativeFixture(unittest.TestCase):
    """Build both demo trees in a temporary folder with a temporary index and compare them with the CrSe benchmark."""

    # Write both trees and import them into an index under a temporary ALTERMAG_HOME.
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        base = Path(cls.tmp.name)
        os.environ["ALTERMAG_HOME"] = str(base / "home")
        for flag, out in (([], "full"), (["--without-nosoc-bands"], "negative")):
            subprocess.run([sys.executable, str(REPO / "demo" / "make_demo.py"), *flag, "--out", str(base / out)],
                           check=True, capture_output=True)
        from altermagagent import benchmark, db, importer
        cls.db, cls.bench = db, benchmark
        importer.import_project(base / "full", "full")
        importer.import_project(base / "negative", "negative")

    # Remove the temporary folder.
    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    # Benchmark results of one project as {row id: result}.
    def results(self, name):
        p = self.db.get_project(name)
        res = self.bench.run_benchmark(REPO / "demo" / "benchmark_crse.json", p["summary"], self.db.get_runs(name))
        return {r["id"]: r["result"] for r in res}

    # Suggestions of the non-SOC branch of one project.
    def nosoc_suggestions(self, name):
        brs = self.db.get_project(name)["summary"]["branches"]
        b = next(b for b in brs if not b["soc"])
        return {st["step"]: st for st in b["suggestions"]["steps"]}, b["suggestions"]["next"]

    # Step (b) "define MSG" computes a magnetic space group tagged pending, writes findsym.in, and the Symmetries point shows it.
    def test_define_msg_computed_pending(self):
        steps, _ = self.nosoc_suggestions("full")
        msg = steps["b"]["msg"]
        self.assertTrue(msg["computed"] and "BNS" in msg["computed"])
        self.assertTrue(Path(msg["findsym_in"]).exists())
        self.assertTrue(any("parsed, pending validation" in e for e in steps["b"]["evidence"]))
        sym = next(p for p in self.db.get_project("full")["summary"]["points"] if p["key"] == "symmetries")
        label = msg["computed"].split(" (BNS")[0]
        self.assertIn(label + " (computed)", sym["value"])
        self.assertNotIn("BNS", sym["value"])
        self.assertNotIn("[", sym["value"])
        self.assertEqual(steps["b"]["state"], "done")
        self.assertIsNone(steps["b"]["action"])

    # The negative tree marks step (d), band structure, as missing in the non-SOC branch and names a band run next.
    def test_negative_step_d_missing(self):
        steps, nxt = self.nosoc_suggestions("negative")
        self.assertEqual(steps["d"]["state"], "missing")
        self.assertIn("ICHARG = 11", nxt)

    # The splitting row prints unknown on the negative tree; every other row passes.
    def test_negative_benchmark_unknown(self):
        res = self.results("negative")
        self.assertEqual(res.pop("spin splitting"), "unknown")
        self.assertEqual(set(res.values()), {"pass"})

    # The full tree passes every row and its non-SOC branch has step (d) done from bands_nosoc.
    def test_full_tree_passes(self):
        self.assertEqual(set(self.results("full").values()), {"pass"})
        steps, _ = self.nosoc_suggestions("full")
        self.assertEqual(steps["d"]["state"], "done")
        self.assertTrue(any("bands_nosoc" in e for e in steps["d"]["evidence"]))

    # Both trees have two branches, non-SOC and SOC, of one collinear realization.
    def test_branches(self):
        for name in ("full", "negative"):
            brs = self.db.get_project(name)["summary"]["branches"]
            self.assertEqual(sorted(b["soc"] for b in brs), [False, True])
            self.assertEqual({b["realization"] for b in brs}, {"collinear"})

    # No suggestion line states that a material is an altermagnet, and no summary point is about spin splitting.
    def test_no_altermagnet_claim(self):
        for name in ("full", "negative"):
            summ = self.db.get_project(name)["summary"]
            text = " ".join(" ".join(st["evidence"] + [st.get("action") or ""]) + b["suggestions"]["next"]
                            for b in summ["branches"] for st in b["suggestions"]["steps"]).lower()
            self.assertNotIn("altermagnet", text)
            self.assertNotIn("splitting", {p["key"] for p in summ["points"]})


if __name__ == "__main__":
    unittest.main()
