"""The chat's read_file_excerpt tool reads small text files of an indexed run and refuses everything else."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from altermagagent import agent


class ReadFileExcerptTest(unittest.TestCase):
    # Build a fake project root with one run folder and patch the index lookups to point at it.
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name) / "proj"
        run = root / "soc" / "bands"
        run.mkdir(parents=True)
        (run / "KPOINTS").write_text("path G-M\n40\nLine-mode\nReciprocal\n")
        (run / "INCAR").write_text("\n".join(f"TAG{i} = {i}" for i in range(200)) + "\n")
        (run / "CHGCAR").write_text("charge\n")
        (run / "blob.dat").write_bytes(b"\x00\x01binary")
        (run / "big.txt").write_bytes(b"a\n" * (1024 * 1024 + 10))
        (root / "secret.txt").write_text("outside the run\n")
        runs = [{"name": "soc/bands", "path": str(run), "details": {"files": []}}]
        self.patches = [mock.patch.object(agent.db, "get_runs", return_value=runs),
                        mock.patch.object(agent.db, "list_projects", return_value=[{"name": "p", "root": str(root)}])]
        for p in self.patches:
            p.start()

    # Stop the patches and remove the fake project.
    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()

    # Call the tool and decode its JSON result.
    def read(self, name, run="soc/bands", **kw):
        return json.loads(agent.read_file_excerpt("p", run, name, **kw))

    # A small KPOINTS file comes back whole, also when the run is named by its last part.
    def test_small_file(self):
        res = self.read("KPOINTS", run="bands")
        self.assertIn("Line-mode", res["text"])
        self.assertEqual(res["run"], "soc/bands")

    # A long file returns at most 60 lines: the first and the last ones.
    def test_line_cap(self):
        lines = [l for l in self.read("INCAR", max_lines=500)["text"].splitlines() if l.startswith("TAG")]
        self.assertEqual(len(lines), 60)
        self.assertEqual(lines[0], "TAG0 = 0")
        self.assertEqual(lines[-1], "TAG199 = 199")

    # CHGCAR, binary files, files over 2 MB and paths outside the run folder are refused.
    def test_refusals(self):
        for name in ("CHGCAR", "blob.dat", "big.txt", "../secret.txt", "/etc/passwd"):
            self.assertTrue(self.read(name)["error"].startswith("refused"), name)

    # An unknown run gives an error naming the indexed runs.
    def test_unknown_run(self):
        self.assertIn("soc/bands", self.read("KPOINTS", run="nope")["error"])


if __name__ == "__main__":
    unittest.main()
