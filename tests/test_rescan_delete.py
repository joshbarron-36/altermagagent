"""Rescan and delete of indexed projects: the index changes, the data folders are only read."""

import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

from _fixtures import demo_tree, temp_index
from altermagagent import cli, db, importer


class RescanDelete(unittest.TestCase):
    """Import the demo tree into a temporary index, then rescan and delete projects."""

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
        temp_index(self, self.tmp.name)

    # Count the rows of one table that belong to a project.
    def rows(self, table, name):
        con = db.connect()
        n = con.execute(f"SELECT COUNT(*) FROM {table} WHERE project = ?", (name,)).fetchone()[0]
        con.close()
        return n

    # Rescan keeps the run count, rebuilds the summary and keeps the user's parent choices.
    def test_rescan_rebuilds_summary(self):
        _, n = importer.import_project(self.demo, "CrSe")
        runs = db.get_runs("CrSe")
        db.set_parent_choice("CrSe", runs[-1]["name"], "")
        db.save_summary("CrSe", {"points": [], "gaps": [], "branches": []})
        name, n2 = importer.rescan_project("CrSe")
        p = db.get_project("CrSe")
        self.assertEqual((name, n2), ("CrSe", n))
        self.assertEqual(len(db.get_runs("CrSe")), n)
        self.assertGreaterEqual(len(p["summary"]["points"]), 8)
        self.assertEqual(p["root"], str(self.demo.resolve()))
        self.assertEqual(db.parent_choices("CrSe"), {runs[-1]["name"]: ""})

    # Rescan of a project whose folder is gone fails and leaves the stored project unchanged.
    def test_rescan_missing_root(self):
        copy = Path(self.tmp.name) / "gone"
        shutil.copytree(self.demo, copy)
        _, n = importer.import_project(copy, "gone")
        before = db.get_project("gone")
        shutil.rmtree(copy)
        with self.assertRaises(FileNotFoundError):
            importer.rescan_project("gone")
        self.assertEqual(db.get_project("gone"), before)
        self.assertEqual(len(db.get_runs("gone")), n)
        with patch("sys.stderr"):
            self.assertEqual(cli.main(["rescan", "gone"]), 1)

    # Rescan whose summary build fails leaves the stored project, its summary and its runs as they were.
    def test_rescan_failed_summary_keeps_index(self):
        _, n = importer.import_project(self.demo, "CrSe")
        before = db.get_project("CrSe")
        with patch("altermagagent.summary.build", side_effect=RuntimeError("boom")):
            with self.assertRaises(RuntimeError):
                importer.rescan_project("CrSe")
        after = db.get_project("CrSe")
        self.assertTrue(after["summary"]["points"])
        self.assertEqual(after["summary"], before["summary"])
        self.assertEqual(after["imported_at"], before["imported_at"])
        self.assertEqual(len(db.get_runs("CrSe")), n)

    # Delete removes the project and its runs, band matches and parent choices, and leaves other projects alone.
    def test_delete_removes_rows(self):
        _, n = importer.import_project(self.demo, "CrSe")
        importer.import_project(self.demo, "other")
        run = db.get_runs("CrSe")[0]["name"]
        db.set_parent_choice("CrSe", run, "")
        db.set_band_match("CrSe", run, "k", {"ok": True})
        self.assertEqual(db.delete_project("CrSe"), n)
        self.assertIsNone(db.get_project("CrSe"))
        for table in db.PROJECT_TABLES + ("runs",):
            self.assertEqual(self.rows(table, "CrSe"), 0, table)
        self.assertEqual(len(db.get_runs("other")), n)
        self.assertTrue(self.demo.is_dir())

    # Delete of an unknown name returns None, and the command exits 1.
    def test_delete_unknown(self):
        self.assertIsNone(db.delete_project("nothing"))
        with patch("sys.stderr"):
            self.assertEqual(cli.main(["delete", "nothing"]), 1)
        with self.assertRaises(ValueError):
            importer.rescan_project("nothing")


if __name__ == "__main__":
    unittest.main()
