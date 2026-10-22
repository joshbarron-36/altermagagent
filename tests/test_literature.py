"""Literature line: never raises, says "literature unavailable" without a token or on 401, caches results, keeps its header."""

import json
import os
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest import mock

from altermagagent import literature as L

PAPER = {"title": "A paper", "authors": ["A. Author", "B. Author"], "source": "arxiv", "url": "https://arxiv.org/abs/1",
         "citation_count": 5, "dates": {"first_available": {"value": "2024-02-01"}}}
RUNS = [{"details": {"formula": "MnTe", "spacegroup": "P6_3/mmc"}}]
SUMMARY = {"diagnosis": {"not_supported": ["non-relativistic spin splitting", "magnetic space group classification"]}}


class FakeResponse:
    """Stand-in for a urllib response with JSON content."""

    # Hold one reply dict.
    def __init__(self, reply):
        self.body = json.dumps(reply).encode()
        self.headers = {"Content-Type": "application/json"}

    # Return the body bytes.
    def read(self):
        return self.body

    # Enter the context manager.
    def __enter__(self):
        return self

    # Leave the context manager.
    def __exit__(self, *a):
        return False


# Fake urlopen that answers initialize, the notification and tools/call with one paper.
def fake_urlopen(req, timeout=None):
    msg = json.loads(req.data)
    if msg.get("method") == "tools/call":
        text = json.dumps({"papers": [PAPER]})
        return FakeResponse({"jsonrpc": "2.0", "id": msg["id"], "result": {"content": [{"type": "text", "text": text}]}})
    if "id" not in msg:
        return FakeResponse({})
    return FakeResponse({"jsonrpc": "2.0", "id": msg["id"], "result": {"protocolVersion": "2025-06-18"}})


class Literature(unittest.TestCase):
    """Each test uses a temporary cache file."""

    # Point the cache at a temporary file.
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cache = mock.patch.object(L, "CACHE_PATH", Path(self.tmp.name) / "literature_cache.json")
        self.cache.start()

    # Restore the cache path and remove the temporary folder.
    def tearDown(self):
        self.cache.stop()
        self.tmp.cleanup()

    # Without a token the search returns an empty list with the note, and the section says unavailable.
    def test_no_token(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            res = L.search_literature("MnTe")
            sec = L.section(SUMMARY, RUNS)
        self.assertEqual(list(res), [])
        self.assertTrue(res.note.startswith("literature unavailable"))
        text = "\n".join(L.section_lines(sec))
        self.assertIn("From the literature, not from these files", text)
        self.assertIn("not evidence from the project folder", text)
        self.assertIn("literature unavailable", text)

    # A 401 reply gives the note instead of an exception.
    def test_401(self):
        err = urllib.error.HTTPError(L.MCP_URL, 401, "Unauthorized", {}, None)
        with mock.patch.dict(os.environ, {"VALENCY_TOKEN": "x"}), mock.patch("urllib.request.urlopen", side_effect=err):
            res = L.search_literature("MnTe")
        self.assertEqual(list(res), [])
        self.assertIn("401", res.note)

    # A working search returns slim papers, and a second call within 24 hours comes from the cache.
    def test_search_and_cache(self):
        with mock.patch.dict(os.environ, {"VALENCY_TOKEN": "x"}), \
                mock.patch("urllib.request.urlopen", side_effect=fake_urlopen) as uo:
            first = L.search_literature("MnTe altermagnet spin splitting")
            calls = uo.call_count
            second = L.search_literature("MnTe altermagnet spin splitting")
            self.assertEqual(uo.call_count, calls)
        self.assertEqual(first[0], {"title": "A paper", "first_author": "A. Author", "year": "2024", "source": "arxiv",
                                    "url": "https://arxiv.org/abs/1", "citation_count": 5})
        self.assertEqual(list(second), list(first))

    # Two queries for MnTe: spin splitting and magnetic order with the space group.
    def test_queries(self):
        qs = [q for q, _ in L.queries("MnTe", "P6_3/mmc", SUMMARY["diagnosis"]["not_supported"])]
        self.assertEqual(qs, ["MnTe altermagnet spin splitting", "MnTe magnetic order P6_3/mmc"])


if __name__ == "__main__":
    unittest.main()
