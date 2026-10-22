"""Measure phase, status and parent accuracy of the importer on seeded synthetic trees with ground truth."""

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "demo"))
sys.path.insert(0, str(REPO))


# Turn the importer's parent fields of one run into ("name", x) or ("candidates", list) or ("unknown", None).
def got_parent(run):
    cands = run["details"].get("parent_candidates") or []
    if run["parent"]:
        return run["parent"], None
    if cands:
        return "candidates", cands
    return "unknown", None


# Return True when the importer's parent answer agrees with the truth under the scoring rule.
def parent_ok(truth, label, cands, note=""):
    if truth["label"] == "candidates":
        return label == "candidates" and truth["true_parent"] in cands
    if truth["label"] == "unknown" and truth["true_parent"] and "restart file not present now" in (note or ""):
        # The restart file was deleted: a link or candidate list marked "restart file not present now" counts when it holds the true parent.
        return label == truth["true_parent"] or (label == "candidates" and truth["true_parent"] in cands)
    return label == truth["label"]


# Generate n trees, import each into a temporary index and compare with the truth; returns (scores, mismatches).
def measure(n=20, first_seed=1, base=None):
    tmp = None
    if base is None:
        tmp = tempfile.TemporaryDirectory()
        base = tmp.name
    base = Path(base)
    os.environ["ALTERMAG_HOME"] = str(base / "home")
    import make_variants
    from altermagagent import db, importer
    db.HOME_DIR, db.DB_PATH = base / "home", base / "home" / "index.sqlite"
    scores = {f: [0, 0] for f in ("phase", "status", "parent")}
    mismatches = []
    try:
        for k in range(n):
            root = base / "trees" / f"tree_{k + 1:02d}"
            truth = make_variants.make_tree(root, first_seed + k)
            importer.import_project(root, root.name)
            runs = {r["name"]: r for r in db.get_runs(root.name)}
            for name, t in truth["runs"].items():
                r = runs.get(name)
                for field in ("phase", "status"):
                    got = r["details"].get("phase") if r and field == "phase" else r["status"] if r else "not found"
                    scores[field][1] += 1
                    if got == t[field]:
                        scores[field][0] += 1
                    else:
                        mismatches.append((root.name, name, field, t[field], got))
                label, cands = got_parent(r) if r else ("not found", None)
                scores["parent"][1] += 1
                if parent_ok(t["parent"], label, cands or [], r["parent_note"] if r else ""):
                    scores["parent"][0] += 1
                else:
                    exp = t["parent"]["label"]
                    if exp == "candidates":
                        exp = f"candidates incl. {t['parent']['true_parent']}"
                    mismatches.append((root.name, name, "parent", exp,
                                       label + (f" {cands}" if cands else "") + f" ({r['parent_note']})" if r else label))
    finally:
        if tmp:
            tmp.cleanup()
    return scores, mismatches


# Print per-field correct / total and percentage, then every mismatch.
def report(scores, mismatches):
    for f, (ok, tot) in scores.items():
        print(f"{f:7s} {ok:4d} / {tot:4d}  {100 * ok / tot:5.1f} %")
    print(f"\n{len(mismatches)} mismatches")
    for tree, run, field, exp, got in mismatches:
        print(f"- {tree} | {run} | {field} | expected: {exp} | got: {got}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--keep", help="write the trees and index under this folder instead of a temporary one")
    a = ap.parse_args()
    report(*measure(a.n, a.seed, a.keep))
