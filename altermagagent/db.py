"""SQLite index stored under ~/.altermagagent, never inside a scanned root."""

import json
import os
import sqlite3
from pathlib import Path

HOME_DIR = Path(os.environ.get("ALTERMAG_HOME", Path.home() / ".altermagagent"))
DB_PATH = HOME_DIR / "index.sqlite"

SCHEMA = """
CREATE TABLE IF NOT EXISTS projects (
    name TEXT PRIMARY KEY,
    root TEXT NOT NULL,
    imported_at TEXT NOT NULL,
    summary_json TEXT,
    details_json TEXT
);
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    project TEXT NOT NULL REFERENCES projects(name) ON DELETE CASCADE,
    seq INTEGER NOT NULL,
    name TEXT NOT NULL,
    path TEXT NOT NULL,
    status TEXT NOT NULL,
    first_time REAL,
    last_time REAL,
    parent TEXT,
    parent_note TEXT,
    details_json TEXT
);
CREATE INDEX IF NOT EXISTS runs_project ON runs(project, seq);
CREATE TABLE IF NOT EXISTS band_match (
    project TEXT NOT NULL,
    run TEXT NOT NULL,
    key TEXT NOT NULL,
    result_json TEXT NOT NULL,
    PRIMARY KEY (project, run)
);
CREATE TABLE IF NOT EXISTS parent_choices (
    project TEXT NOT NULL,
    run TEXT NOT NULL,
    parent TEXT,
    PRIMARY KEY (project, run)
);
"""
UNKNOWN_PARENT = ""


# Open the index database, creating the folder and tables when needed.
def connect():
    HOME_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB_PATH)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    con.executescript(SCHEMA)
    return con


# Replace one project and all its runs in a single transaction.
def save_project(name, root, imported_at, runs, summary=None, details=None):
    con = connect()
    with con:
        con.execute("DELETE FROM runs WHERE project = ?", (name,))
        con.execute("DELETE FROM projects WHERE name = ?", (name,))
        con.execute("INSERT INTO projects VALUES (?, ?, ?, ?, ?)",
                    (name, root, imported_at, json.dumps(summary), json.dumps(details or {})))
        for r in runs:
            con.execute(
                "INSERT INTO runs (project, seq, name, path, status, first_time, last_time, parent, parent_note, details_json)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (name, r["seq"], r["name"], r["path"], r["status"], r["first_time"], r["last_time"],
                 r.get("parent"), r.get("parent_note"), json.dumps(r["details"])))
    con.close()


# Store a freshly built summary for an existing project.
def save_summary(name, summary):
    con = connect()
    with con:
        con.execute("UPDATE projects SET summary_json = ? WHERE name = ?", (json.dumps(summary), name))
    con.close()


# Return all projects as plain dicts with run counts.
def list_projects():
    con = connect()
    rows = con.execute(
        "SELECT p.name, p.root, p.imported_at, COUNT(r.id) AS n_runs FROM projects p"
        " LEFT JOIN runs r ON r.project = p.name GROUP BY p.name ORDER BY p.name").fetchall()
    con.close()
    return [dict(r) for r in rows]


# Return one project dict with its summary, or None.
def get_project(name):
    con = connect()
    row = con.execute("SELECT * FROM projects WHERE name = ?", (name,)).fetchone()
    con.close()
    if row is None:
        return None
    d = dict(row)
    d["summary"] = json.loads(d.pop("summary_json") or "null")
    d["details"] = json.loads(d.pop("details_json") or "{}")
    return d


# Return the runs of a project in historical order, details decoded.
def get_runs(name):
    con = connect()
    rows = con.execute("SELECT * FROM runs WHERE project = ? ORDER BY seq", (name,)).fetchall()
    con.close()
    out = []
    choices = parent_choices(name)
    for row in rows:
        d = dict(row)
        d["details"] = json.loads(d.pop("details_json") or "{}")
        d["parent_auto"], d["parent_auto_note"], d["parent_user"] = d["parent"], d["parent_note"], d["name"] in choices
        if d["parent_user"]:
            auto = f'{d["parent"]} ({d["parent_auto_note"]})' if d["parent"] else d["parent_auto_note"]
            d["parent"] = choices[d["name"]] or None
            d["parent_note"] = f"set by the user in the interface; automatic: {auto}"
        out.append(d)
    return out


# Return the user's parent choices of a project as {run: parent or ""}; "" means "parent unknown".
def parent_choices(name):
    con = connect()
    rows = con.execute("SELECT run, parent FROM parent_choices WHERE project = ?", (name,)).fetchall()
    con.close()
    return {r["run"]: r["parent"] or UNKNOWN_PARENT for r in rows}


# Store the user's parent choice for one run (None removes it, "" means unknown); survives a re-import.
def set_parent_choice(project, run, parent):
    con = connect()
    with con:
        if parent is None:
            con.execute("DELETE FROM parent_choices WHERE project = ? AND run = ?", (project, run))
        else:
            con.execute("INSERT OR REPLACE INTO parent_choices VALUES (?, ?, ?)", (project, run, parent))
    con.close()


# Cached Wannier band-match result of one run when its cache key still matches, else None.
def get_band_match(project, run, key):
    con = connect()
    row = con.execute("SELECT key, result_json FROM band_match WHERE project = ? AND run = ?", (project, run)).fetchone()
    con.close()
    return json.loads(row["result_json"]) if row and row["key"] == key else None


# Store the Wannier band-match result of one run in the index (never in a data folder).
def set_band_match(project, run, key, result):
    con = connect()
    with con:
        con.execute("INSERT OR REPLACE INTO band_match VALUES (?, ?, ?, ?)", (project, run, key, json.dumps(result)))
    con.close()
