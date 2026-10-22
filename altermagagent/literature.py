"""Published papers found by a Valency search: a line kept apart from the evidence read from the project folders."""

import json
import os
import time
import urllib.error
import urllib.request

from . import db

MCP_URL = "https://labs.valency.io/mcp"
TIMEOUT_S = 10
CACHE_PATH = db.HOME_DIR / "literature_cache.json"
CACHE_SECONDS = 24 * 3600
HEADER = "From the literature, not from these files"
HEADER_NOTE = ("These are published papers found by a Valency search, not evidence from the project folder; "
               "the suggestions above do not change because of them.")
PROTOCOL = "2025-06-18"


class LiteratureUnavailable(Exception):
    """Raised inside this module only; callers get an empty Papers list with a note instead."""


class Papers(list):
    """A list of paper dicts with a note: empty string when the search worked, else "literature unavailable: ..."."""

    # Build the list from papers and a note.
    def __init__(self, papers=(), note=""):
        super().__init__(papers)
        self.note = note


# Parse an MCP response body that is either plain JSON or a text/event-stream with "data:" lines.
def parse_mcp(content_type, body):
    if "text/event-stream" in (content_type or ""):
        for line in body.splitlines():
            if line.startswith("data:"):
                return json.loads(line[5:].strip())
        raise LiteratureUnavailable("event stream had no data line")
    return json.loads(body)


class Client:
    """A minimal MCP client over HTTP POST with JSON-RPC 2.0, using urllib only."""

    # Keep the token and the headers; nothing is sent until initialize.
    def __init__(self, token, url=MCP_URL, timeout=TIMEOUT_S):
        self.url, self.timeout, self.next_id = url, timeout, 1
        self.headers = {"Authorization": "Bearer " + token, "Content-Type": "application/json",
                        "Accept": "application/json, text/event-stream"}

    # Send one JSON-RPC message; return the parsed reply, or None for a notification.
    def post(self, msg):
        req = urllib.request.Request(self.url, data=json.dumps(msg).encode(), headers=self.headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                ctype, body, sid = resp.headers.get("Content-Type"), resp.read().decode("utf-8", "replace"), \
                    resp.headers.get("Mcp-Session-Id")
        except urllib.error.HTTPError as exc:
            if exc.code == 401:
                raise LiteratureUnavailable("the Valency token was rejected (HTTP 401); renew it with valency_login.py --refresh")
            raise LiteratureUnavailable(f"HTTP {exc.code} from {self.url}")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise LiteratureUnavailable(f"no answer from {self.url} within {self.timeout} s ({type(exc).__name__})")
        if sid:
            self.headers["Mcp-Session-Id"] = sid
        if "id" not in msg or not body.strip():
            return None
        reply = parse_mcp(ctype, body)
        if "error" in reply:
            raise LiteratureUnavailable(f"MCP error: {reply['error'].get('message', reply['error'])}")
        return reply.get("result", {})

    # Send one request with a fresh id and return its result.
    def request(self, method, params=None):
        msg = {"jsonrpc": "2.0", "id": self.next_id, "method": method}
        self.next_id += 1
        if params is not None:
            msg["params"] = params
        return self.post(msg)

    # Open the MCP session: initialize, remember the protocol version, then send notifications/initialized.
    def initialize(self):
        res = self.request("initialize", {"protocolVersion": PROTOCOL, "capabilities": {},
                                          "clientInfo": {"name": "AltermagAgent", "version": "1.0"}})
        self.headers["MCP-Protocol-Version"] = res.get("protocolVersion", PROTOCOL)
        self.post({"jsonrpc": "2.0", "method": "notifications/initialized"})
        return res

    # Call one tool and return the JSON decoded from its first text content block.
    def call_tool(self, name, arguments):
        res = self.request("tools/call", {"name": name, "arguments": arguments})
        if res.get("isError"):
            raise LiteratureUnavailable(f"tool {name} reported an error")
        text = next((c.get("text", "") for c in res.get("content", []) if c.get("type") == "text"), "")
        try:
            return json.loads(text)
        except ValueError:
            raise LiteratureUnavailable(f"tool {name} returned no JSON")


# Reduce one Valency paper object to title, first author, year, source, url and citation count.
def slim(p):
    dates = p.get("dates") or {}
    first = (dates.get("first_available") or {}).get("value") or p.get("first_submitted") or p.get("datestamp") or ""
    authors = p.get("authors") or []
    return {"title": " ".join((p.get("title") or "").split()), "first_author": authors[0] if authors else "",
            "year": first[:4], "source": p.get("source", ""), "url": p.get("url", ""),
            "citation_count": p.get("citation_count")}


# Read the query cache, or an empty dict when it is missing or unreadable.
def load_cache():
    try:
        return json.loads(CACHE_PATH.read_text())
    except (OSError, ValueError):
        return {}


# Store the papers of one query in the cache with the current time.
def save_cache(key, papers):
    cache = load_cache()
    cache[key] = {"time": time.time(), "papers": papers}
    try:
        CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        CACHE_PATH.write_text(json.dumps(cache, indent=1, ensure_ascii=False))
    except OSError:
        pass


# Search published papers for a query; returns Papers (possibly empty with a note), never raises.
def search_literature(query, limit=3, client=None):
    token = os.environ.get("VALENCY_TOKEN", "").strip()
    if not token:
        return Papers(note="literature unavailable: VALENCY_TOKEN is not set")
    key = f"{query}|{limit}"
    hit = load_cache().get(key)
    if hit and time.time() - hit.get("time", 0) < CACHE_SECONDS:
        return Papers(hit["papers"])
    try:
        if client is None:
            client = Client(token)
            client.initialize()
        res = client.call_tool("semantic_search_papers", {"query": query, "limit": limit, "include_abstract": False})
        papers = [slim(p) for p in (res.get("papers") or [])[:limit]]
    except LiteratureUnavailable as exc:
        return Papers(note=f"literature unavailable: {exc}")
    except Exception as exc:
        return Papers(note=f"literature unavailable: {type(exc).__name__}: {exc}")
    save_cache(key, papers)
    return Papers(papers)


# Formula and space group of the project (or of one branch), from the first run that has them.
def project_material(runs, branch=None):
    if branch:
        names = set(branch.get("runs") or [])
        runs = [r for r in runs if r["name"] in names] or runs
    for r in runs:
        d = r["details"]
        if d.get("formula"):
            return d["formula"], d.get("spacegroup")
    return None, None


# Open questions of one branch: title and evidence of every suggestion step that is not done.
def open_steps(branch):
    out = []
    for st in ((branch or {}).get("suggestions") or {}).get("steps", []):
        if st["state"] in ("missing", "partial"):
            out += [st["title"]] + list(st.get("evidence") or [])
    return out


# One or two queries for the formula; the second only when a "not supported" claim names the magnetic order.
def queries(formula, spacegroup, not_supported):
    text = " ".join(not_supported or []).lower()
    out = [(f"{formula} altermagnet spin splitting", "spin splitting")]
    if "magnetic" in text or "order" in text:
        out.append((f"{formula} magnetic order {spacegroup or ''}".strip(), "magnetic order and its classification"))
    return out


# Build the literature section for one project or branch: header, note and up to three papers per query; never raises.
def section(summary, runs, limit=3, branch=None):
    formula, sg = project_material(runs, branch)
    sec = {"header": HEADER, "note": HEADER_NOTE, "queries": [], "unavailable": ""}
    if not formula:
        sec["unavailable"] = "literature unavailable: no formula was read from the project files"
        return sec
    unsupported = (summary.get("diagnosis") or {}).get("not_supported") or open_steps(branch)
    for q, topic in queries(formula, sg, unsupported):
        papers = search_literature(q, limit=limit)
        if papers.note:
            sec["unavailable"] = papers.note
            break
        sec["queries"].append({"query": q, "topic": topic, "papers": list(papers)})
    return sec


# The literature section as plain text lines for the command line.
def section_lines(sec):
    lines = [sec["header"], f"  {sec['note']}"]
    if sec.get("unavailable"):
        lines.append(f"  {sec['unavailable']}")
    for q in sec["queries"]:
        lines.append(f'  Search "{q["query"]}" (relates to: {q["topic"]}):')
        for p in q["papers"] or []:
            lines.append(f'    - {p["title"]} ({p["first_author"]}, {p["year"]}) {p["url"]}')
        if not q["papers"]:
            lines.append("    - no papers found")
    return lines
