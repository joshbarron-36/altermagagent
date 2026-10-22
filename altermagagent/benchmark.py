"""Compare the indexed values of a project with a benchmark JSON file and give pass, fail or unknown per row."""

import json

from .summary import norm, splitting_point

QUANTITIES = ("net_moment", "site_moment", "spacegroup_number", "lattice_a", "lattice_c", "splitting_max",
              "symmops", "total_energy")


# Find one run by name, or None.
def run_by_name(runs, name):
    return next((r for r in runs if r["name"].rstrip("/") == (name or "").rstrip("/")), None)


# Actual value of one benchmark row from the index as (value, source text); value None means unknown.
def actual(row, summary, runs):
    q = row.get("quantity")
    pts = {p["key"]: p for p in summary.get("points", [])}
    mag = (pts.get("magnetic") or {}).get("data") or {}
    run = run_by_name(runs, row.get("run") or mag.get("run")) if q in ("net_moment", "site_moment") else run_by_name(runs, row.get("run"))
    if q in ("spacegroup_number", "lattice_a", "lattice_c") and run is None:
        run = next((r for r in runs if r["details"].get("structure_file")), None)
    if q == "net_moment":
        mom = run and run["details"].get("moments")
        if not mom:
            return None, "no magnetization block in OUTCAR"
        net = norm(mom["net_vector"]) if mom.get("net_vector") else abs(mom["net"])
        return net, f'{run["name"]}/OUTCAR'
    if q == "site_moment":
        mom = run and run["details"].get("moments")
        if not mom:
            return None, "no magnetization block in OUTCAR"
        vals = [abs(s["moment"]) for s in mom["sites"] if s.get("element") == row.get("element")]
        return (max(vals), f'{run["name"]}/OUTCAR') if vals else (None, f'no {row.get("element")} site in {run["name"]}/OUTCAR')
    if q in ("spacegroup_number", "lattice_a", "lattice_c"):
        if run is None:
            return None, "no readable structure"
        key = {"spacegroup_number": "spacegroup_number", "lattice_a": "a", "lattice_c": "c"}[q]
        return run["details"].get(key), f'{run["name"]}/{run["details"].get("structure_file")}'
    if q == "splitting_max":
        sd = splitting_point(runs).get("data") or {}
        if sd.get("split") is None:
            return None, sd.get("reason", "no non-SOC eigenvalues read")
        return sd["split"], f'{sd["run"]}/vasprun.xml, {sd["nk"]} k-points {sd["sampled"]}'
    if q == "symmops":
        for r in runs:
            for e in r["details"].get("steps", {}).get("symmetrisation", []):
                if e.get("symmops"):
                    return e["symmops"], f'{r["name"].rstrip("/")}/{e["folder"]}/mysymmop.dat'
        return None, "no mysymmop.dat with listed operations"
    if q == "total_energy":
        if run is None:
            return None, f'no run named {row.get("run")}'
        return run["details"].get("energy"), f'{run["name"]}/{run["details"].get("energy_source")}'
    return None, f"unknown quantity {q}; known: {', '.join(QUANTITIES)}"


# Compare every row of a benchmark dict and return result dicts with pass, fail or unknown.
def compare(bench, summary, runs):
    out = []
    for row in bench.get("rows", []):
        val, src = actual(row, summary, runs)
        exp, tol = row.get("expected"), row.get("tolerance", 0)
        if exp is None:
            res, why = "unknown", "expected value not filled in"
        elif val is None:
            res, why = "unknown", "value not read from these files"
        else:
            diff = abs(float(val) - float(exp))
            res, why = ("pass" if diff <= float(tol) else "fail"), f"|actual − expected| = {diff:.2g}, tolerance {tol}"
        out.append({"id": row.get("id", row.get("quantity")), "quantity": row.get("quantity"), "expected": exp,
                    "actual": val, "tolerance": tol, "units": row.get("units", ""), "result": res, "reason": why,
                    "source": src})
    return out


# Load a benchmark file and compare it with an indexed project.
def run_benchmark(path, summary, runs):
    with open(path, encoding="utf-8") as fh:
        bench = json.load(fh)
    return compare(bench, summary, runs)


# Format one value for the result table.
def fmt(v):
    if v is None:
        return "—"
    return f"{v:.6g}" if isinstance(v, float) else str(v)


# Render results as plain text lines, one row per benchmark entry, plus a count line.
def as_text(results):
    lines = []
    for r in results:
        lines.append(f'{r["result"].upper():7} {r["id"]}: expected {fmt(r["expected"])}, actual {fmt(r["actual"])} {r["units"]}; '
                     f'{r["reason"]} ({r["source"]})')
    n = {k: sum(1 for r in results if r["result"] == k) for k in ("pass", "fail", "unknown")}
    lines.append(f'{n["pass"]} pass, {n["fail"]} fail, {n["unknown"]} unknown')
    return "\n".join(lines)
