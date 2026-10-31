"""K-path distances, segment breaks and high-symmetry labels for band plots (stage 2 of the specification).

Distance convention: the parser's rec_basis is VASP's reciprocal basis in 1/Angstrom without the factor 2 pi.
Here it is multiplied by 2 pi, so b_i . a_j = 2 pi delta_ij and distances are in 1/Angstrom as in pymatgen and VASPKIT.
Segment boundaries come from a line-mode KPOINTS file next to vasprun.xml, else the vasprun's own line-mode divisions,
else exactly repeated consecutive k-points, else the path is one continuous segment; no distance is added across a boundary.
Labels come from --high-symmetry, else KLABELS, else <kpoints_labels>, else KPOINTS comments, else none.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from altermagagent.bands import build_ticks, klabels_ticks, pretty_label, read_klabels

from .options import KPathError, split_list
from .parser import VasprunData

TWO_PI = 2.0 * math.pi
# KPOINTS corners are typed by hand with fewer decimals than vasprun.xml prints, so they are compared within this.
CORNER_TOL = 1e-5
# Consecutive k-points closer than this (fractional units) are one repeated point that ends a segment.
REPEAT_TOL = 1e-8


@dataclass
class KpointsLineMode:
    """A line-mode KPOINTS file: points per segment, coordinate kind and the corner points with their labels."""

    path: Path
    divisions: int
    reciprocal: bool
    corners: list[tuple[np.ndarray, str]]

    # Number of straight segments the corners describe.
    @property
    def nsegments(self) -> int:
        return len(self.corners) // 2


@dataclass
class KPath:
    """Horizontal axis of a band plot: distances, drawable branches, ticks, and where each piece of it came from."""

    x: np.ndarray
    branches: list[tuple[int, int]]
    segment_starts: list[int]
    ticks: list[tuple[float, str]]
    segment_source: str
    label_source: str
    notes: list[str] = field(default_factory=list)


# Read a line-mode KPOINTS file; return None with a note when it is absent, not line mode or unreadable.
def read_kpoints_line_mode(path: Path) -> tuple[KpointsLineMode | None, str | None]:
    if not path.is_file():
        return None, None
    try:
        lines = path.read_text(errors="replace").splitlines()
        divisions = int(lines[1].split()[0])
        if not lines[2].strip().lower().startswith("l"):
            return None, f"{path.name} is not in line mode; its labels are not used"
        reciprocal = lines[3].strip().lower()[:1] not in ("c", "k")
        corners = []
        for line in lines[4:]:
            bits = line.split("!")[0].split()
            if len(bits) < 3:
                continue
            coords = np.array([float(b) for b in bits[:3]])
            comment = line.split("!", 1)[1] if "!" in line else " ".join(bits[3:])
            corners.append((coords, comment.strip()))
    except (IndexError, ValueError, OSError) as error:
        return None, f"{path.name} could not be read as a line-mode KPOINTS file ({error})"
    if len(corners) % 2 or not corners or divisions < 2:
        return None, f"{path.name} has {len(corners)} corner lines and {divisions} points per segment; it is not used"
    return KpointsLineMode(path=path, divisions=divisions, reciprocal=reciprocal, corners=corners), None


# Cumulative distance in 1/Angstrom (2 pi convention) along the k-point list; no distance is added before a break index.
def cumulative_distance(kpoints: np.ndarray, rec_basis: np.ndarray, breaks: Sequence[int] = ()) -> np.ndarray:
    cartesian = TWO_PI * np.asarray(kpoints, dtype=float) @ np.asarray(rec_basis, dtype=float)
    steps = np.linalg.norm(np.diff(cartesian, axis=0), axis=1)
    # Indices below 1 are ignored: index 0 starts the path anyway, and index - 1 would wrap round to the last step.
    index = np.asarray(breaks, dtype=int)
    steps[index[index >= 1] - 1] = 0.0
    return np.concatenate([[0.0], np.cumsum(steps)])


# True when every corner of a reciprocal line-mode KPOINTS file sits at the matching segment end of the k-point list.
def _corners_match(kpoints_file: KpointsLineMode, k: np.ndarray) -> bool:
    n = kpoints_file.divisions
    ends = np.array([j * n + (n - 1 if end else 0) for j in range(kpoints_file.nsegments) for end in (0, 1)])
    corners = np.array([coords for coords, _ in kpoints_file.corners])
    return bool(np.allclose(corners, k[ends], atol=CORNER_TOL, rtol=0.0))


# Segment start indices (0-based, excluding 0), their source, and the KPOINTS file when it was the accepted source.
def find_segments(data: VasprunData, kpoints_file: KpointsLineMode | None,
                  notes: list[str]) -> tuple[list[int], str, KpointsLineMode | None]:
    nk, k = data.nkpoints, data.kpoints
    if kpoints_file is not None:
        n = kpoints_file.divisions
        if nk != kpoints_file.nsegments * n:
            notes.append(f"{kpoints_file.path.name} describes {kpoints_file.nsegments} segments of {n} points, "
                         f"but the file has {nk} k-points; it is not used")
        elif kpoints_file.reciprocal and not _corners_match(kpoints_file, k):
            notes.append(f"the corners of {kpoints_file.path.name} do not match the k-points of the vasprun; it is not used")
        else:
            return list(range(n, nk, n)), kpoints_file.path.name, kpoints_file
    n = data.kpoint_divisions
    if n and n >= 2 and nk % n == 0 and (data.kpoint_generation == "listgenerated" or labels_from_vasprun(data)):
        return list(range(n, nk, n)), "vasprun line-mode divisions", None
    starts = (np.flatnonzero(np.all(np.abs(np.diff(k, axis=0)) <= REPEAT_TOL, axis=1)) + 1).tolist()
    if starts:
        return starts, "repeated k-points", None
    notes.append("no KPOINTS file, line-mode divisions or repeated k-points mark segment ends; the path is drawn as one line")
    return [], "none (one continuous segment)", None


# Parse --high-symmetry "1:G,41:M" into (0-based index, label) pairs, checking each 1-based number against the k-point count.
def parse_label_option(text: str, nkpoints: int) -> list[tuple[int, str]]:
    pairs = []
    for item in split_list(text):
        match = re.fullmatch(r"(\d+)\s*:\s*(\S.*)", item)
        if not match:
            raise KPathError(f"label {item!r} must have the form NUMBER:NAME, for example 1:G")
        number = int(match.group(1))
        if not 1 <= number <= nkpoints:
            raise KPathError(f"label {item!r}: k-point {number} is outside 1..{nkpoints} (numbering starts at 1)")
        pairs.append((number - 1, match.group(2).strip()))
    if not pairs:
        raise KPathError("--high-symmetry is empty")
    return pairs


# K-point index (0-based) of corner line n (1-based) of a line-mode path with the given points per segment.
def corner_to_index(corner: int, divisions: int) -> int:
    segment, end = divmod(corner - 1, 2)
    return segment * divisions + (divisions - 1 if end else 0)


# (index, label) pairs from <kpoints_labels> corner numbers, or None when they do not fit the k-point list.
def labels_from_vasprun(data: VasprunData) -> list[tuple[int, str]] | None:
    # The segment count is nkpoints / divisions, so the block may list only the corners that carry a label.
    entries, n, nk = data.kpoint_label_entries, data.kpoint_divisions, data.nkpoints
    if not entries or not n or nk % n:
        return None
    if not all(1 <= number <= 2 * (nk // n) for number, _ in entries):
        return None
    return [(corner_to_index(number, n), name) for number, name in entries if name] or None


# (index, label) pairs from the corner comments of an accepted line-mode KPOINTS file, or None without labels.
def labels_from_kpoints(kpoints_file: KpointsLineMode) -> list[tuple[int, str]] | None:
    n = kpoints_file.divisions
    pairs = [(corner_to_index(i, n), label) for i, (_, label) in enumerate(kpoints_file.corners, start=1) if label]
    return pairs or None


# Tick positions and printed labels for (index, label) pairs; labels sharing one position are joined as "A|B".
def ticks_from_indices(pairs: Sequence[tuple[int, str]] | None, x: np.ndarray) -> list[tuple[float, str]] | None:
    if not pairs:
        return None
    return build_ticks([(float(x[i]), pretty_label(label, keep_numbers=True)) for i, label in pairs], float(x[-1]))


# Ticks from a KLABELS file, snapped to the nearest segment end when within 1 percent of the path length.
def ticks_from_klabels(path: Path, x: np.ndarray, starts: Sequence[int], kpoints_file: KpointsLineMode | None,
                       notes: list[str]) -> list[tuple[float, str]] | None:
    rows = read_klabels(path) if path.is_file() else None
    if not rows:
        return None
    length = float(x[-1])
    entries = []
    if kpoints_file is not None:
        entries = [(float(x[corner_to_index(i, kpoints_file.divisions)]), label)
                   for i, (_, label) in enumerate(kpoints_file.corners, start=1)]
    ticks = klabels_ticks(rows, length, entries, notes)
    ends = np.array(sorted({0.0, length} | {float(x[s]) for s in starts}))
    snapped = []
    for position, label in ticks:
        nearest = float(ends[np.argmin(np.abs(ends - position))])
        snapped.append((nearest if abs(nearest - position) <= 0.01 * length else position, label))
    return snapped


# Build the k-path of a parsed file: distances, branches and ticks, using files next to the vasprun when present.
def build_kpath(data: VasprunData, labels: str | None = None, folder: Path | None = None) -> KPath:
    notes: list[str] = []
    folder = Path(data.path).parent if folder is None else folder
    kpoints_file, note = read_kpoints_line_mode(folder / "KPOINTS")
    if note:
        notes.append(note)
    starts, segment_source, accepted = find_segments(data, kpoints_file, notes)
    if data.nkpoints < 2:
        raise KPathError("the k-point list has fewer than two points, so no band path can be drawn")
    x = cumulative_distance(data.kpoints, data.reciprocal_lattice, starts)
    if x[-1] <= 0:
        raise KPathError("the k-point list has zero length, so no band path can be drawn")
    bounds = [0, *starts, data.nkpoints]
    branches = [(a, b - 1) for a, b in zip(bounds[:-1], bounds[1:]) if b - 1 > a]
    resolvers: list[tuple[str, Callable[[], list[tuple[float, str]] | None]]] = [
        ("--high-symmetry option", lambda: ticks_from_indices(parse_label_option(labels, data.nkpoints), x) if labels else None),
        ("KLABELS", lambda: ticks_from_klabels(folder / "KLABELS", x, starts, accepted, notes)),
        ("vasprun <kpoints_labels>", lambda: ticks_from_indices(labels_from_vasprun(data), x)),
        (f"{accepted.path.name} comments" if accepted else "KPOINTS comments",
         lambda: ticks_from_indices(labels_from_kpoints(accepted), x) if accepted else None),
    ]
    found = ((name, ticks) for name, resolve in resolvers if (ticks := resolve()))
    label_source, ticks = next(found, ("none", [(float(x[s]), "") for s in [0, *starts]] + [(float(x[-1]), "")]))
    return KPath(x=x, branches=branches, segment_starts=starts, ticks=ticks, segment_source=segment_source,
                 label_source=label_source, notes=notes)


# Split a --label-map text at commas outside $...$ mathtext and outside braces, so "X=X$_{1,2}$" stays one entry.
def split_label_map(text: str) -> list[str]:
    items, current, in_math, depth = [], [], False, 0
    for char in text:
        if char == "$":
            in_math = not in_math
        elif char == "{":
            depth += 1
        elif char == "}":
            depth = max(depth - 1, 0)
        if char == "," and not in_math and depth == 0:
            items.append("".join(current))
            current = []
        else:
            current.append(char)
    items.append("".join(current))
    return [item.strip() for item in items if item.strip()]


# Parse --label-map "G=Γ,C0=C$_0$" into (key, new label) pairs; a value may use Matplotlib mathtext between $ signs.
def parse_label_map(text: str) -> list[tuple[str, str]]:
    pairs = []
    for item in split_label_map(text):
        key, sep, value = item.partition("=")
        if not sep or not key.strip() or not value.strip():
            raise KPathError(f"--label-map entry {item!r} must have the form OLD=NEW, for example G=Γ")
        pairs.append((key.strip(), value.strip()))
    if not pairs:
        raise KPathError("--label-map is empty")
    for _, value in pairs:
        _check_mathtext(value)
    return pairs


# Raise KPathError when a label's $...$ parts are not valid Matplotlib mathtext, so the error comes before drawing.
def _check_mathtext(value: str) -> None:
    if "$" not in value:
        return
    if value.count("$") % 2:
        raise KPathError(f"label {value!r} has an odd number of $ signs; mathtext needs them in pairs, as in C$_0$")
    from matplotlib.mathtext import MathTextParser  # only here, so the k-path module stays free of Matplotlib otherwise
    try:
        MathTextParser("path").parse(value)
    except ValueError as error:
        raise KPathError(f"label {value!r} is not valid mathtext ({str(error).splitlines()[0]})") from None


# True when a printed label part is named by a --label-map key, written as printed ("Γ") or as typed ("G", "Gamma").
def _label_matches(part: str, key: str) -> bool:
    return part == key or part == pretty_label(key, keep_numbers=True)


# The k-path with its tick labels renamed by --label-map pairs; a key that names no label raises KPathError.
def rename_labels(kpath: KPath, pairs: Sequence[tuple[str, str]]) -> KPath:
    present = [part for _, label in kpath.ticks for part in label.split("|") if part]
    unused = [key for key, _ in pairs if not any(_label_matches(part, key) for part in present)]
    if unused:
        raise KPathError(f"--label-map names {', '.join(unused)}, which is not on the path; "
                         f"labels on the path: {', '.join(dict.fromkeys(present)) or 'none'}")

    # New text of one label part: the value of the first matching key, else the part unchanged.
    def renamed(part: str) -> str:
        return next((value for key, value in pairs if _label_matches(part, key)), part)

    ticks = [(position, "|".join(renamed(part) for part in label.split("|"))) for position, label in kpath.ticks]
    return replace(kpath, ticks=ticks, label_source=f"{kpath.label_source}, renamed by --label-map")


# The k-path with every tick label removed; the vertical lines at the ticks stay.
def without_labels(kpath: KPath) -> KPath:
    return replace(kpath, ticks=[(position, "") for position, _ in kpath.ticks], label_source="none (--no-labels)")


# 0-based k-points kept by --kpoint-stride: every stride-th point of each segment from its start, plus each segment's end.
def stride_indices(nkpoints: int, segment_starts: Sequence[int], stride: int) -> np.ndarray:
    if stride < 1:
        raise KPathError(f"--kpoint-stride must be 1 or more, not {stride}")
    bounds = [0, *segment_starts, nkpoints]
    kept: set[int] = set()
    for start, stop in zip(bounds[:-1], bounds[1:]):
        if stop > start:
            kept.update(range(start, stop, stride))
            kept.add(stop - 1)
    return np.array(sorted(kept), dtype=int)


# The k-path restricted to the kept k-points (from stride_indices), with branches and segment starts renumbered.
def subsample(kpath: KPath, kept: np.ndarray) -> KPath:
    position = {int(k): i for i, k in enumerate(kept)}
    branches = [(position[a], position[b]) for a, b in kpath.branches]
    starts = [position[s] for s in kpath.segment_starts]
    return replace(kpath, x=kpath.x[kept], branches=branches, segment_starts=starts)
