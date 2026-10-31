"""Stage 1 of the orbital-projected band plotter: a streaming reader of vasprun.xml.

The reader walks the file once with lxml.etree.iterparse and frees every block after reading it.
It never uses PROCAR, EIGENVAL or the energy-integrated <dos><partial> block in place of band projections.

Array conventions of the returned VasprunData (all indices 0-based inside Python):

- eigenvalues, occupations: shape (spin, kpoint, band); spin has 1 entry, or 2 for collinear ISPIN = 2.
- projections: shape (component, kpoint, band, ion, orbital).
- component is ("total",) without spin, ("up", "down") for ISPIN = 2, ("total", "mx", "my", "mz") for noncollinear runs.
- The number of spin or component blocks is counted in the file and must match the labels of the detected layout.
- The <dimension> entries of a band-data <array> may come in any order; each block is transposed to the order above.
- orbital labels come from the <field> entries of <projected><array>, in file order, with padding stripped.
- Ion n in VASP's 1-based numbering is index n - 1; use VasprunData.ion_index for an explicit, checked conversion.
- reciprocal_lattice is VASP's rec_basis in 1/Angstrom without the factor 2 pi; rows are b1, b2, b3.
- kpoint_label_entries are the raw (number, name) pairs of <kpoints_labels>.
- With several <calculation> blocks, band data, Fermi energy and partial-DOS flag all come from the last one alone.
- In the VASP 6.3.2 band files checked, number n is the n-th corner line of a line-mode KPOINTS file, not a k-point index.
- With with_projections=False the <projected> numbers are freed unread: projections is None, has_projected_block says
  whether the block exists, and the orbital labels are still read from its <field> entries.
"""

from __future__ import annotations

import contextlib
import gzip
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any

import numpy as np
from lxml import etree

LAYOUT_NONMAGNETIC = "nonmagnetic"
LAYOUT_COLLINEAR = "collinear"
LAYOUT_NONCOLLINEAR = "noncollinear"

# Projection component labels of each layout; the number of blocks itself is read from the file and checked against these.
COMPONENT_LABELS = {
    LAYOUT_NONMAGNETIC: ("total",),
    LAYOUT_COLLINEAR: ("up", "down"),
    LAYOUT_NONCOLLINEAR: ("total", "mx", "my", "mz"),
}

# Band-data arrays: state attribute, name used in messages, and the canonical axis order before the field axis.
# The <dimension> names of the file may come in any order; each block is transposed to this order.
_BAND_ARRAYS = {
    "eigen": ("<eigenvalues>", ("spin", "kpoint", "band")),
    "proj_eigen": ("<projected><eigenvalues>", ("spin", "kpoint", "band")),
    "proj": ("<projected><array>", ("spin", "kpoint", "band", "ion")),
}
# Words for each axis in messages: (position word, block noun).
_AXIS_WORDS = {"spin": ("spin block", "spin"), "kpoint": ("k-point", "k-point"), "band": ("band", "band"), "ion": ("ion", "ion")}
PARAMETERS_OF_INTEREST = ("ISPIN", "LSORBIT", "LNONCOLLINEAR", "LORBIT", "NBANDS")
STRUCTURE_PREFERENCE = ("finalpos", "calculation", "initialpos")
PROGRESS_MIN_BYTES = 50_000_000

# Tags that drive the stream; <r>, <c>, <v> and <i> are read from their parents and never reported one by one.
_STREAM_TAGS = ("generator", "incar", "kpoints", "parameters", "atominfo", "structure", "calculation",
                "eigenvalues", "projected", "dos", "array", "set", "eigenvalues_kpoints_opt", "projected_kpoints_opt")
_KPOINTS_OPT_TAGS = ("eigenvalues_kpoints_opt", "projected_kpoints_opt")
# Array kind of a <projected> array that is freed without being read (parse_vasprun with with_projections=False).
SKIPPED = "skipped"


# Number of eigenvalue spin blocks of a layout: two for collinear ISPIN = 2, otherwise one.
def eigen_spin_blocks(layout: str) -> int:
    return 2 if layout == LAYOUT_COLLINEAR else 1


class VasprunError(Exception):
    """Base class of every error raised while reading a vasprun.xml."""


class VasprunParseError(VasprunError):
    """The file is not well-formed XML, is truncated, or lacks a block the tool needs."""


class VasprunLayoutError(VasprunError):
    """A block has dimensions or a number of entries that disagree with the rest of the file."""


class ProjectionsMissingError(VasprunError):
    """The file has no <projected> block, so band-resolved orbital weights are not available."""


@dataclass
class VasprunData:
    """Everything stage 1 extracts from one vasprun.xml; see the module docstring for the array shapes."""

    path: str
    vasp_version: str | None
    parameters: dict[str, Any]
    incar: dict[str, Any]
    layout: str
    species: list[str]
    lattice: np.ndarray
    reciprocal_lattice: np.ndarray
    positions: np.ndarray
    structure_source: str
    kpoints: np.ndarray
    kpoint_weights: np.ndarray
    kpoint_label_entries: list[tuple[int, str]]
    kpoint_generation: str | None
    kpoint_divisions: int | None
    efermi: float | None
    eigenvalues: np.ndarray
    occupations: np.ndarray
    projections: np.ndarray | None
    projection_components: tuple[str, ...]
    orbital_labels: list[str]
    projected_eigenvalues: np.ndarray | None = None
    has_partial_dos: bool = False
    has_projected_block: bool = False
    notes: list[str] = field(default_factory=list)

    # Number of ions in the structure.
    @property
    def nions(self) -> int:
        return len(self.species)

    # Number of k-points in the file's k-point list.
    @property
    def nkpoints(self) -> int:
        return int(self.kpoints.shape[0])

    # Number of bands per k-point and spin.
    @property
    def nbands(self) -> int:
        return int(self.eigenvalues.shape[2])

    # Convert a 1-based VASP ion number to the 0-based array index, with a bounds check.
    def ion_index(self, vasp_index: int) -> int:
        if isinstance(vasp_index, bool) or not isinstance(vasp_index, (int, np.integer)):
            raise TypeError(f"ion number must be an integer, got {vasp_index!r}")
        if not 1 <= vasp_index <= self.nions:
            raise IndexError(f"ion {vasp_index} is outside 1..{self.nions} (VASP numbering starts at 1)")
        return int(vasp_index) - 1

    # 1-based VASP ion numbers of every ion of one chemical element.
    def ions_of_element(self, element: str) -> list[int]:
        found = [i + 1 for i, name in enumerate(self.species) if name == element]
        if not found:
            raise KeyError(f"no ion of element {element!r}; elements present: {sorted(set(self.species))}")
        return found

    # 0-based column of an orbital label exactly as written in the file's <field> entries.
    def orbital_index(self, label: str) -> int:
        if label not in self.orbital_labels:
            raise KeyError(f"orbital {label!r} not in this file; available: {', '.join(self.orbital_labels)}")
        return self.orbital_labels.index(label)

    # Plain-text summary of the parsed file, used by the "inspect" command.
    def summary_lines(self) -> list[str]:
        lines = [f"file: {self.path}",
                 f"VASP version: {self.vasp_version or 'not stated'}",
                 f"layout: {self.layout} (ISPIN {self.parameters.get('ISPIN')}, LSORBIT {self.parameters.get('LSORBIT')}, "
                 f"LNONCOLLINEAR {self.parameters.get('LNONCOLLINEAR')}, LORBIT {self.parameters.get('LORBIT')})",
                 f"ions: {self.nions} ({', '.join(f'{i + 1}:{s}' for i, s in enumerate(self.species))})",
                 f"structure from: {self.structure_source}",
                 f"k-points: {self.nkpoints}, bands: {self.nbands}, eigenvalue spin blocks: {self.eigenvalues.shape[0]}",
                 f"Fermi energy: {self.efermi if self.efermi is not None else 'not found'} eV"]
        if self.kpoint_generation:
            divisions = f", divisions {self.kpoint_divisions}" if self.kpoint_divisions is not None else ""
            lines.append(f"k-point generation: {self.kpoint_generation}{divisions}")
        if self.kpoint_label_entries:
            lines.append("k-point label entries (corner number:name): "
                         + ", ".join(f"{k}:{name}" for k, name in self.kpoint_label_entries))
        if self.projections is None and self.has_projected_block:
            shape = (len(self.projection_components), self.nkpoints, self.nbands, self.nions, len(self.orbital_labels))
            lines.append(f"projections: present, not read; shape {shape} = (component, kpoint, band, ion, orbital)")
            lines.append(f"components: {', '.join(self.projection_components)}")
            lines.append(f"orbitals: {', '.join(self.orbital_labels)}")
        elif self.projections is None:
            lines.append("projections: none (no <projected> block)")
        else:
            lines.append(f"projections: shape {self.projections.shape} = (component, kpoint, band, ion, orbital)")
            lines.append(f"components: {', '.join(self.projection_components)}")
            lines.append(f"orbitals: {', '.join(self.orbital_labels)}")
            total = self.projections[0]
            lines.append(f"weights of component '{self.projection_components[0]}': min {total.min():.4f}, max {total.max():.4f}")
        lines.extend(f"note: {note}" for note in self.notes)
        return lines


class _ProgressReader:
    """Stream wrapper that prints to stderr, every ten percent, how far the file on disk has been read."""

    # Wrap the stream to parse; disk is the file object on disk (the same object for uncompressed input) and size its length.
    def __init__(self, stream: IO[bytes], disk: IO[bytes], size: int, label: str) -> None:
        self.stream = stream
        self.disk = disk
        self.size = max(size, 1)
        self.label = label
        self.next_step = 10

    # Read a chunk and print progress when the position on disk crosses a new ten-percent step.
    def read(self, n: int = -1) -> bytes:
        chunk = self.stream.read(n)
        percent = 100 * self.disk.tell() // self.size
        while percent >= self.next_step and self.next_step <= 100:
            print(f"reading {self.label}: {self.next_step}%", file=sys.stderr, flush=True)
            self.next_step += 10
        return chunk


# Convert the text of an <i> or <v> parameter to bool, int, float or a list of numbers.
def _convert_value(element: etree._Element) -> Any:
    text = (element.text or "").strip()
    kind = element.get("type", "")
    if kind == "logical":
        return text.upper().startswith("T")
    if kind == "string":
        return text
    if element.tag == "v":
        values = text.split()
        return [int(v) for v in values] if kind == "int" else [float(v) for v in values]
    if kind == "int":
        return int(text)
    try:
        return float(text)
    except ValueError:
        return text


# Read every <i> and <v> name of interest below an element, the last occurrence winning.
def _read_parameters(element: etree._Element, recursive: bool) -> dict[str, Any]:
    found: dict[str, Any] = {}
    children = element.iter("i", "v") if recursive else element.iterchildren("i", "v")
    for child in children:
        name = child.get("name")
        if name:
            found[name] = _convert_value(child)
    return found


# Parse the rows of a <varray> into a float array of shape (rows, columns).
def _varray(element: etree._Element) -> np.ndarray:
    rows = []
    for number, v in enumerate(element.iterchildren("v"), start=1):
        if not (v.text or "").strip():
            raise VasprunParseError(f"<varray name=\"{element.get('name', '')}\"> row {number} is an empty <v></v>")
        rows.append(np.array(v.text.split(), dtype=float))
    return np.array(rows) if rows else np.zeros((0, 0))


# Convert the texts of many <r> rows to one flat float array; a token that is not a number raises ValueError.
def _texts_to_array(texts: list[str]) -> np.ndarray:
    return np.fromstring(" ".join(texts), dtype=float, sep=" ")


# Delete an element's content and its already-handled earlier siblings to keep memory flat.
def _free(element: etree._Element) -> None:
    element.clear()
    parent = element.getparent()
    if parent is not None:
        while element.getprevious() is not None:
            del parent[0]


# True when any ancestor of the element has one of the given tags.
def _inside(element: etree._Element, tags: tuple[str, ...]) -> bool:
    return any(True for _ in element.iterancestors(*tags))


# Decide the layout from the authoritative <parameters> block (not <incar>, which VASP may override).
def detect_layout(parameters: dict[str, Any]) -> str:
    if parameters.get("LNONCOLLINEAR") or parameters.get("LSORBIT"):
        return LAYOUT_NONCOLLINEAR
    if parameters.get("ISPIN", 1) == 2:
        return LAYOUT_COLLINEAR
    return LAYOUT_NONMAGNETIC


# Orbital labels of a <projected> array, refusing an empty <field> so later name matching never meets a blank label.
def _orbital_labels(fields: list[str], name: str) -> list[str]:
    if "" in fields:
        raise VasprunParseError(f"{name} has an empty <field> orbital label (column {fields.index('') + 1})")
    return fields


# Dimension names (dim 1 first) and field labels of an <array>, read from its already parsed header children.
def _array_header(array: etree._Element) -> tuple[list[str], list[str]]:
    dims = sorted(((int(d.get("dim", "0")), (d.text or "").strip()) for d in array.iterchildren("dimension")))
    fields = [(f.text or "").strip() for f in array.iterchildren("field")]
    return [name for _, name in dims], fields


@dataclass
class _ArrayPlan:
    """How the band-data <array> being streamed maps onto its output array."""

    kind: str
    name: str
    order: tuple[str, ...]
    sizes: dict[str, int]
    rows: etree.XPath
    inner_shape: tuple[int, ...]
    perm: tuple[int, ...] | None
    axis_pos: tuple[int, int]
    index: list[int] = field(default_factory=lambda: [-1, -1])


class _StreamState:
    """Mutable state of one pass over the file; kept apart from VasprunData so partial results never leak."""

    # Start with nothing read.
    def __init__(self) -> None:
        self.version: str | None = None
        self.incar: dict[str, Any] = {}
        self.parameters: dict[str, Any] = {}
        self.layout: str | None = None
        self.species: list[str] = []
        self.structures: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        self.kpoints: np.ndarray | None = None
        self.weights: np.ndarray | None = None
        self.kpoint_labels: list[tuple[int, str]] = []
        self.generation: str | None = None
        self.divisions: int | None = None
        self.efermi: float | None = None
        self.eigen: np.ndarray | None = None
        self.proj_eigen: np.ndarray | None = None
        self.proj: np.ndarray | None = None
        self.orbitals: list[str] = []
        self.has_partial_dos = False
        self.notes: list[str] = []
        # Whether <projected> arrays are stored, and whether the current <calculation> has a <projected> block at all.
        self.with_projections = True
        self.saw_projected = False
        # Kind of the band-data array being streamed (None outside one), its plan once the header is read, and the <set> depth.
        self.array_kind: str | None = None
        self.plan: _ArrayPlan | None = None
        self.set_depth = 0
        self.dims: tuple[int, int, int] | None = None
        # Number of <calculation> blocks started, and whether an earlier one (since discarded) had a <projected> block.
        self.calculations = 0
        self.earlier_projections = False

    # Start a new <calculation>: drop the per-calculation data of the earlier one so the result never mixes two.
    def start_calculation(self) -> None:
        self.calculations += 1
        self.earlier_projections = self.earlier_projections or self.saw_projected
        self.saw_projected = False
        self.eigen = self.proj_eigen = self.proj = None
        self.orbitals = []
        self.efermi = None
        self.has_partial_dos = False

    # Number of k-points, bands and ions known before the band data starts; raise when a header block is missing.
    def sizes(self) -> tuple[int, int, int]:
        if self.dims is not None:
            return self.dims
        if self.kpoints is None:
            raise VasprunParseError("no <kpoints><varray name=\"kpointlist\"> before the band data")
        if not self.species:
            raise VasprunParseError("no <atominfo> ion list before the band data")
        nbands = self.parameters.get("NBANDS")
        if not isinstance(nbands, int):
            raise VasprunParseError("NBANDS missing from <parameters>")
        self.dims = (self.kpoints.shape[0], nbands, len(self.species))
        return self.dims


# Kind of band-data array an <array> starts, from its parent and grandparent tags; None for every other array.
def _array_kind(element: etree._Element, with_projections: bool = True) -> str | None:
    parent = element.getparent()
    if parent is None:
        return None
    if parent.tag == "projected":
        return "proj" if with_projections else SKIPPED
    grand = parent.getparent()
    if parent.tag != "eigenvalues" or grand is None:
        return None
    kind = {"projected": "proj_eigen", "calculation": "eigen"}.get(grand.tag)
    return SKIPPED if kind == "proj_eigen" and not with_projections else kind


# Sentence end that states how many blocks of an axis the rest of the file implies.
def _expected_phrase(axis: str, count: int, state: _StreamState) -> str:
    if axis == "spin":
        return f"a {state.layout} run has {count}"
    return {"kpoint": f"{count} k-points are listed", "band": f"NBANDS is {count}",
            "ion": f"<atominfo> lists {count} ions"}[axis]


# Position of the current block for messages, such as " spin block 1, k-point 3".
def _where(plan: _ArrayPlan, levels: int) -> str:
    parts = [f"{_AXIS_WORDS[plan.order[i]][0]} {plan.index[i] + 1}" for i in range(levels)]
    return (" " + ", ".join(parts)) if parts else ""


# Handle the start of an <array> or of a <set> inside a band-data array: plan the array and advance the block counters.
def _on_start(element: etree._Element, state: _StreamState) -> None:
    if element.tag == "calculation":
        state.start_calculation()
        return
    if element.tag == "array":
        state.array_kind = _array_kind(element, state.with_projections)
        state.plan = None
        state.set_depth = 0
        return
    if state.array_kind is None:
        return
    state.set_depth += 1
    if state.set_depth == 1 and state.array_kind == SKIPPED:
        # Without projections only the orbital labels of <projected><array> are kept; its numbers are freed unread.
        if element.getparent().getparent().tag == "projected":
            state.saw_projected = True
            state.orbitals = _orbital_labels(_array_header(element.getparent())[1], "the <projected> array")
    elif state.array_kind == SKIPPED:
        return
    elif state.set_depth == 1:
        state.plan = _plan_array(element.getparent(), state.array_kind, state)
    elif state.set_depth in (2, 3):
        plan = state.plan
        level = state.set_depth - 2
        plan.index[level] += 1
        if level == 0:
            plan.index[1] = -1
        axis = plan.order[level]
        limit = plan.sizes[axis]
        if plan.index[level] >= limit:
            # More blocks than the rest of the file allows: report rather than guess which block is which.
            raise VasprunLayoutError(f"{plan.name}{_where(plan, level)} has more than {limit} {_AXIS_WORDS[axis][1]} blocks; "
                                     f"{_expected_phrase(axis, limit, state)}")


# Check an array header, allocate its output array and return the plan that maps file blocks onto it.
def _plan_array(array: etree._Element, kind: str, state: _StreamState) -> _ArrayPlan:
    if state.layout is None:
        raise VasprunParseError("no <parameters> block before the band data")
    nk, nb, nion = state.sizes()
    name, canonical = _BAND_ARRAYS[kind]
    dims, fields = _array_header(array)
    if sorted(dims) != sorted(canonical):
        raise VasprunLayoutError(f"{name} dimensions {dims}, expected the names {list(canonical)} in any order")
    if kind == "proj":
        if not fields:
            raise VasprunLayoutError(f"{name} has no <field> orbital labels")
        state.orbitals = _orbital_labels(fields, name)
        state.saw_projected = True
        blocks = len(COMPONENT_LABELS[state.layout])
    else:
        if fields[:2] != ["eigene", "occ"]:
            raise VasprunLayoutError(f"{name} fields {fields}, expected eigene and occ first")
        blocks = eigen_spin_blocks(state.layout)
    sizes = {"spin": blocks, "kpoint": nk, "band": nb, "ion": nion}
    order = tuple(reversed(dims))
    inner = order[2:]
    remaining = [axis for axis in canonical if axis not in order[:2]] + ["field"]
    perm = tuple((*inner, "field").index(axis) for axis in remaining)
    setattr(state, kind, np.zeros(tuple(sizes[axis] for axis in canonical) + (len(fields),)))
    return _ArrayPlan(kind=kind, name=name, order=order, sizes=sizes,
                      rows=etree.XPath("set/" * (len(inner) - 1) + "r/text()"),
                      inner_shape=tuple(sizes[axis] for axis in inner) + (len(fields),),
                      perm=None if perm == tuple(range(len(perm))) else perm,
                      axis_pos=(canonical.index(order[0]), canonical.index(order[1])))


# Store one finished block of the second-outermost dimension into the output array, checking its size.
def _store_block(element: etree._Element, state: _StreamState) -> None:
    plan = state.plan
    values = _texts_to_array(plan.rows(element))
    expected = int(np.prod(plan.inner_shape))
    if values.size != expected:
        shape = " x ".join(f"{n} {axis}" for n, axis in zip(plan.inner_shape, (*plan.order[2:], "field")))
        raise VasprunLayoutError(f"{plan.name}{_where(plan, 2)}: {values.size} values, expected {shape}")
    block = values.reshape(plan.inner_shape)
    if plan.perm is not None:
        block = block.transpose(plan.perm)
    target: np.ndarray = getattr(state, plan.kind)
    index: list[Any] = [slice(None)] * (target.ndim - 1)
    index[plan.axis_pos[0]], index[plan.axis_pos[1]] = plan.index
    target[tuple(index)] = block


# Check that a finished set of the outermost or second-outermost dimension held exactly the expected number of blocks.
def _check_count(state: _StreamState, level: int) -> None:
    plan = state.plan
    axis = plan.order[level]
    count, limit = plan.index[level] + 1, plan.sizes[axis]
    if count != limit:
        raise VasprunLayoutError(f"{plan.name}{_where(plan, level)} has {count} {_AXIS_WORDS[axis][1]} blocks, "
                                 f"but {_expected_phrase(axis, limit, state)}")


# Read one finished <structure>: lattice, reciprocal lattice and fractional positions.
def _read_structure(element: etree._Element) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    crystal = element.find("crystal")
    if crystal is None:
        raise VasprunParseError("<structure> without <crystal>")
    arrays = {v.get("name"): _varray(v) for v in crystal.iterchildren("varray")}
    positions = element.find("varray[@name='positions']")
    if "basis" not in arrays or "rec_basis" not in arrays or positions is None:
        raise VasprunParseError("<structure> lacks basis, rec_basis or positions")
    return arrays["basis"], arrays["rec_basis"], _varray(positions)


# Read the k-point list, weights and optional labels from the top-level <kpoints> block.
def _read_kpoints(element: etree._Element, state: _StreamState) -> None:
    for varray in element.iterchildren("varray"):
        if varray.get("name") == "kpointlist":
            state.kpoints = _varray(varray)
        elif varray.get("name") == "weights":
            state.weights = _varray(varray).reshape(-1)
    generation = element.find("generation")
    if generation is not None:
        state.generation = generation.get("param")
        divisions = generation.find("i[@name='divisions']")
        if divisions is not None:
            state.divisions = int(divisions.text)
    labels = element.find("kpoints_labels")
    if labels is not None:
        state.kpoint_labels = [(int(i.text), i.get("name", "")) for i in labels.iterchildren("i")]
    if state.kpoints is None:
        raise VasprunParseError("<kpoints> has no kpointlist")


# Read element symbols per ion from <atominfo>.
def _read_atominfo(element: etree._Element, state: _StreamState) -> None:
    for array in element.iterchildren("array"):
        if array.get("name") == "atoms":
            rows = array.find("set")
            if rows is None:
                raise VasprunParseError("<atominfo><array name=\"atoms\"> has no <set> of ions")
            state.species = [(rc[0].text or "").strip() for rc in rows.iterchildren("rc")]
    declared = element.find("atoms")
    if declared is not None and not (declared.text or "").strip():
        raise VasprunParseError("<atominfo><atoms> is empty")
    if declared is not None and int(declared.text) != len(state.species):
        raise VasprunLayoutError(f"<atominfo> declares {declared.text.strip()} ions but lists {len(state.species)}")


# Handle the end of a streamed element; return True when the element may be freed.
def _on_end(element: etree._Element, state: _StreamState) -> bool:
    tag = element.tag
    if tag == "set" and state.array_kind == SKIPPED:
        state.set_depth -= 1
        return True
    if tag == "set" and state.array_kind is not None:
        depth = state.set_depth
        if depth == 3:
            _store_block(element, state)
        elif depth == 2:
            _check_count(state, 1)
        state.set_depth -= 1
        # Sets below the stored block are read with it and freed with it.
        return depth <= 3
    if tag == "array" and state.array_kind is not None:
        if state.plan is not None and state.array_kind != SKIPPED:
            _check_count(state, 0)
        state.array_kind = state.plan = None
        return True
    parent = element.getparent()
    parent_tag = parent.tag if parent is not None else ""
    if tag in _KPOINTS_OPT_TAGS:
        state.notes.append(f"<{tag}> block present (KPOINTS_OPT run); it is not read in stage 1")
    elif tag == "generator":
        version = element.find("i[@name='version']")
        state.version = (version.text or "").strip() if version is not None else None
    elif tag == "incar":
        state.incar = _read_parameters(element, recursive=False)
    elif tag == "parameters":
        state.parameters = _read_parameters(element, recursive=True)
        state.layout = detect_layout(state.parameters)
    elif tag == "kpoints" and parent_tag == "modeling":
        _read_kpoints(element, state)
    elif tag == "atominfo":
        _read_atominfo(element, state)
    elif tag == "structure":
        name = element.get("name") or ("calculation" if parent_tag == "calculation" else "")
        if name in STRUCTURE_PREFERENCE:
            state.structures[name] = _read_structure(element)
    elif tag == "dos" and parent_tag == "calculation":
        efermi = element.find("i[@name='efermi']")
        if efermi is not None:
            state.efermi = float(efermi.text)
        if element.find("partial") is not None:
            state.has_partial_dos = True
    # <atominfo> reads its own <array> and <set> children when it ends, so those are kept until then.
    return tag not in ("array", "set") or not _inside(element, ("atominfo",))


# Open a .xml or .xml.gz file; return the stream to parse and the file on disk, whose position measures progress.
def _open(path: Path, stack: contextlib.ExitStack) -> tuple[IO[bytes], IO[bytes]]:
    disk = stack.enter_context(path.open("rb"))
    if path.suffix == ".gz":
        return stack.enter_context(gzip.GzipFile(fileobj=disk, mode="rb")), disk
    return disk, disk


# Stream-parse a vasprun.xml (or .xml.gz) and return its structure, k-points, eigenvalues and projections.
def parse_vasprun(path: str | os.PathLike[str], require_projections: bool = True, progress: bool | None = None,
                  with_projections: bool = True) -> VasprunData:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no file at {path}")
    size = path.stat().st_size
    report = size >= PROGRESS_MIN_BYTES if progress is None else progress
    state = _StreamState()
    state.with_projections = with_projections
    with contextlib.ExitStack() as stack:
        stream, disk = _open(path, stack)
        # For .xml.gz input the percentage counts compressed bytes read from disk.
        label = f"{path.name} (compressed bytes)" if path.suffix == ".gz" else path.name
        source = _ProgressReader(stream, disk, size, label) if report else stream
        try:
            for event, element in etree.iterparse(source, events=("start", "end"), tag=_STREAM_TAGS, huge_tree=True):
                if event == "start":
                    _on_start(element, state)
                elif _on_end(element, state):
                    _free(element)
        except etree.XMLSyntaxError as error:
            raise VasprunParseError(f"{path} is not complete, well-formed XML ({error}); the run may have stopped early, "
                                    "or the file was cut while copying") from error
        except (ValueError, IndexError, AttributeError, TypeError) as error:
            raise VasprunParseError(f"{path}: a number or row could not be read ({error})") from error
    return _assemble(path, state, require_projections)


# Check the collected pieces for consistency and build the VasprunData result.
def _assemble(path: Path, state: _StreamState, require_projections: bool) -> VasprunData:
    nk, nb, nion = state.sizes()
    if state.eigen is None or state.layout is None:
        raise VasprunParseError(f"{path} has no <calculation><eigenvalues> block")
    layout = state.layout
    if not state.saw_projected and state.earlier_projections:
        raise VasprunParseError(f"{path} has {state.calculations} <calculation> blocks; an earlier one has a <projected> block "
                                "but the last one has none, so its eigenvalues and those projections would not match")
    if not state.saw_projected and require_projections:
        hint = " The <dos><partial> block is energy-integrated DOS and is not used in its place." if state.has_partial_dos else ""
        raise ProjectionsMissingError(f"{path} has no <projected> block, so there are no band-resolved orbital weights. "
                                      f"Rerun the band calculation with LORBIT = 11 (LORBIT is {state.parameters.get('LORBIT')})." + hint)
    if state.weights is None or state.weights.shape[0] != nk:
        raise VasprunLayoutError("k-point weights missing or not one per k-point")
    source = next((name for name in STRUCTURE_PREFERENCE if name in state.structures), None)
    if source is None:
        raise VasprunParseError(f"{path} has no usable <structure> block")
    lattice, reciprocal, positions = state.structures[source]
    if positions.shape != (nion, 3):
        raise VasprunLayoutError(f"structure has {positions.shape[0]} positions for {nion} ions")
    notes = list(state.notes)
    if state.saw_projected and not state.with_projections:
        notes.append("<projected> block present; its numbers were not read")
    if state.calculations > 1:
        notes.append(f"{state.calculations} calculation blocks; last one used")
    if state.incar.get("ISPIN") is not None and state.incar.get("ISPIN") != state.parameters.get("ISPIN"):
        notes.append(f"INCAR ISPIN = {state.incar['ISPIN']} but VASP used ISPIN = {state.parameters.get('ISPIN')}; "
                     "the <parameters> value is the one applied")
    if state.efermi is None:
        notes.append("no Fermi energy in <dos>; supply one by hand when plotting")
    proj_eigen = state.proj_eigen[..., 0] if state.proj_eigen is not None else None
    if proj_eigen is not None:
        gap = float(np.max(np.abs(proj_eigen - state.eigen[..., 0])))
        if gap > 0:
            notes.append(f"<projected><eigenvalues> differ from <eigenvalues> by up to {gap:.2e} eV")
    components = COMPONENT_LABELS[layout] if state.saw_projected else ()
    return VasprunData(path=str(path), vasp_version=state.version, parameters=state.parameters, incar=state.incar,
                       layout=layout, species=state.species, lattice=lattice, reciprocal_lattice=reciprocal,
                       positions=positions, structure_source=source, kpoints=state.kpoints, kpoint_weights=state.weights,
                       kpoint_label_entries=state.kpoint_labels, kpoint_generation=state.generation,
                       kpoint_divisions=state.divisions, efermi=state.efermi, eigenvalues=state.eigen[..., 0].copy(),
                       occupations=state.eigen[..., 1].copy(), projections=state.proj, projection_components=components,
                       orbital_labels=state.orbitals, projected_eigenvalues=proj_eigen, has_partial_dos=state.has_partial_dos,
                       has_projected_block=state.saw_projected, notes=notes)
