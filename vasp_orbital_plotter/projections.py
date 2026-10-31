"""Atom and orbital selection and the summed projection weight of each state (stage 2 of the specification).

Users name ions in VASP's 1-based numbering; the conversion to 0-based array indices goes through VasprunData.ion_index.
Orbital names are matched against the file's own <field> labels, which are the authority, after one normalising rule.
The weight of a state is the plain sum of the stored values over the selected (ion, orbital) pairs; nothing is squared.
Stage 3 chooses the block: --spin up, down or both for collinear ISPIN = 2 files (each spin with its own eigenvalues),
--component total, mx, my or mz for noncollinear files; mx, my and mz are signed and never treated as spin channels.
Band window, k-point stride and selected ions are cut before the weight sum, so its cost scales with what is drawn.
COMPATIBILITY is the one table of which flags each (layout, channel) pair accepts, with every refusal text.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np

from .options import MAX_SERIES, split_list
from .parser import LAYOUT_COLLINEAR, LAYOUT_NONCOLLINEAR, LAYOUT_NONMAGNETIC, VasprunData


class SelectionError(ValueError):
    """A user selection names an ion, element or orbital that the file does not have, or is written wrongly."""


SPINS = ("up", "down", "both")
COMPONENTS = ("total", "mx", "my", "mz")
MAGNETISATION = ("mx", "my", "mz")

SHELLS = "spdf"
# Number of components a complete s, p, d or f shell has; a shell group is offered only when all are in the file.
SHELL_SIZES = {"s": 1, "p": 3, "d": 5, "f": 7}
# The only hand-written groups: cubic-field subsets of the d shell, in VASP's spelling of the file labels.
HAND_GROUPS: dict[str, tuple[str, ...]] = {"t2g": ("dxy", "dxz", "dyz"), "eg": ("dz2", "x2-y2")}


# Normalised spelling of an orbital name: lower case without spaces, "_", "^" or "-", so "d_{x^2-y^2}" and "dx2y2" agree.
def normalise(name: str) -> str:
    return re.sub(r"[\s_^{}\-]+", "", str(name)).lower()


# Name without its leading shell letter (s, p, d, f), or None when there is none or nothing would remain.
def _without_shell(key: str) -> str | None:
    return key[1:] if len(key) > 1 and key[0] in SHELLS else None


# 0-based column of one orbital name among the file's labels, or None; a leading shell letter is optional on either side.
def orbital_column(name: str, labels: Sequence[str]) -> int | None:
    keys = [normalise(label) for label in labels]
    key = normalise(name)
    for candidates, wanted in ((keys, key), (keys, _without_shell(key)), ([_without_shell(k) for k in keys], key)):
        if wanted is not None and wanted in candidates:
            return candidates.index(wanted)
    return None


# Columns of each shell in the file: a label with a shell letter starts or continues that shell; one without continues the last.
def shell_columns(labels: Sequence[str]) -> dict[str, list[int]]:
    shells: dict[str, list[int]] = {}
    current = None
    for column, label in enumerate(labels):
        key = normalise(label)
        current = key[0] if key[:1] in SHELL_SIZES else current
        if current is not None:
            shells.setdefault(current, []).append(column)
    return shells


# Every group name with its columns in the file, or with the missing members when the group is incomplete.
def _groups(labels: Sequence[str]) -> dict[str, tuple[list[int], list[str]]]:
    groups: dict[str, tuple[list[int], list[str]]] = {}
    shells = shell_columns(labels)
    for shell, size in SHELL_SIZES.items():
        columns = shells.get(shell, [])
        missing = [] if len(columns) == size else [f"{size - len(columns)} of the {size} {shell} components"]
        groups[shell] = (columns, missing)
    for group, members in HAND_GROUPS.items():
        found = [(member, orbital_column(member, labels)) for member in members]
        groups[group] = ([c for _, c in found if c is not None], [m for m, c in found if c is None])
    return groups


# Groups whose every member is among the file's orbital labels: s, p, d, f first, then t2g and eg.
def available_groups(labels: Sequence[str]) -> list[str]:
    return [group for group, (_, missing) in _groups(labels).items() if not missing]


# Text listing the file's orbital labels and usable groups, for error messages.
def _available(labels: Sequence[str]) -> str:
    return f"available orbitals: {', '.join(labels)}; groups: {', '.join(available_groups(labels) + ['all'])}"


# 0-based orbital columns of one name or group; a name that is not in the file raises with the available labels.
def resolve_orbital(name: str, labels: Sequence[str]) -> list[int]:
    key = normalise(name)
    if key == "all":
        return list(range(len(labels)))
    groups = _groups(labels)
    if key in groups and (key not in SHELL_SIZES or orbital_column(name, labels) is None):
        columns, missing = groups[key]
        if missing:
            raise SelectionError(f"orbital group {name!r} is incomplete in this file (missing: {', '.join(missing)}); "
                                 f"{_available(labels)}")
        return columns
    column = orbital_column(name, labels)
    if column is None:
        raise SelectionError(f"orbital {name!r} is not in this file; {_available(labels)}")
    return [column]


# 0-based orbital columns of several names or groups, without duplicates, in first-named order.
def resolve_orbitals(names: Iterable[str], labels: Sequence[str]) -> list[int]:
    columns: list[int] = []
    for name in names:
        columns += [column for column in resolve_orbital(name, labels) if column not in columns]
    return columns


# 1-based VASP ion numbers named by numbers, ranges ("3-5") or element symbols, checked against the file.
def parse_atoms(items: Iterable[str], data: VasprunData) -> list[int]:
    ions: list[int] = []
    for item in items:
        match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", item)
        if match:
            first, last = int(match.group(1)), int(match.group(2))
            if first > last:
                raise SelectionError(f"atom range {item!r} runs backwards; write it as {last}-{first}")
            numbers = list(range(first, last + 1))
        elif item.isdigit():
            numbers = [int(item)]
        else:
            numbers = ions_of_element(item, data)
        for number in numbers:
            try:
                data.ion_index(number)
            except IndexError as error:
                raise SelectionError(f"atom {number}: {error}") from None
            if number not in ions:
                ions.append(number)
    return ions


# 1-based numbers of every ion of an element, with a SelectionError that lists the elements present.
def ions_of_element(element: str, data: VasprunData) -> list[int]:
    try:
        return data.ions_of_element(element)
    except KeyError:
        raise SelectionError(f"no ion of element {element!r} in this file; elements present: "
                             f"{', '.join(dict.fromkeys(data.species))}") from None


@dataclass
class Selection:
    """A set of (ion, orbital) pairs: ions as 1-based VASP numbers, orbitals as 0-based columns of the file."""

    pairs: list[tuple[int, int]]
    label: str


# Short text for a set of ions: an element name when it covers exactly that element's ions, else the numbers.
def _ions_text(ions: Sequence[int], data: VasprunData) -> str:
    elements = list(dict.fromkeys(data.species[n - 1] for n in ions))
    if len(elements) == 1 and sorted(ions) == data.ions_of_element(elements[0]):
        return elements[0]
    if sorted(ions) == list(range(1, data.nions + 1)):
        return "all atoms"
    return " ".join(f"{data.species[n - 1]}{n}" for n in ions)


# Build one selection term from atom items and orbital names; empty atoms mean all ions, empty orbitals all orbitals.
def _term(atom_items: Sequence[str], orbital_names: Sequence[str], data: VasprunData) -> tuple[list[tuple[int, int]], str]:
    ions = parse_atoms(atom_items, data) if atom_items else list(range(1, data.nions + 1))
    columns = resolve_orbitals(orbital_names, data.orbital_labels) if orbital_names else list(range(len(data.orbital_labels)))
    orbitals_text = "+".join(orbital_names) if orbital_names else "all orbitals"
    return [(ion, col) for ion in ions for col in columns], f"{_ions_text(ions, data)} {orbitals_text}"


# Build the selection from --atoms, --element, --orbitals and repeated --select "ATOMS:ORBITALS" terms.
def build_selection(data: VasprunData, atoms: str | None = None, elements: str | None = None, orbitals: str | None = None,
                    select: Sequence[str] | None = None) -> Selection:
    if not data.orbital_labels:
        raise SelectionError("this file has no orbital labels, so nothing can be selected")
    terms: list[tuple[list[tuple[int, int]], str]] = []
    plain_atoms = split_list(atoms) + split_list(elements)
    if select:
        if plain_atoms or orbitals:
            raise SelectionError("use either --select terms or --atoms/--element/--orbitals, not both")
        for text in select:
            if text.count(":") != 1:
                raise SelectionError(f"--select {text!r} must have the form ATOMS:ORBITALS, for example Mn:d or 1,3:pz")
            atom_part, orbital_part = text.split(":")
            terms.append(_term(split_list(atom_part), split_list(orbital_part), data))
    else:
        terms.append(_term(plain_atoms, split_list(orbitals), data))
    pairs = list(dict.fromkeys(pair for term_pairs, _ in terms for pair in term_pairs))
    if not pairs:
        raise SelectionError("the selection is empty")
    return Selection(pairs=pairs, label=" + ".join(label for _, label in terms))


_SPIN_NONMAGNETIC = "--spin needs a collinear ISPIN = 2 file; this file is nonmagnetic and has one block"
_SPIN_NONCOLLINEAR = ("--spin needs a collinear ISPIN = 2 file; this file is noncollinear, whose blocks are total, mx, my "
                      "and mz, not spin channels: use --component total, mx, my or mz")


@dataclass(frozen=True)
class ChannelRule:
    """One (layout, channel) entry: projection blocks drawn (one panel each), signed or not, refused flags with their text.

    An entry with `error` set is a channel the layout does not have; --spin or --component asking for it gets that text.
    """

    blocks: tuple[str, ...] = ()
    signed: bool = False
    refuses: dict[str, str] = field(default_factory=dict)
    error: str | None = None


# Error texts of a signed magnetisation channel for --normalize and --overlay.
def _signed_refusals(name: str) -> dict[str, str]:
    return {"--normalize": f"--normalize divides by the total of each state and is not defined for the signed {name} block; "
                           "drop --normalize",
            "--overlay": f"--overlay is not available for the signed {name} block, whose two sign colours would clash "
                         "with the series colours; plot one selection at a time"}


# The one compatibility table of layouts, channels (--spin or --component values) and flags; every refusal text lives here.
COMPATIBILITY: dict[tuple[str, str], ChannelRule] = {
    (LAYOUT_NONMAGNETIC, "total"): ChannelRule(("total",)),
    **{(LAYOUT_NONMAGNETIC, name): ChannelRule(error=f"--component {name} needs a noncollinear file; this file is "
                                                     "nonmagnetic, whose only block is total") for name in MAGNETISATION},
    **{(LAYOUT_NONMAGNETIC, name): ChannelRule(error=_SPIN_NONMAGNETIC) for name in SPINS},
    # --component total is accepted on collinear files and changes nothing; the panels follow --spin.
    (LAYOUT_COLLINEAR, "total"): ChannelRule(("up", "down")),
    (LAYOUT_COLLINEAR, "both"): ChannelRule(("up", "down")),
    (LAYOUT_COLLINEAR, "up"): ChannelRule(("up",)),
    (LAYOUT_COLLINEAR, "down"): ChannelRule(("down",)),
    **{(LAYOUT_COLLINEAR, name): ChannelRule(error=f"--component {name} needs a noncollinear file; this file is collinear "
                                                   "ISPIN = 2, whose blocks are spin up and spin down: use --spin up, down "
                                                   "or both") for name in MAGNETISATION},
    (LAYOUT_NONCOLLINEAR, "total"): ChannelRule(("total",)),
    **{(LAYOUT_NONCOLLINEAR, name): ChannelRule((name,), signed=True, refuses=_signed_refusals(name))
       for name in MAGNETISATION},
    **{(LAYOUT_NONCOLLINEAR, name): ChannelRule(error=_SPIN_NONCOLLINEAR) for name in SPINS},
}
# Channel drawn when neither --spin nor --component is given.
DEFAULT_CHANNEL = {LAYOUT_NONMAGNETIC: "total", LAYOUT_COLLINEAR: "both", LAYOUT_NONCOLLINEAR: "total"}


# Refusal text of a flag with a channel, or None; without a layout, every layout that can draw the channel is consulted.
def flag_refusal(channel: str, flag: str, layout: str | None = None) -> str | None:
    for (layout_name, name), rule in COMPATIBILITY.items():
        if name == channel and rule.error is None and layout in (None, layout_name) and flag in rule.refuses:
            return rule.refuses[flag]
    return None


@dataclass
class Block:
    """Source of one panel: the projection block index and the panel heading (None for one panel)."""

    projection: int
    heading: str | None


@dataclass
class BlockChoice:
    """The panels to draw and how to name them: label for text, tag for output names, quantity for the axis and legend."""

    blocks: list[Block]
    label: str
    tag: str
    quantity: str
    signed: bool
    collinear: bool
    caption: str | None
    note: str | None

    # Eigenvalue block of a panel: each spin has its own for collinear files, otherwise there is only block 0.
    def eigen(self, block: Block) -> int:
        return block.projection if self.collinear else 0


# Blocks for --spin and --component, checked against COMPATIBILITY together with --normalize and --overlay.
def choose_blocks(data: VasprunData, spin: str | None = None, component: str | None = None, normalize: bool = False,
                  overlay: bool = False) -> BlockChoice:
    if data.projections is None:
        raise SelectionError("this file has no <projected> block; rerun the band calculation with LORBIT = 11")
    layout = data.layout
    for flag, value, names in (("--spin", spin, SPINS), ("--component", component, COMPONENTS)):
        if value is None:
            continue
        rule = COMPATIBILITY.get((layout, value))
        if rule is None:
            raise SelectionError(f"{flag} {value!r} must be one of {', '.join(names)}")
        if rule.error:
            raise SelectionError(rule.error)
    collinear = layout == LAYOUT_COLLINEAR
    channel = (spin if collinear else component) or DEFAULT_CHANNEL[layout]
    rule = COMPATIBILITY[(layout, channel)]
    for flag, given in (("--normalize", normalize), ("--overlay", overlay)):
        if given and flag in rule.refuses:
            raise SelectionError(rule.refuses[flag])
    labels = list(data.projection_components)
    several = len(rule.blocks) > 1
    blocks = [Block(labels.index(name), f"spin {name}" if several else None) for name in rule.blocks]
    if collinear:
        note = None if spin else ("collinear ISPIN = 2 file: spin up and spin down are drawn as two panels; "
                                  "choose one with --spin up or --spin down")
        label = "spin up and down" if several else f"spin {channel}"
        return BlockChoice(blocks, label, channel, "projection weight", False, True, None, note)
    caption = None
    if layout == LAYOUT_NONCOLLINEAR:
        value = data.parameters.get("LSORBIT")
        lsorbit = "not stated" if value is None else ("T" if value else "F")
        if rule.signed:
            caption = f"{channel} block as read from vasprun.xml, sign convention unverified; LSORBIT = {lsorbit}"
        elif lsorbit != "T":
            caption = f"noncollinear run, LSORBIT = {lsorbit}: total charge block, plotted like a spin-orbit run"
    return BlockChoice(blocks, channel, channel if rule.signed else "", channel if rule.signed else "projection weight",
                       rule.signed, False, caption, None)


# Contiguous band slice of every band that enters [emin, emax] at some k-point in any of the (kpoint, band) arrays.
def bands_reaching(energies: Sequence[np.ndarray], emin: float, emax: float) -> slice:
    reach = np.zeros(energies[0].shape[1], dtype=bool)
    for es in energies:
        reach |= (es.max(axis=0) >= emin) & (es.min(axis=0) <= emax)
    found = np.flatnonzero(reach)
    return slice(int(found[0]), int(found[-1]) + 1) if found.size else slice(0, 0)


# One selection per --select term for --overlay, labelled by the term or by --label-series in the same order.
def build_series(data: VasprunData, select: Sequence[str] | None, labels: Sequence[str] | None = None) -> list[Selection]:
    if not select:
        raise SelectionError("--overlay needs --select ATOMS:ORBITALS terms, one per series, for example "
                             "--select Mn:d --select Te:p")
    if len(select) > MAX_SERIES:
        raise SelectionError(f"--overlay draws at most {MAX_SERIES} series, but {len(select)} --select terms were given")
    if labels and len(labels) != len(select):
        raise SelectionError(f"{len(labels)} --label-series names for {len(select)} --select terms; give one name per term")
    series = [build_selection(data, select=[term]) for term in select]
    for selection, label in zip(series, labels or []):
        selection.label = label
    return series


# Increasing indices as basic slices of constant step, in order; a stride of each segment plus its end gives a few slices.
def _runs(indices: np.ndarray) -> list[slice]:
    runs: list[slice] = []
    i, n = 0, len(indices)
    while i < n:
        if i + 1 == n:
            runs.append(slice(int(indices[i]), int(indices[i]) + 1))
            break
        j, step = i + 1, int(indices[i + 1] - indices[i])
        while j + 1 < n and indices[j + 1] - indices[j] == step:
            j += 1
        runs.append(slice(int(indices[i]), int(indices[j]) + 1, step))
        i = j + 1
    return runs


# Views of block `component` at the kept k-points (one per run of _runs) inside the band window and the ion range.
def _views(data: VasprunData, component: int, kpoints: np.ndarray | None, bands: slice | None,
           ions: slice = slice(None)) -> list[np.ndarray]:
    """Views, never copies, so the k-point stride, the band window and the ion range all cut the work of the sums."""
    block = data.projections[component]
    rows = [slice(None)] if kpoints is None else _runs(np.asarray(kpoints))
    return [block[r, slice(None) if bands is None else bands, ions] for r in rows]


# Total of every state over all ions and orbitals, shape (kpoint, band), for --normalize; computed once per panel.
def state_totals(data: VasprunData, component: int, kpoints: np.ndarray | None = None,
                 bands: slice | None = None) -> np.ndarray:
    return np.concatenate([view.sum(axis=(2, 3)) for view in _views(data, component, kpoints, bands)])


# Summed weight of every state, shape (kpoint, band), from projection block `component`, over the given k-points and bands.
def select_weights(data: VasprunData, selection: Selection, component: int, normalize: bool = False,
                   kpoints: np.ndarray | None = None, bands: slice | None = None,
                   total: np.ndarray | None = None) -> np.ndarray:
    """The einsum runs on views cut to the band window, the kept k-points and the range of selected ions; `total` is reused if given."""
    indices = [data.ion_index(ion) for ion, _ in selection.pairs]
    first = min(indices)
    mask = np.zeros((max(indices) - first + 1, data.projections.shape[-1]), dtype=data.projections.dtype)
    for index, (_, orbital) in zip(indices, selection.pairs):
        mask[index - first, orbital] = 1.0
    views = _views(data, component, kpoints, bands, slice(first, first + mask.shape[0]))
    weights = np.concatenate([np.einsum("kbio,io->kb", view, mask) for view in views])
    if normalize:
        total = state_totals(data, component, kpoints, bands) if total is None else total
        weights = np.divide(weights, total, out=np.zeros_like(weights), where=np.abs(total) > 1e-12)
    return weights


# Lines listing every ion with its 1-based VASP number, element and fractional position.
def atom_lines(data: VasprunData) -> list[str]:
    lines = ["number  element  fractional position"]
    for i, (name, pos) in enumerate(zip(data.species, data.positions), start=1):
        lines.append(f"{i:6d}  {name:7s}  {pos[0]:10.6f} {pos[1]:10.6f} {pos[2]:10.6f}")
    return lines


# Lines listing the file's orbital labels by 1-based column, then the usable groups with their members.
def orbital_lines(data: VasprunData) -> list[str]:
    labels = data.orbital_labels
    lines = ["column  label in file"] + [f"{i:6d}  {label}" for i, label in enumerate(labels, start=1)]
    groups = _groups(labels)
    lines.append("groups: " + ", ".join(f"{g} = {'+'.join(labels[c] for c in groups[g][0])}" for g in available_groups(labels))
                 + ", all = every column")
    lines.append("names match the labels above without case, spaces, '_', '^' or '-', and with or without the shell letter")
    return lines
