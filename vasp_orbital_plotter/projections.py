"""Atom and orbital selection and the summed projection weight of each state (stage 2 of the specification).

Users name ions in VASP's 1-based numbering; the conversion to 0-based array indices goes through VasprunData.ion_index.
Orbital names are matched against the file's own <field> labels, which are the authority, after one normalising rule.
The weight of a state is the plain sum of the stored values over the selected (ion, orbital) pairs; nothing is squared.
Stage 2 uses one projection block: "total" for nonmagnetic and noncollinear files, spin up for collinear ISPIN = 2 files.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from .options import split_list
from .parser import LAYOUT_COLLINEAR, LAYOUT_NONCOLLINEAR, VasprunData


class SelectionError(ValueError):
    """A user selection names an ion, element or orbital that the file does not have, or is written wrongly."""


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


# Projection block used in stage 2: its index, its label and a note for the user (None when nothing needs saying).
def select_component(data: VasprunData) -> tuple[int, str, str | None]:
    if data.projections is None:
        raise SelectionError("this file has no <projected> block; rerun the band calculation with LORBIT = 11")
    components = list(data.projection_components)
    if data.layout == LAYOUT_NONCOLLINEAR:
        return components.index("total"), "total", None
    if data.layout == LAYOUT_COLLINEAR:
        return 0, components[0], ("collinear ISPIN = 2 file: stage 2 plots the spin-up bands and weights only; "
                                  "spin choice and two panels come in stage 3")
    return 0, components[0], None


# Summed weight of every state, shape (kpoint, band), from projection block `component`; optionally normalised per state.
def select_weights(data: VasprunData, selection: Selection, component: int, normalize: bool = False) -> np.ndarray:
    block = data.projections[component]
    mask = np.zeros(block.shape[2:])
    for ion, orbital in selection.pairs:
        mask[data.ion_index(ion), orbital] = 1.0
    weights = np.tensordot(block, mask, axes=([2, 3], [0, 1]))
    if normalize:
        total = block.sum(axis=(2, 3))
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
