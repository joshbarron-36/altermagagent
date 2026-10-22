Written by Diana Golovanova for her comparison module

# Compare a Wannier Hamiltonian with VASP bands

The comparison command computes bands from `seedname_hr.dat` and overlays them on a verified reference k-path.
It uses the indexed calculation lineage to locate the reference.
It also locates existing band figures and records their paths for the researcher.
Existing figures and calculation folders are never changed.

## Run

Use the existing project environment; no new dependency is required.

```bash
altermag import /path/to/calculations --name CrSe
altermag compare-wannier CrSe path/to/wannier90_hr.dat --emin -2 --emax 2
```

The Hamiltonian path can be absolute or relative to the indexed project root.
You can also select its indexed run:

```bash
altermag compare-wannier CrSe --run-name coplanar/2-SOC/3.1-wannier
```

The chat tool is named `compare_wannier`.
Ask: "Compare the SOC Wannier Hamiltonian in coplanar/2-SOC/3.1-wannier with its ancestor bands."
Both chat providers receive this tool definition.

## Reference discovery

The nearest recorded band ancestor is preferred.
Band and Wannier runs often both descend directly from the same self-consistent field (SCF) calculation.
If there is no direct band ancestor, an earlier band run sharing a recorded SCF ancestor is eligible.
The reference must have the same structure, spin–orbit coupling (SOC), spin axis and input magnetic texture.
Multiple eligible references require an explicit choice:

```bash
altermag compare-wannier CrSe --run-name coplanar/2-SOC/3.1-wannier \
  --reference-run coplanar/2-SOC/2-bands
```

Unknown parents, broken links and cycles are reported rather than resolved by folder names or dates.
Use `altermag parent` to correct a parent only when the provenance is known.
The audit retains the original parent notes, including inferred or user-selected links.
An indexed parent link is evidence, not independent proof that a file was copied from that parent.

## Required inputs and checks

The owning Wannier run and reference band run need `POSCAR` and `INCAR`.
Structures must use the same ordered species, lattice basis and atomic coordinates within 0.000001 angstrom.
Equivalent cells requiring a rotation or basis transformation are rejected because their fractional k-paths need transformation too.

Strict mode requires byte-identical `POSCAR`, `INCAR`, `CHGCAR` and `CHG` files.
It records SHA-256 cryptographic hashes and reports parsed `INCAR` differences.
Charge files are read only in bounded chunks for the equality check explicitly requested for this feature.
They are never parsed as densities or modified; the ordinary importer still records only their metadata.
Matching inputs do not prove magnetic or electronic convergence.

Band and Wannier `INCAR` files often differ because the two tasks need different controls.
An explicit preliminary mode allows such differences and missing charge files:

```bash
altermag compare-wannier CrSe --run-name coplanar/2-SOC/3.1-wannier \
  --structure-only --emin -2 --emax 2
```

This mode still requires the same structure and magnetic configuration.
The figure and audit identify the overlay as preliminary whenever required inputs fail equality.
Inspect all differences, especially physical controls such as constraints or interaction parameters.
The chat must not enable this mode without an explicit preliminary-comparison request.

The reference needs reciprocal-coordinate line-mode `KPOINTS` and one of these numerical outputs:

- `EIGENVAL`, preferred to avoid parsing large projected XML files.
- `vasprun.xml`, when `EIGENVAL` is absent.
- `REFORMATTED_BAND.dat` or `REFORMATED_BAND.dat`, plus `FERMI_ENERGY`.

Output k-point coordinates must match the expanded `KPOINTS` coordinates within 0.0000001 in fractional units.
For VASPKIT tables, the tool checks row counts and normalized path distances instead.
Such tables do not independently verify output k-coordinates, and the audit states this limitation.
Spin-polarised VASPKIT tables require `EIGENVAL` or XML to identify channels.
A saved image alone is insufficient for numerical comparison.

## Hamiltonian interpolation

The command reads the real-space Hamiltonian and its degeneracies directly.
It diagonalizes the Fourier transform at the reference k-points, including repeated segment endpoints.
Disconnected segments are drawn separately.

For Wigner–Seitz distance interpolation, provide the corresponding `seedname_wsvec.dat` alongside the Hamiltonian.
Pair-dependent translations are folded once into effective real-space matrices before evaluating the path.
This avoids storing millions of translation objects or repeating their processing at each k-point.
Exports declaring `write_ndegen_applied=true` already contain these weights and shifts; they are not applied twice.
Without a translation file, an accompanying `.win` or export header must explicitly declare `use_ws_distance=false`.
The default distance convention is not silently replaced by another interpolation method.

The `.win` file must belong to this Hamiltonian; `num_wann` is cross-checked when present.
The command rejects non-finite data and Hamiltonians with a Hermiticity residual above 0.00001 electronvolt (eV).
This numerical check does not establish physical accuracy.

Format and convention references:

- [Wannier90 output formats](https://wannier90.readthedocs.io/en/stable/user_guide/wannier90/files/).
- [Wannier90 interpolation notes](https://wannier90.readthedocs.io/en/stable/user_guide/wannier90/notes_interpolations/).

## Energy reference and errors

Both datasets use one absolute Fermi energy.
The command prefers a shared SCF ancestor's `OUTCAR` or `FERMI_ENERGY`.
If unavailable, it uses the band calculation's recorded Fermi energy and identifies that fallback.
You may provide a known common value with `--fermi-eV`.
VASPKIT's recorded Fermi shift is undone before the common value is subtracted.
No fitted energy shift is applied.

For a reference containing two spin channels, select `--spin up` or `--spin down` for the corresponding Hamiltonian.
SOC spinor eigenvalues are compared as a single channel.

Error statistics are optional because disentangled Wannier bands need not map to one contiguous VASP band range.
Only provide `--band-start N` if Wannier band 1 corresponds to VASP band N and all following bands correspond in order.
The command does not infer or optimize that mapping.
For that mapping it reports mean absolute error, root mean square error, the 95th percentile and maximum error in millielectronvolts.
It also reports the sample count and count above 20 millielectronvolts within the requested energy window.
Visual agreement alone does not validate crossings, topology or three-dimensional interpolation accuracy.

## Outputs

Every run creates a unique output directory under `~/.altermagagent/figures/<project>/`.
Use `--out-dir` for a different directory outside all indexed calculation roots.

- `overlay.png`, `overlay.pdf` and `overlay.svg`: the comparison figure.
- `overlay.csv`: fractional k-points, path distances and both sets of energies relative to the common Fermi energy.
- `overlay.json`: sources, original figures, lineage, hashes, differences, interpolation convention and optional errors.

For limited CPU use on a shared machine, set `OPENBLAS_NUM_THREADS=2` before running the command.
