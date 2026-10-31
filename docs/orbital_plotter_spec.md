# Orbital-projected band plotter: specification

This is the specification written by Diana on 26 October 2026, her text, lightly reformatted: one long sentence split, one table rebuilt, five glued line breaks restored.
The duplicated sections 9 to 12 are removed.

```text
Development task: VASP orbital-projected band structure plotting tool
1. Objective
This tool can be run after a band structure calculation. Develop a standalone Python tool that reads vasprun.xml directly and generates publication-quality orbital-projected electronic band structures from VASP calculations.
The tool should provide functionality similar to VASPKIT's orbital-projected band plots (fat-band plots), but without requiring VASPKIT, PROCAR, EIGENVAL, or other VASP output files.
The user should be able to select:

A specific atom by its index in the VASP structure.
Multiple atoms or all atoms of a chemical element.
A specific orbital, such as s, px, py, pz, dxy, dxz, dyz, dz2, or dx2-y2.
Orbital groups such as p, d, t2g, or eg.
Multiple atom–orbital combinations.
The resulting plots should display electronic band energies, with orbital contributions represented through marker sizes or colors.
The tool is intended for actual density functional theory (DFT) research, including magnetic and spin–orbit-coupled materials.
Priorities: numerical correctness, robust parsing, computational efficiency, and publication-quality visualization.
2. Scientific definition
For each electronic state (|\psi_{n\mathbf{k}}\rangle), the orbital-projection weight is obtained from the site- and orbital-resolved information provided by VASP.
For a selected atom (a) and orbital (\alpha), define
[
w_{n\mathbf{k}}^{a,\alpha}
]
as the corresponding projection weight stored in vasprun.xml.
For multiple selected atoms and orbitals, the total weight is
\sum_{a\in A}
\sum_{\alpha\in O}
w_{n\mathbf{k}}^{a,\alpha}.
]
Here:

(A) is the set of selected atomic indices.
(O) is the set of selected orbitals.
(n) is the band index.
(\mathbf{k}) is the reciprocal-space point.
Important:

Use VASP's stored projection weights directly. Do not square them again.
Do not automatically normalize projection weights to unity.
Do not assume that summing all local orbital projections gives exactly one.
Allow optional normalization only as an explicitly selected visualization mode.
Preserve the distinction between orbital-projected charge and orbital-resolved spin magnetization.
3. XML parsing
Implement a robust parser using Python, preferably lxml.etree.iterparse for memory efficiency.
Extract the following information directly from vasprun.xml:
```

| Information | Required |
| --- | --- |
| Number of atoms and their chemical species | Yes |
| Atom ordering and indices | Yes |
| Reciprocal lattice vectors | Yes |
| K-point coordinates and ordering | Yes |
| Band eigenvalues | Yes |
| Fermi energy | Yes |
| Orbital-resolved band weights | Yes |
| Spin configuration | Yes |
| Available orbital-channel labels | Yes |
| High-symmetry point labels | If available |

```text
Investigate the actual XML structure, including:

<atominfo>
<kpoints>
<eigenvalues>
<projected>
Relevant INCAR parameters
The final electronic calculation block
The parser must identify orbital columns from the XML <field> labels instead of assuming that column positions are universal.
Check carefully how different VASP versions represent spin-resolved orbital projections.
Distinguish the band-projection data in <projected> from the energy-integrated information in <dos><partial>.
If the required projections are missing, raise an informative error suggesting an appropriate VASP calculation, commonly using LORBIT = 11.
Do not silently substitute projected density-of-states data for band projections.
Efficiency requirements
Real VASP XML files may be several gigabytes.

Avoid loading the entire XML tree into memory.
Consider a streaming parser or a two-pass extraction strategy.
Minimize unnecessary copies of large arrays.
Support .xml and optionally .xml.gz.
Display clear progress information for large files.
Handle incomplete or corrupted XML with informative error messages.
Keep the parser separate from the plotting functionality.
4. Atom and orbital selection
Provide an intuitive command-line interface.
For example:
python plot_orbitals.py \
    --vasprun vasprun.xml \
    --atoms 1,3 \
    --orbitals dxy,dxz \
    --emin -2 \
    --emax 2 \
    --output bands.png

Support selection by chemical element:
python plot_orbitals.py \
    --vasprun vasprun.xml \
    --element Cr \
    --orbitals d

Provide inspection commands:
python plot_orbitals.py --vasprun vasprun.xml --list-atoms

python plot_orbitals.py --vasprun vasprun.xml --list-orbitalsAtom indices exposed to users should follow VASP's conventional 1-based indexing.
Internally, Python arrays may use 0-based indexing, but the conversion must be explicit and carefully tested.
The program should support orbital groups:

s
p = px + py + pz
d = dxy + dyz + dz2 + dxz + dx2-y2
t2g = dxy + dxz + dyz
eg = dz2 + dx2-y2
Only support a group if all the necessary orbital components are present in the XML.
Use the actual orbital names in the XML as the authoritative source and implement aliases where needed.
If an orbital is unavailable, print the available orbitals and exit gracefully.
5. Band-structure plotting
Use Matplotlib for the initial implementation.
Implement two main visualization modes.
Mode A: Fat-band plot
Plot all electronic bands as thin gray lines.
Overlay markers at each ((k,n)), with marker area proportional to the selected orbital weight:
[
S_{n\mathbf{k}}=S_{\max}f(W_{n\mathbf{k}}).
]
By default, use a linear size mapping.
Allow the user to adjust marker scaling and transparency.
Mode B: Colored orbital projection
Plot the band energies with colors representing orbital weights.
Include:

An appropriate colormap.
A colorbar indicating projection weight.
Adjustable color limits.
A neutral background showing the complete band structure.
Prefer a representation that preserves band connectivity and does not create misleading interpolation between unrelated electronic states.
General plot options
Provide:

Configurable figure dimensions.
Adjustable energy range.
Fermi-energy subtraction, (E-E_F).
Optional manually specified Fermi energy.
Horizontal reference line at zero energy.
High-symmetry vertical lines.
Configurable high-symmetry labels.
Linewidth, marker size, transparency, and font settings.
Optional legend.
Support for PNG, PDF, and SVG output.
Resolution control for raster images.
A clean, publication-ready default style.
Example:
python plot_orbitals.py \
    --vasprun vasprun.xml \
    --atom 2 \
    --orbital pz \
    --mode fatband \
    --emin -3 \
    --emax 3 \
    --marker-scale 80 \
    --output atom2_pz.pdf6. K-path reconstruction
Reconstruct the horizontal axis using the k-point coordinates and reciprocal lattice vectors extracted from vasprun.xml.
For consecutive k-points, calculate
\left|
\mathbf{B}
(\mathbf{k}_{i+1}-\mathbf{k}_i)
\right|,
]
where (\mathbf{B}) transforms fractional reciprocal coordinates into Cartesian reciprocal-space coordinates.
Then
[
x_i=\sum_{j<i}\Delta k_j.
]
Ensure that the reciprocal-space units and lattice convention are handled correctly.
Important requirements:

Preserve the original ordering of the k-points.
Identify repeated endpoints and possible discontinuities between path segments.
Do not draw artificial lines connecting disconnected branches.
Avoid treating unrelated k-points as neighboring points on a continuous path.
High-symmetry labels may not be recoverable from vasprun.xml alone.
Therefore, implement an optional interface for supplying labels manually, for example:
--high-symmetry "1:G,41:M,81:K,121:G"

Interpret these as 1-based k-point indices.
Optionally support reading a VASP KPOINTS file for high-symmetry labels, but do not make this file mandatory.
7. Magnetic and spin–orbit-coupled calculations
This is important because the tool will be used for magnetic and topological materials.
Support three categories of VASP calculations.
A. Nonmagnetic calculations
Read and plot the available charge projection weights.
B. Collinear spin-polarized calculations
Support separate spin-up and spin-down channels.
Example:
python plot_orbitals.py \
    --vasprun vasprun.xml \
    --element Cr \
    --orbitals d \
    --spin up

Also allow separate panels for spin-up and spin-down bands.
C. Noncollinear and spin–orbit-coupled calculations
Correctly interpret the corresponding XML projection structure.
Distinguish between:

Orbital-projected total charge.
Orbital-resolved magnetization along x.
Orbital-resolved magnetization along y.
Orbital-resolved magnetization along z.
Do not interpret noncollinear magnetization components as independent spin-up or spin-down band channels.
For ordinary orbital projections, use the total charge contribution.
For future spin-texture analysis, design the internal data structures so the magnetization components can also be extracted.
Because magnetization components may have positive or negative values, use an appropriate signed-color representation rather than treating them as nonnegative marker weights.
First make total orbital-charge plotting reliable for spin–orbit-coupled calculations. Add magnetization-resolved plotting only after the underlying XML conventions have been verified.
8. Software architecture
Organize the implementation into modules with distinct responsibilities.
Suggested structure:
vasp_orbital_plotter/
    __init__.py
    parser.py
    projections.py
    kpath.py
    plotting.py
    cli.py
tests/
    test_parser.py
    test_projections.py
    test_kpath.py
    test_plotting.py
README.md
pyproject.toml

Use NumPy arrays for numerical data and Matplotlib for visualization.
Create an internal data representation that includes at least:

Eigenvalues indexed by spin channel, k-point, and band.
Projection weights indexed by spin or charge component, k-point, band, atom, and orbital.
Orbital names.
Chemical species and atom indices.
K-point coordinates.
Reciprocal lattice vectors.
Fermi energy.
Relevant calculation metadata.
Do not rely on positional assumptions about XML dimensions without validating them.
Avoid unnecessary dependencies.
pymatgen may be used as an independent reference for validation, but the primary parser should read the XML directly.
9. Validation and testing
Numerical validation is mandatory.
Create synthetic XML fixtures with known eigenvalues and projection weights to test:

Atom indexing.
Orbital indexing.
Summation over multiple atoms.
Summation over orbital groups.
Nonmagnetic calculations.
Collinear spin-polarized calculations.
Spin–orbit-coupled charge projections.
Missing orbital-projection data.
Incorrect XML dimensions.
Repeated or disconnected k-path segments.
For numerical validation, compare the extracted data against the existing pymatgen.io.vasp.outputs.Vasprun parser with projected-eigenvalue parsing enabled.
Specifically compare:
[
E_{n\mathbf{k}}^{\mathrm{tool}}
\quad\text{and}\quad
E_{n\mathbf{k}}^{\mathrm{pymatgen}},
]
and
[
W_{n\mathbf{k}}^{\mathrm{tool}}
\quad\text{and}\quad
W_{n\mathbf{k}}^{\mathrm{pymatgen}}.
]
The numerical differences should remain within the precision of the values stored in the XML.
Also compare at least one actual fat-band plot against VASPKIT or another independently generated reference plot, if reference output files are available.
Do not claim numerical agreement without performing the comparison.
10. Development stages
Stage 1 — Parsing and data validation
Implement the XML parser and verify the extracted eigenvalues, atomic indices, orbital labels, and projection weights.
Before starting visualization, demonstrate that the parser returns physically and dimensionally consistent arrays.
Stage 2 — Functional orbital-projected plotting
Implement the command-line interface and the two plotting modes.
Initially support:

Nonmagnetic calculations.
One or multiple selected atoms.
Individual orbitals and orbital groups.
Energy-window selection.
Automatic Fermi-energy alignment.
Fat-band and color-projected plots.
Stage 3 — Advanced functionality
Add:

Collinear spin polarization.
Spin–orbit-coupled calculations.
Optional magnetization-resolved plotting.
Multiple projection overlays.
Custom high-symmetry labels.
Publication-quality formatting.
Performance optimization for large XML files.
11. Final deliverables
Provide:

Complete, executable Python source code.
A dependency specification.
A README with installation and usage instructions.
Examples of plotting individual atoms, orbitals, and orbital groups.
Unit tests and numerical validation tests.
At least one generated example figure using test data.
A short explanation of the internal XML structure and how the code interprets the orbital-projection weights.
The tool must run directly from the command line without modification of its source code.
12. Instructions for your development process
Follow this workflow:

Inspect first: Examine the available project files and actual vasprun.xml structure before writing the parser. If no real XML file is available, use synthetic fixtures and explicitly state which features remain unverified.
Implement incrementally: Complete and test the parser before starting advanced plotting.
Test continuously: Run the relevant tests after every substantial implementation step.
Avoid unsupported assumptions: If VASP changes the XML layout between spin configurations or versions, identify the difference rather than forcing a single interpretation.
Report limitations honestly: Clearly distinguish implemented and verified functionality from features that require additional testing.
Keep the code maintainable: Use functions with clear responsibilities, type hints, descriptive names, and informative exceptions.
Do not stop at pseudocode: Deliver working software, a usable command-line interface, and tests.
Acceptance criterion: I should be able to place the tool in a directory containing a valid orbital-projected VASP vasprun.xml and execute a single command specifying an atom and orbital. I should then obtain an accurate, publication-quality orbital-projected band structure without using VASPKIT or PROCAR.
```
