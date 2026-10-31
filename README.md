# AltermagAgent

AltermagAgent reads a tree of VASP calculation folders and never changes them.
It rebuilds the project storyline: the runs in order, how each run came from the one before, and eight summary points with the files behind each one.
It also shows the branches of the project with suggestions for the next calculations.
It shows the storyline in a three-panel browser interface with a chat that answers only from what it found in your folders.

## Screenshots

- `screenshot_v1.png`: the synthetic CrSe demo project in the interface.
- `screenshot_v2_MnTe_real.png`: a real MnTe altermagnet project.

![CrSe demo](screenshot_v1.png)

![Real MnTe project](screenshot_v2_MnTe_real.png)

## Reproduce the demo in five commands

Run these after the install below, with the virtual environment active.

```bash
python demo/make_demo.py
altermag import demo/CrSe_demo --name CrSe
altermag summary CrSe
altermag bands CrSe bands_nosoc
altermag gui
```

Then open http://localhost:8080/?project=CrSe in your browser.

## Install

You need Python 3.12 or newer and git.
The steps below use `python3.14`; use `python3.12` or `python3.13` if that is what you have.
A first install downloads about 90 packages and takes about three minutes; a repeat install takes under a minute.

```bash
git clone https://github.com/joshbarron-36/altermagagent.git
cd altermagagent
python3.14 -m venv .venv
```

Activate the virtual environment, using the line for your system.

```bash
source .venv/bin/activate          # Linux and Mac
```

```
.venv\Scripts\activate             # Windows, Command Prompt or PowerShell
```

Then install the packages and the `altermag` command.

```bash
pip install -r requirements.txt
pip install -e .
```

On Windows, use `py -3.12 -m venv .venv` if `python3.14` is not found.

## The demo project

`python demo/make_demo.py` writes `demo/CrSe_demo`, a synthetic CrSe project with seven runs, so you need no VASP to try the tool.
It contains a relaxation, SOC and non-SOC scf and band runs, a failed run, a run that may still be running and a Wannier run.

```bash
python demo/make_demo.py
altermag import demo/CrSe_demo --name CrSe
altermag summary CrSe
```

## Commands

```bash
altermag import <folder> --name <project>    # read-only scan of a folder tree into the index
altermag rescan <project>                    # re-import a project from its stored folder; rebuilds runs and summary
altermag delete <project>                    # remove a project from the index; folders on disk are not touched
altermag summary <project>                   # eight summary points, branches, gaps, runs in order, procedure checks
altermag graph <project>                     # branches (structure, magnetic realization, SOC) and their runs as text
altermag suggest <project> <branch>          # suggestions for the next DFT calculations of one branch
altermag literature <project>                # published papers found by search, not evidence from the files
altermag bands <project> <run> [--emin -3] [--emax 3] [--style paper]
altermag structure <project> <run> [--supercell 1 1 1]
altermag compare-wannier <project> [path/to/wannier90_hr.dat] [--run-name <run>]
altermag spin-projection <project> [hr_file] [--run-name <run>] [--component x|y|z|all] [--order udud|uudd]
altermag gui [--port 8080]                   # start the browser interface
```

### import

`altermag import` walks the folder tree and stores one entry per VASP run in the index.
Importing again under the same name replaces the project.
Large outputs such as CHGCAR, WAVECAR and PROCAR are never opened; only their size and date are kept.

### summary

`altermag summary` prints the eight points: structures, magnetism, symmetries, SOC, Hamiltonian, transport, status and last update.
The Symmetries point gives the magnetic space group once per magnetic realization, with its magnetic point group; the crystal space group is in the Structures point.
The Magnetic point lists every magnetic realization found, for example "collinear (3 runs); coplanar (10 runs)", with net and site moments.
The spin splitting scan is no longer a summary point; it stays in the band figure output of `altermag bands`.
It then lists the gaps, for example a Wannier Hamiltonian without a transport step.
Each run has a phase (relax, scf, band, wannier, symmetrisation, projection, transport or other) and procedure checks.
The checks come from a written DFT procedure in `docs/detection_rules_from_procedure.md`, for example EDIFF and the number of job slots.

### graph and magnetic realizations

`altermag graph` groups the runs into branches and prints each branch with its runs in historical order and their parent links.
A branch is one structure, one magnetic realization and one SOC setting.
The structure match uses pymatgen StructureMatcher without volume rescaling and with a 1 percent lattice tolerance.
The magnetic realization comes from INCAR MAGMOM, or from the converged OUTCAR moments when MAGMOM is missing.
It is classified as collinear, coplanar, non-coplanar or non-magnetic, and two runs share it when their moment patterns agree up to one global rotation within 0.1 μB.
A branch is named after its first folder plus noSOC or SOC, such as `coplanar-SOC`, when that is unambiguous, otherwise `realization-1-SOC`.

### suggest

`altermag suggest <project> <branch>` prints the suggestions for the next DFT calculations of one branch.
It follows the physicist's loop in six steps, each marked done, partial, missing or not applicable with its evidence file:
(a) structure and moments, (b) symmetries, (c) converged SCF, (d) band structure and projections, (e) Wannierisation, (f) transport.
Step (b) runs "define MSG": spglib computes the magnetic space group of the branch's first converged run.
It uses the converged OUTCAR moments, or INCAR MAGMOM when the run has no OUTCAR, and tags the label "parsed, pending validation".
The magnetic space group counts as done when a findsym.log lies in the folder of the branch or above it.
Its label is then shown "from findsym.log", with the computed label next to it and "agree" or "differ".
The tool always writes the FINDSYM input (structure and moment vectors) to `~/.altermagagent/msg/<project>/<branch>/findsym.in`.
In the interface the same action is the "define MSG" button under the suggestions of the selected branch.
The suggestions end with the section "From the literature, not from these files"; without `VALENCY_TOKEN` it says "literature unavailable".
The first missing step becomes the recommended next calculation, written as one sentence with the INCAR tags to change.
The tool never states that a material is an altermagnet.

### bands

`altermag bands` reads `vasprun.xml` and `KPOINTS` of a band run and draws the bands relative to the Fermi level.
A labelled line-mode `KPOINTS` gives the high-symmetry labels; otherwise the x axis is the k-point index.
Spin-polarised runs without SOC show spin up in blue and spin down in red, with the largest splitting and its k-point.
SOC runs are drawn in one colour.
The `paper` style is white, serif and 300 dpi; the `screen` style is a larger figure for slides.

### structure

`altermag structure` reads POSCAR (or CONTCAR) and the MAGMOM line of INCAR and draws the atoms with the moments as arrows.
The arrows are the starting moments from INCAR, not the converged moments from OUTCAR.
The interactive HTML view loads 3Dmol.js from the internet; the PNG and SVG work offline.

### compare-wannier

`altermag compare-wannier` (alias `wannier-compare`) overlays a Wannier Hamiltonian on the VASP bands of its ancestor band run.
Without required input equality (POSCAR, INCAR, CHGCAR, CHG) it refuses, unless `--structure-only` asks for a preliminary overlay.
Errors are computed only with `--band-start N`, where Wannier band 1 is VASP band N, within -1 to +1 eV by default (`--window`).
They are MAE, RMSE, the 95th percentile, the maximum and the count of samples above 20 meV.
`--band-start auto` picks the offset by smallest error, and the figure, the JSON and step (e) then say "not verified".
Without seedname_wsvec.dat or a use_ws_distance=false declaration, the direct-HR fallback is named in the figure and the JSON.
Options: `--reference-run` (alias `--reference`), `--emin`, `--emax`, `--fermi-eV`, `--spin up|down`, `--out-dir` and `--no-ws-distance`.
Step (e) runs the comparison in preliminary mode with the automatic offset, so its line starts with "preliminary:".
The module was written by Diana Golovanova.

For Hamiltonian overlays, see the [Wannier comparison README](docs/wannier_comparison/README.md).
It documents ancestor discovery, matching inputs, preliminary comparisons and the generated figures.

### gui

`altermag gui` starts the interface on http://localhost:8080.
The left panel lists projects and has an Import folder button.
The middle panel shows the eight points, then a graph of the branches drawn like a git history, then the selected branch.
The graph has one lane per branch, one node per run in time order coloured by phase, and lines for parent links.
Click a node or a lane to select its branch: its header with the 3D structure view, its run cards in order and its suggestions.
Click a summary point to see its evidence files, and click a run card to see its files, parameters, moments and procedure checks.
Band runs have a Plot pill with energy inputs and a Redraw button, and runs with a structure file have a Structure pill.
The right panel is the chat.

## Orbital-projected bands from vasprun.xml

The package `vasp_orbital_plotter` draws fat-band and colour-projected band structures straight from a VASP `vasprun.xml`.
It needs no VASPKIT, PROCAR or EIGENVAL; the band run must have LORBIT = 11 so that the file holds the `<projected>` block.
It is installed with the steps above, as the command `vasp-orbital-bands`; `python -m vasp_orbital_plotter` works the same way.
Diana's specification is in `docs/orbital_plotter_spec.md`; stages 1 to 3 are built.

In a folder with `vasprun.xml`, one command plots the d bands of every Mn ion:

```bash
vasp-orbital-bands --element Mn --orbitals d --output mn_d.png
```

`--vasprun` defaults to `vasprun.xml` in the current folder.
`--output` is needed there, because the plotter never writes next to a vasprun.xml unless you name the file yourself.
Run elsewhere, it writes a name built from the selection, such as `bands_Mn_d_fatband.png`, into the current folder.

Selection examples (ions use VASP's numbering, starting at 1):

```bash
vasp-orbital-bands --vasprun run/vasprun.xml --list-atoms
vasp-orbital-bands --vasprun run/vasprun.xml --list-orbitals
vasp-orbital-bands --vasprun run/vasprun.xml --atoms 2 --orbitals pz --mode fatband --emin -3 --emax 3 --marker-scale 80 --output atom2_pz.pdf
vasp-orbital-bands --vasprun run/vasprun.xml --atoms 1,3-4 --orbitals dxy,dxz --output bands.svg
vasp-orbital-bands --vasprun run/vasprun.xml --element Te --orbitals p --mode color --output te_p.png
vasp-orbital-bands --vasprun run/vasprun.xml --select Mn:eg --select 3:pz --normalize --output mix.png
vasp-orbital-bands --vasprun run/vasprun.xml --element Mn --orbitals d --high-symmetry "1:G,41:M,81:K,121:G" --output mn_d.png
vasp-orbital-bands inspect --vasprun run/vasprun.xml
```

- Options use their full names: `--atoms`, `--element`, `--orbitals`, `--high-symmetry`, `--output`; abbreviations are refused.
- `--atom` and `--orbital` are accepted as other spellings of `--atoms` and `--orbitals`, so the specification's `--atom 2 --orbital pz` works.
- `--mode` is `fatband` or `color`.
- Orbital names follow the file's own labels, ignoring case, spaces, `_`, `^` and `-`; the shell letter is optional, so `dx2-y2` finds `x2-y2`.
- Groups: `s`, `p`, `d`, `f` (all components of that shell in the file), `t2g`, `eg` and `all`; an incomplete group is refused.
- The weight is the sum of VASP's stored values over the chosen ions and orbitals; `--normalize` divides by the total of each state.
- Energies are E − E_F by default; `--absolute` keeps absolute energies and `--efermi` overrides the file's Fermi level.
- `--emin` and `--emax` are relative to the Fermi level unless `--absolute` is given; the default window is −3 to 3 eV.
- If no band reaches the window, the figure is still written with empty axes, and a `note:` line says so.
- K-path distances are in 1/Å with the factor 2π, as in pymatgen and VASPKIT; separate path branches are never joined.
- Segments come from a line-mode KPOINTS file, else the vasprun's line-mode divisions, else repeated k-points, else one line.
- Labels come from `--high-symmetry "1:G,41:M"`, else KLABELS, else the vasprun's `<kpoints_labels>`, else KPOINTS comments.
- `inspect`, `--list-atoms` and `--list-orbitals` skip the projection numbers and never load Matplotlib, so they answer quickly.
- `--list-atoms` also works on a file without a `<projected>` block; `--list-orbitals` then says the file has no orbital labels.
- The format follows the extension: PNG, SVG or PDF; `--dpi`, `--figsize`, `--font-size`, `--cmap`, `--cmin` and `--cmax` adjust the look.
- A bad option, atom or orbital stops with one `error:` line and exit code 2; an unknown orbital lists the available ones.

### Spin, magnetisation, overlays and formatting (stage 3)

```bash
vasp-orbital-bands --vasprun cr/vasprun.xml --element Cr --orbitals d --spin up --output cr_d_up.png
vasp-orbital-bands --vasprun cr/vasprun.xml --element Cr --orbitals d --spin both --mode color --output cr_d_spins.png
vasp-orbital-bands --vasprun run/vasprun.xml --element Mn --orbitals d --component mz --mode color --output mn_d_mz.png
vasp-orbital-bands --vasprun run/vasprun.xml --select Mn:d --select Te:p --overlay --output overlay.png
vasp-orbital-bands --vasprun run/vasprun.xml --select Mn:d --select Te:p --overlay --label-series "Mn 3d" --label-series "Te 5p" --mode color --output overlay_colour.pdf
vasp-orbital-bands --vasprun run/vasprun.xml --element Mn --orbitals d --label-map 'C0=C$_0$,S=X' --style talk --despine --output talk.svg
vasp-orbital-bands --vasprun big/vasprun.xml --element Cr --orbitals d --kpoint-stride 2 --float32 --output big.png
```

- `--spin up|down|both` is for collinear ISPIN = 2 files; each spin uses its own weights and its own eigenvalues.
- The default for collinear files is `both`: two panels side by side sharing the energy axis, one colour scale.
- Without `--figsize` each panel is 6 by 4.5 inches, so two panels make a 12 inch wide figure; an explicit `--figsize` is the whole figure.
- `--component total|mx|my|mz` is for noncollinear and SOC files; the default is `total`, the charge projection.
- mx, my and mz are magnetisation blocks, never spin channels; they are signed and drawn with signed colours.
- Colour mode uses a diverging map (`--cmap-signed`, default RdBu_r) with a grey middle at zero, symmetric about the largest |value|.
- Fat-band mode draws red markers for positive and blue for negative values, with area proportional to |value|.
- Each mx, my or mz figure has a caption saying "sign convention unverified" and the run's LSORBIT; the sign is as read from the file.
- Noncollinear runs without SOC (LSORBIT = F) are plotted like SOC runs, with a caption giving LSORBIT; this is an assumption until Diana confirms.
- `--spin` on a nonmagnetic or noncollinear file, or `--component mx` on a nonmagnetic or collinear file, stops with an `error:` line naming the right option.
- `--normalize` and `--overlay` are refused for mx, my and mz.
- Without `--overlay`, repeated `--select` terms are summed; with it, each term is its own series, colour and legend entry (at most 6).
- In colour mode an overlay colours each state by its dominant series, faded to grey as that weight falls to zero; there is no colour bar.
- `--label-series` (repeatable, one per `--select`) names the series in the legend.
- `--label-map 'G=Γ,C0=C$_0$'` renames path labels (single quotes keep the shell off the `$` signs); keys match the printed or the typed name, and mathtext renders.
- Entries are split at commas, except inside `$...$` mathtext or braces, so `'X=X$_{1,2}$,G=Γ'` renames two labels.
- `--no-labels` keeps the high-symmetry lines and drops their names; an unknown `--label-map` key lists the labels on the path.
- `--style paper` (default) is the stage 2 look; `--style talk` scales fonts by 1.4, lines by 1.5 and marker area by 2.
- `--rcparams-file style.json` applies a JSON object of Matplotlib rcParams; `--despine` removes the top and right frame lines.
- `--title` replaces the title and `--no-title` removes it; ticks point inwards on the bands and the colour bar.
- SVG and PDF output is fully vector, with editable TrueType text in PDF; `--rasterize` embeds markers and lines as an image instead.
- `--kpoint-stride N` plots every Nth k-point of each segment and always keeps both segment ends.
- Weights are summed only for bands that reach the energy window and only at the kept k-points; `--float32` halves the memory of the projections.

Mn d fat bands and Te p colour projection of Diana's MnTe `9.2-band` run (SOC, total charge):

![Mn d fat bands of MnTe](docs/figures/orbital_plotter/mnte_Mn_d_fatband.png)

![Te p colour projection of MnTe](docs/figures/orbital_plotter/mnte_Te_p_color.png)

Signed mz projection of the Mn d states (sign convention unverified), and an overlay of Mn d and Te p, from the same run:

![Mn d mz projection of MnTe](docs/figures/orbital_plotter/mnte_Mn_d_mz_color.png)

![Mn d and Te p overlay of MnTe](docs/figures/orbital_plotter/mnte_Mn_d_Te_p_overlay.png)

## Open a project directly

Add `?project=NAME` to the address, for example http://localhost:8080/?project=MnTe.
Without it, the interface opens the first project in the index.

## The chat

The chat answers only from the indexed summary and runs, and names the run and file behind each fact.
It can draw a band figure when you ask for one.
It picks its model in this order: an Anthropic key, then a local Ollama server, then a fixed notice.

### With an Anthropic API key

```bash
export ANTHROPIC_API_KEY=sk-ant-...     # Linux and Mac
```

```
set ANTHROPIC_API_KEY=sk-ant-...        # Windows Command Prompt
$env:ANTHROPIC_API_KEY="sk-ant-..."     # Windows PowerShell
```

Set the key in the same terminal before `altermag gui`.
The key is read from the environment only and is never stored in the code.

### With Ollama, no key needed

1. Install Ollama from https://ollama.com/download.
2. Download a model once: `ollama pull llama3.1:8b`.
3. Keep `ollama serve` running; the Ollama desktop app does this for you.
4. Start `altermag gui` and ask in the chat panel.

The chat finds the server at `OLLAMA_HOST` (default `http://localhost:11434`).
To use another model, pull it and set `ALTERMAG_OLLAMA_MODEL`, for example `export ALTERMAG_OLLAMA_MODEL=qwen2.5:7b`.
An 8 billion parameter model needs about 6 GB of free memory.
Drawing figures from the chat works with llama3.1 and qwen models; other models still answer and point you to the Redraw button.

## Literature, kept apart from the files

A separate section, "From the literature, not from these files", follows the suggestions of the selected branch in the interface.
It also follows the branches in `altermag summary` and the suggestions in `altermag suggest`; `altermag literature <project>` prints only that section.
It lists published papers found by a search of the Valency literature service (https://labs.valency.io/mcp), not evidence read from your folder.
The suggestions never change because of these papers.
For the project formula it runs up to two searches, for example "MnTe altermagnet spin splitting" and "MnTe magnetic order P6_3/mmc".
Each search shows up to three papers with title, first author, year and link.
Results are cached per search for 24 hours in `~/.altermagagent/literature_cache.json`.
The chat can call the same search as a tool; its answers are labelled "from the literature" and cite title and year.

The search needs a Valency access token in the environment variable `VALENCY_TOKEN`.
Without it, or when the service fails or refuses the token, the section prints "literature unavailable" with the reason, and everything else works as before.

How to get a token:
1. Valency uses OAuth with a browser login; there is no key to copy from a web page.
2. Register an OAuth client through the service's registration endpoint, listed at https://labs.valency.io/.well-known/oauth-authorization-server.
3. Open the authorization URL in a browser with PKCE, log in, and let it redirect to a small local callback server.
4. Exchange the returned code at the token endpoint for an access token and a refresh token.
5. Export the access token as `VALENCY_TOKEN`, for example from a private file such as `~/.altermagagent/env`, then `source` that file.
6. When the token expires, exchange the refresh token for a new one at the same token endpoint.

On Josh's laptop a small login script outside this repository does these steps and writes `~/.altermagagent/env`; run its refresh option when the token expires.

## Notebook walkthrough with marimo

`notebooks/altermagagent_walkthrough.py` is a marimo notebook that shows the whole demo and analysis in one file.
It shows the summary points, the per-branch suggestions, the runs in order, the band and structure figures and the benchmark.
Two sliders set the energy window of the band figure, and a dropdown picks the project; no API key is needed.
If the index is empty, the notebook imports `demo/CrSe_demo` first.

```bash
marimo edit notebooks/altermagagent_walkthrough.py   # editable view with the code
marimo run notebooks/altermagagent_walkthrough.py    # app view without the code
```

## Where the index and figures live

Everything the tool writes goes to `~/.altermagagent`:

- `index.sqlite`: the index of all imported projects.
- `figures/<project>/`: band and structure figures as PNG, SVG, HTML and JSON.
- `agent_calls.log`: one line per chat call with provider, model and token counts.

Set `ALTERMAG_HOME` to use another folder, for example `export ALTERMAG_HOME=/tmp/altermag_test`.

## Your data folders are never written

The scan, the figures and the chat only read your VASP folders.
Nothing is ever written inside a scanned folder, and the import refuses a folder that contains the index.
You can point the tool at a shared or backed-up project tree without risk to it.

## Known limits

- Measured accuracy on 20 synthetic CrSe trees: phase 100 %, status 100 %, parent 98.5 %; files copied with `cp -p` can still reorder runs, see `docs/accuracy.md`.
- Magnetic space groups and spin groups are not computed yet; the summary says "pending".
- Reading folders on a cluster over SSH is not available yet; copy or mount the folders locally.
- Small Ollama models (3 billion parameters) give weak answers and sometimes draw figures nobody asked for.
- Structure arrows show the INCAR starting moments, not the OUTCAR moments.
- The 3D structure view needs an internet connection.

## Software and references

AltermagAgent reads results that were computed earlier and adds its own analysis on top.
The references below are split so that a reader can tell which is which.

How to cite:

- Using AltermagAgent to inspect existing results: cite AltermagAgent and the original computational work.
- Publishing results that rely on a specific method or package: also cite that method or package below.

### Scientific data provenance

Software and methods that produced the CrSe and MnTe results; AltermagAgent only reads their outputs.

The density functional theory (DFT) calculations for CrSe and MnTe were performed by Diana Golovanova with the Vienna ab initio Simulation Package (VASP) [1, 2].
The electron-ion interaction was described with the projector augmented-wave (PAW) method [2].
Exchange and correlation were treated in the generalized gradient approximation with the Perdew-Burke-Ernzerhof (PBE) functional [3].
Maximally localized Wannier functions [4] were constructed with Wannier90 3.1 [5], with the electronic-structure data passed from VASP through its VASP2WANNIER90 interface.
The resulting Wannier tight-binding Hamiltonians were used for the later analysis.

1. G. Kresse, J. Furthmüller, Phys. Rev. B 54, 11169 (1996), https://doi.org/10.1103/PhysRevB.54.11169
2. G. Kresse, D. Joubert, Phys. Rev. B 59, 1758 (1999), https://doi.org/10.1103/PhysRevB.59.1758
3. J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996), https://doi.org/10.1103/PhysRevLett.77.3865
4. N. Marzari, D. Vanderbilt, Phys. Rev. B 56, 12847 (1997), https://doi.org/10.1103/PhysRevB.56.12847
5. G. Pizzi et al., J. Phys.: Condens. Matter 32, 165902 (2020), https://doi.org/10.1088/1361-648X/ab51ff

Wannier Hamiltonian symmetrisation, used to generate input data, not run by AltermagAgent:

- wannhr_symm_Mag: A tool for symmetrization of magnetic WannierTB, https://github.com/quanshengwu/wannier_tools/tree/master/utility/wannhr_symm_Mag, Changming Yue.
- wannhr_symm: A tool for symmetrization of non-magnetic WannierTB, https://github.com/quanshengwu/wannier_tools/tree/master/utility/wannhr_symm/, Changming Yue.
- Both are part of WannierTools: Q. Wu, S. Zhang, H.-F. Song, M. Troyer, A. A. Soluyanov, "WannierTools: An open-source software package for novel topological materials", Comput. Phys. Commun. 224, 405-416 (2018), https://doi.org/10.1016/j.cpc.2017.09.033

Magnetic symmetry files (findsym.in, findsym.log), generated with Symmetr, not run by AltermagAgent:

- Symmetr: J. Železný, "Symmetr: a Python package for determining symmetry properties of crystals", Comput. Phys. Commun. 328, 110309 (2026), https://doi.org/10.1016/j.cpc.2026.110309, code at https://github.com/zeleznyj/symmetr (preprint arXiv:2602.21034)
- Symmetr internally uses FINDSYM, developed by Harold Stokes and Dorian Hatch, to determine crystal symmetry:
  H. T. Stokes and D. M. Hatch, "FINDSYM: program for identifying the space-group symmetry of a crystal", J. Appl. Cryst. 38, 237-238 (2005), https://doi.org/10.1107/S0021889804031528
  H. T. Stokes, D. M. Hatch, and B. J. Campbell, FINDSYM, ISOTROPY Software Suite, iso.byu.edu.

### Analysis dependencies

Packages whose analysis AltermagAgent runs itself.

- spglib (crystal and magnetic space groups): A. Togo, K. Shinohara, I. Tanaka, Sci. Technol. Adv. Mater.: Methods 4, 2384822 (2024), https://doi.org/10.1080/27660400.2024.2384822
- spglib magnetic space group search: K. Shinohara, A. Togo, I. Tanaka, Acta Cryst. A 79, 390 (2023), https://doi.org/10.1107/S2053273323005016
- spinspg (spin groups of collinear non-SOC runs): K. Shinohara et al., Acta Cryst. A 80, 94 (2024), https://doi.org/10.1107/S2053273323009257
- pymatgen (structures, band files and matching): S. P. Ong et al., Comput. Mater. Sci. 68, 314 (2013), https://doi.org/10.1016/j.commatsci.2012.10.028
- Wannier band comparison: an internal AltermagAgent module written by Diana Golovanova, see `docs/wannier_comparison`.

### Interface and external services

Libraries and services for display, interaction and literature search.

- NiceGUI: F. Schindler, R. Trappe, "NiceGUI: Web-based user interfaces with Python. The nice way.", https://doi.org/10.5281/zenodo.7785516
- marimo: A. Agrawal, M. Scolnick, "marimo - an open-source reactive notebook for Python", https://doi.org/10.5281/zenodo.12735329
- 3Dmol.js (3D structure view): N. Rego, D. Koes, Bioinformatics 31, 1322 (2015), https://doi.org/10.1093/bioinformatics/btu829
- Matplotlib: J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering 9(3), 90-95 (2007), https://doi.org/10.1109/MCSE.2007.55
- Plotly: Plotly Technologies Inc., Collaborative Data Science, Montréal, QC (2015), https://plot.ly
- Valency: the Valency MCP platform for the literature search, https://valency.io

## License

See `LICENSE`.
