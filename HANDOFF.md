# Handoff note

## Orbital-projected band plotter, stage 3, 30 October 2026

- Branch `orbital-plotter-stage3`, off main e6800b5, in the worktree `worktrees/orbital-plotter-stage3`. Committed on the branch after review and tests; merge to main pending.
- Plan: `orbital-plotter-plan_v3.html` next to the repository (block choice, signed plots, overlays, labels, styles, large files, assumptions).
- Nothing in `altermagagent/` changed; no new dependency.
- `--spin up|down|both` for collinear ISPIN = 2 files: weights from the chosen spin block, eigenvalues from the same spin block.
- The collinear default is now `both` (two panels side by side, shared energy axis, one colour scale), with a note.
- Without `--figsize` the figure is 6 by 4.5 inches per panel; an explicit `--figsize` is the whole figure.
- `--component total|mx|my|mz` for noncollinear files, default total; mx, my and mz are never treated as spin channels.
- Colour mode for mx, my, mz: diverging map (`--cmap-signed`, default RdBu_r) with a grey middle, centred on zero, symmetric about the largest |value|.
- The grey middle keeps zero-magnetisation states visible on white paper; the colour bar reads "mz projection (sign convention unverified)".
- Fat-band mode for mx, my, mz: red markers for positive, blue for negative, area proportional to |value|, legend per sign.
- Colour mode for mx, my, mz draws no legend entry, even with `--legend`: a flat-colour line would hide the sign, and the colour bar names the quantity.
- Every mx, my, mz figure has a caption: block as read from vasprun.xml, "sign convention unverified", and LSORBIT; it is also printed as a note.
- Refusals with one `error:` line, exit 2: `--spin` on nonmagnetic or noncollinear files (the latter points to `--component`).
- Also refused: `--component` other than total on nonmagnetic or collinear files, and `--normalize` or `--overlay` with mx, my, mz.
- `--overlay`: each `--select` term is its own series with its own colour and legend entry; `--label-series` (repeatable) names them; at most 6.
- Series colours in fixed order: the stage 2 blue, vermillion, green, orange, purple, sky blue; the order passes the colour-vision separation check.
- Colour-mode overlays colour each state by its dominant series, faded to grey as that weight falls to zero, with a legend and no colour bar.
- Judgment call: dominant colour rather than refusal, because "which selection dominates" is the usual question of a character plot.
- `--label-map 'G=Γ,C0=C$_0$'` renames labels (keys match the printed or typed name); unknown keys and bad mathtext are refused before parsing.
- Commas inside `$...$` mathtext or inside braces do not split `--label-map` entries, so `'X=X$_{1,2}$,G=Γ'` gives two entries.
- `--no-labels` keeps the lines and drops the names. C₀ and Σ₀ from KLABELS still render as subscripts.
- `--style paper` (default, the stage 2 sizes) or `talk` (fonts × 1.4, lines × 1.5, marker area × 2); scaling keeps an explicit `--font-size` meaningful.
- `--rcparams-file` applies a JSON object of rcParams, checked by Matplotlib; `--despine` removes the top and right frame lines.
- Ticks point inwards on the bands and the colour bar; PDF and PS embed TrueType text; SVG and PDF are fully vector unless `--rasterize`.
- Figures now use Matplotlib's constrained layout, so the margins differ slightly; the two stage 2 docs figures were regenerated for that.
- Legends have a light frame and take the best free corner; long captions wrap to the figure width.
- `--kpoint-stride N` keeps every Nth k-point of each segment plus both segment ends; branches and segment starts are renumbered.
- Weights are summed with `numpy.einsum` on views cut to the band window, the kept k-points and the range of selected ions; nothing is copied.
- The kept k-points are read as a few constant-step slices (one per segment plus the ends), so `--kpoint-stride` cuts the sum itself.
- `--float32` stores the projections in single precision (parser option `float32=True`); the default stays float64.
- The parser's row XPath now returns plain strings (`smart_strings=False`), about 5 to 10 percent faster on large files.
- New file `tests/test_orbital_stage3.py`, 38 tests, including CLI exit-code tests for every new refusal (20 cases).
- Full check: `PYTHONPATH=. python -m unittest discover -s tests` runs 206 tests, OK, one skip (168 before stage 3, 201 before the review fixes).
- Verified on MnTe `9.2-band`, read-only: total per state, summed over ions and orbitals, lies in 0.303 to 0.996.
- There, |m| = √(mx² + my² + mz²) per state never exceeds the total (largest excess 0.0); per (ion, orbital) entry |m| equals the total within 1.2e-4.
- That is what a projection of one spinor state gives, so the four blocks behave as charge plus a magnetisation vector.
- Summed over occupied states along the path, Mn1 gives (mx, my, mz) = (3.88, 2.24, 0.001) and Mn2 the opposite; MAGMOM is (3.98, 2.30, 0) and (−3.98, −2.30, 0), SAXIS 0 0 1.
- So the sign and axes of mx, my agree with MAGMOM on this run; that is a consistency check, not a proof of VASP's convention.
- Timings with a warm file cache on this laptop, CrSe `2-SOC/2-bands` (354 MB, 640 k-points, 120 bands, Cr d): parse 6.3 to 6.4 s before, 5.6 to 5.8 s after.
- Parse peak memory 293 MB before and after, 1.10 times the 265 MB projection array; with `--float32` the array is 133 MB and the peak 166 MB.
- Select, k-path and plot: 1.02 s before, 1.00 s after, 0.67 s with `--kpoint-stride 4`; the weight sum itself takes 2 to 4 ms.
- After the code-quality pass (parse cached, best of 8 runs, same laptop): 0.95 s, and 0.53 s with `--kpoint-stride 4`.
- Before that pass, measured the same way in the same session: 0.90 s and 0.51 s, so the difference is within the noise of this laptop.
- The weight sum of Cr d alone: 2.1 ms before and 1.9 ms after; with `--kpoint-stride 4`, 2.2 ms before and 0.4 ms after.
- With `--normalize` and `--kpoint-stride 4`: 5.0 ms before and 0.7 ms after; the state totals are summed once per panel.
- The plot phase adds about 70 MB, mostly Matplotlib itself (55 MB on a tiny file).
- The `--float32` figure differs from the float64 one in 9 of 2.43 million pixels, by 1/255 each.
- CrSe `2-SOC/1-scf` (785 MB, parse only): 14.0 s before, 13.0 s after, peak 614 MB for a 597 MB array; with `--float32` 13.0 s, 299 MB array, 329 MB peak.
- The first cold-cache parse of the 354 MB file took 9.5 s; parsing stays bound by lxml and number conversion, which a new dependency would be needed to beat.
- New figures at 150 dpi from `9.2-band`: `docs/figures/orbital_plotter/mnte_Mn_d_mz_color.png` (198 KB) and `mnte_Mn_d_Te_p_overlay.png` (273 KB).
- Regenerated: `mnte_Mn_d_fatband.png` (229 KB) and `mnte_Te_p_color.png` (155 KB). README section extended with the stage 3 options and examples.
- After the code-quality pass all four figures were regenerated with the same commands; they are byte-identical to the files above.
- Code-quality pass: the stage 2 entry points `select_component` and `plot_bands` are removed; their tests use `choose_blocks` and `plot_panels`.
- One table, `COMPATIBILITY` in `projections.py`, keyed by (layout, channel), holds the blocks, the sign and every refusal text.
- `choose_blocks` (now also given `normalize` and `overlay`) and the command line's early check both read it; the CLI checks before parsing.
- `BlockChoice` carries `tag` (output name), `quantity` (colour bar and legend text) and `eigen(block)`; `Block` keeps only the projection index and heading.
- The command line cuts the panels to the band window and shifts energies once; `plot_panels` takes energies relative to `zero` and draws them as given.
- Review fixes: every refusal that needs no file now comes before the parse, in the command line's `_check_flags`; tests mock the parser to prove it.
- Those are the default output in the vasprun's own folder without `--output`, `--kpoint-stride` below 1, and `--overlay` without `--select`.
- Also before the parse: more than 6 `--select` terms with `--overlay`, and a `--label-series` count that differs from the `--select` count.
- So the earlier regression, where a stride below 1 was refused only after parsing, is resolved; `stride_indices` and `build_series` keep their checks for library callers.
- `plot_panels` refuses an empty panel list, a panel without series and a panel whose k-point count differs from the path, before drawing.
- The default figure size is computed only in `plot_panels`; `_figure` takes the final size. All four docs figures re-render byte-identical.
- Every function in `vasp_orbital_plotter/` has exactly one header comment line; longer explanations moved into docstrings.
- Styles are one dict per preset (`scale` and `rc`) in `options.py`; `validate` reads the rcParams file once and keeps it on the options.
- `SERIES_COLOURS` moved to `options.py` and `MAX_SERIES` is its length; series colours are assigned only with `--overlay`.
- Assumption: mx, my, mz are used as read from the file; captions say "sign convention unverified".
- Assumption: noncollinear runs with LSORBIT = F are plotted like SOC runs (same four blocks); their total-charge figures carry a caption giving LSORBIT.
- Assumption: collinear plots are verified on the synthetic fixture only; no real ISPIN = 2 band file is in the data.
- Unverified: the mx, my, mz convention against VASP's documentation, and its behaviour with a SAXIS other than 0 0 1.
- Unverified: collinear and nonmagnetic plots on real files; fat-band weights against VASPKIT (that needs PROCAR, which the rules forbid).
- Question for Diana: is marker area proportional to the raw weight what you want, or should the default follow VASPKIT's scaling?
- Question for Diana: should the default output in a run folder be refused as now, or written there, which breaks the read-only rule for data folders?
- Question for Diana: for collinear runs, is two panels (both spins) the right default, or should one spin be the default?
- Question for Diana: do you prefer another default energy window, figure size or colour map for papers?
- Question for Diana: can you confirm the sign and axis convention of mx, my, mz, for example relative to SAXIS?
- Question for Diana: should noncollinear runs without SOC (LSORBIT = F) be treated like SOC runs, as assumed here?
- Question for Diana: do you have a real collinear ISPIN = 2 band run with LORBIT = 11 we may use for checks?

## Orbital-projected band plotter, stage 2, 30 October 2026

- Branch `orbital-plotter-stage2`, off main c54d570, in the worktree `worktrees/orbital-plotter-stage2`. Merged into main on 30 October 2026 after review and tests.
- Plan: `orbital-plotter-plan_v2.html` next to the repository (selection grammar, label sources, plot modes, validation, open points).
- Nothing in `altermagagent/` changed; the plotter reuses helpers of `altermagagent/bands.py` (pretty_label, build_ticks, read_klabels, klabels_ticks, mpl_label, segments, slug).
- No new dependency: Matplotlib, NumPy and lxml were already pinned.
- `projections.py`: ions by 1-based number, range or element; the 0-based conversion goes through the parser's `ion_index`.
- Orbitals are matched to the file's `<field>` labels by one rule: case, spaces, `_`, `^` and `-` are ignored, and the shell letter is optional.
- So `dx2-y2` finds VASP's `x2-y2` and `z2` finds `dz2`; the old alias table (`dx2`, `d3z2-r2`, `dzx`) is gone.
- Groups s, p, d, f come from the file labels (a label without a shell letter, such as `x2-y2`, continues the shell before it).
- A shell group needs all 1, 3, 5 or 7 components; t2g and eg are the only hand-written groups; `all` is every column.
- An unknown orbital stops with one `error:` line listing the file's orbitals and usable groups, exit code 2.
- `--select ATOMS:ORBITALS` (repeatable) sums several combinations; a pair named twice is counted once.
- Weight = plain sum of the stored values over the chosen (ion, orbital) pairs; raw by default, `--normalize` divides by the state's total over all ions and orbitals.
- Projection block: "total" for nonmagnetic and noncollinear files; spin up (block 0, eigenvalue block 0) for collinear files, with a printed note.
- `kpath.py`: distance = 2π × rec_basis, in 1/Å, the same convention as pymatgen's reciprocal_lattice used by `wannier_compare.reference_path`.
- Segments from a matching line-mode KPOINTS next to vasprun.xml, else the vasprun's line-mode divisions, else exactly repeated k-points.
- Otherwise the path is one continuous line with a note; the former jump heuristic (steps over ten times the median) is removed.
- At a segment start the distance does not grow and no line joins the two branches.
- Labels: `--high-symmetry "1:G,41:M"`, else KLABELS, else `<kpoints_labels>` converted from corner numbers, else KPOINTS comments, else none.
- The label sources are one ordered list; the first one that gives ticks wins, so the reported source always matches the ticks.
- `plotting.py`: fat-band mode (marker area = scale × weight, default 60 points² per unit weight) and colour mode with a colour bar.
- Colour mode colours each half-step of a band by its nearest computed state, so no two states' colours are blended; the default map is Blues without its white end.
- Colour mode draws no grey underlay, since the coloured lines cover it completely; the regenerated Te p figure is pixel-identical without it.
- Fat-band markers are drawn only for states inside the energy window with a positive weight.
- Defaults: E − E_F, window −3 to 3 eV, dashed zero line, 6 × 4.5 inches, 300 dpi, font size 12; PNG, SVG or PDF by extension.
- `cli.py`: `inspect` (`--vasprun`, `--list-atoms`, `--list-orbitals`), `plot`, and the plot flags directly on the main command; `python -m vasp_orbital_plotter` works too.
- `--atom` and `--orbital` are accepted as other spellings of `--atoms` and `--orbitals`, matching the specification's example.
- No other aliases (`--elements`, `--normalise`, `--labels`, `-o`) and no prefix abbreviations.
- `--mode` accepts `fatband` or `color` only; `inspect` takes the file with `--vasprun` only, no longer as a positional argument.
- Defaults live once, in `PlotOptions` (`options.py`); the command line reads them, and `plot_bands` runs `PlotOptions.validate` once.
- `validate` also refuses a font size of 0 or less, a negative line width, and `--cmin` not below `--cmax`.
- `--emin` and `--emax` are relative to the Fermi level unless `--absolute`; a window no band reaches gives empty axes, a `note:` line and exit 0.
- `inspect`, `--list-atoms` and `--list-orbitals` parse with `with_projections=False`: orbital labels are kept, projection numbers are freed unread.
- The listings do not require a `<projected>` block; `--list-orbitals` on such a file gives one `error:` line saying it has no orbital labels.
- Those paths never import Matplotlib, `plotting.py` or `kpath.py`; on the MnTe run the light parse takes about half the time.
- `--vasprun` defaults to `vasprun.xml`. Without `--output`, the file goes to the current folder, but never next to the vasprun.xml: there the command stops and asks for `--output`.
- Exit codes: 0 for success, 2 for every user error, printed as one `error:` line without a traceback (argparse errors included).
- Only the project's own errors and a missing, folder or unreadable input path become `error:` lines; any other error keeps its traceback.
- An empty `<field></field>` orbital label in `<projected>` stops the parser with a VasprunParseError naming the column.
- Verified on MnTe `9.2-band`, read-only: band energies agree with VASPKIT's BAND.dat to 5.0e-5 eV at most (100 k-points × 60 bands).
- That is the rounding of vasprun.xml, which stores 4 decimals, while BAND.dat comes from EIGENVAL with 6; REFORMATTED_BAND.dat (4 decimals) agrees to 9.6e-5 eV.
- K-distances agree with BAND.dat to 4.9e-6 1/Å (its 5-decimal rounding); path length 4.274183 here and 4.27418 there, so VASPKIT also uses 1/Å with 2π.
- The KLABELS positions (0.000, 0.875, 1.380, 2.390, 3.400, 4.274) sit at our segment ends to within their 3-decimal rounding.
- Labels on that run come from KLABELS: Γ, Y, C₀, Σ₀, Γ, S; segments from KPOINTS, five of 20 points, continuous, no breaks.
- Example figures from that run: `worktrees/orbital-plotter-stage2-figures/` (Mn d fat bands, Te p colour, Mn eg wide window as PDF, Mn1 t2g normalised as SVG).
- Two of them are in the repository: `docs/figures/orbital_plotter/mnte_Mn_d_fatband.png` (218 KB) and `mnte_Te_p_color.png` (146 KB), at 150 dpi; the README shows both.
- README: new section "Orbital-projected bands from vasprun.xml". Spec: the glued `--list-orbitals` line is split; the header now says six line breaks restored.
- Tests: `test_orbital_projections.py` (20), `test_orbital_kpath.py` (21, 4 on the real run), `test_orbital_plotting.py` (17, 1 on the real run).
- Full check: `PYTHONPATH=. python -m unittest discover -s tests` runs 168 tests, OK, one skip (108 before stage 2). Tests write only into temporary folders.
- Unverified: the fat-band weights against VASPKIT, because that needs PROCAR, which the repository rule forbids reading; only energies and distances were compared.
- Unverified: nonmagnetic and collinear plots on real files (none in the data), and real KPOINTS paths with breaks (the MnTe paths are continuous).
- Unverified: Cartesian line-mode KPOINTS files; their corners are not checked against the k-points, only the count is.
- Question for Diana: is marker area proportional to the raw weight what you want, or should the default follow VASPKIT's scaling?
- Question for Diana: should the default output in a run folder be refused as now, or written there, which breaks the read-only rule for data folders?
- Question for Diana: for collinear runs, is spin up a sensible stage 2 default until stage 3 adds the spin choice and two panels?
- Question for Diana: do you prefer another default energy window, figure size or colour map for papers?

## Orbital-projected band plotter, stage 1, 30 October 2026

- Branch `orbital-plotter`, off main af4227f. Not committed yet; Josh's session reviews before any commit.
- New package `vasp_orbital_plotter/`, a sibling of `altermagagent/`; nothing in `altermagagent/` changed.
- Diana's specification is now in `docs/orbital_plotter_spec.md`: her text, lightly reformatted (one long sentence split, one table rebuilt, five glued line breaks restored).
- The duplicated sections 9 to 12 of the specification are removed.
- Stage 1 only: the parser in `vasp_orbital_plotter/parser.py` and the command `vasp-orbital-bands inspect <vasprun.xml>`.
- `projections.py`, `kpath.py` and `plotting.py` are stubs that raise NotImplementedError naming stage 2.
- `pyproject.toml` lists the new package and the console script; `lxml==6.1.3` is pinned there and in `requirements.txt`.
- The parser streams the file once with lxml iterparse and frees each block after reading it; .xml.gz is accepted.
- It reads ISPIN, LSORBIT, LNONCOLLINEAR, LORBIT and NBANDS from `<parameters>`, the values VASP applied.
- It reads the Fermi energy, structure, reciprocal lattice, k-points and weights, eigenvalues, occupations and `<projected>` weights.
- Projection array shape: (component, kpoint, band, ion, orbital); eigenvalues and occupations: (spin, kpoint, band).
- Components: "total" without spin; "up", "down" for ISPIN = 2; "total", "mx", "my", "mz" for noncollinear runs.
- Orbital labels come from the `<field>` entries of the file, in file order; nothing is hard-coded.
- Users give ions in VASP's 1-based numbering; `ion_index(n)` returns n - 1 and refuses numbers out of range.
- Without `<projected>`, the parser stops and asks for LORBIT = 11; it never uses `<dos><partial>` in its place.
- Spin-block counts are read from the file and checked against the layout's labels; `<dimension>` order may vary and is transposed.
- Wrong band, k-point or spin-block counts, unknown dimension names and truncated XML stop with a message that names the place.
- With several `<calculation>` blocks, all band data comes from the last one only, and a note says so.
- If only an earlier `<calculation>` has `<projected>`, the parser stops instead of pairing those weights with other eigenvalues.
- Parse errors never escape as tracebacks: empty `<v>` rows, a missing ion `<set>` and similar faults raise VasprunParseError.
- The command prints such errors as one `error:` line and exits with code 2.
- Progress also works for .xml.gz files; there the percentage counts compressed bytes read from disk.
- `ion_index` refuses True and False, which Python would otherwise accept as the integers 1 and 0.
- The real data, read-only: 13 vasprun.xml files, 19 MB to 785 MB, all with `<projected>` and LORBIT = 11.
- All 13 are noncollinear (LNONCOLLINEAR = T), including the CrSe "noSOC" runs, which have LSORBIT = F.
- Their INCAR says ISPIN = 2, but VASP applied ISPIN = 1; the parser uses the applied value and adds a note.
- Every file has 4 projection blocks (total, mx, my, mz) and 1 eigenvalue block.
- Orbital fields in every file: s, py, pz, px, dxy, dyz, dz2, dxz, x2-y2 (VASP 5.4.4 and 6.3.2 write the same order).
- VASP 6.3.2 also writes eigenvalues inside `<projected>`; the parser reads them and notes any difference from the main block.
- The MnTe band runs have `<kpoints_labels>`; their numbers are KPOINTS corner lines, not k-point indices, so they are kept raw.
- Verified against pymatgen 2026.9.24 `Vasprun(parse_projected_eigen=True)`: MnTe `9.2-band` and `9.4-band`, full files.
- Result there: eigenvalues, occupations, Fermi energy, k-points, k-weights, species, total weights and mx, my, mz all differ by 0.0.
- Also verified: CrSe `2-SOC/2-bands` and `1-noSOC/030_bands`, sliced to their first 40 k-points outside the data folder; differences 0.0.
- The slice agrees exactly with the first 40 k-points of the full CrSe SOC band file read by the parser.
- Full CrSe `2-SOC/2-bands` (354 MB, 640 k-points, 120 bands, 12 ions): 6.7 to 7.0 s and 291 MB peak memory, about the size of the array.
- Largest file, CrSe `2-SOC/1-scf` (785 MB, 1440 k-points): 20.3 s and 626 MB peak memory for a 597 MB weight array.
- Tests: `tests/test_orbital_parser.py`, 32 tests, on synthetic fixtures in `tests/fixtures/orbital/` and on the MnTe `9.2-band` file.
- The fixtures cover nonmagnetic, ISPIN = 2 (with f orbitals) and SOC layouts, plus missing, short, extra and truncated blocks and a permuted dimension order.
- `mnte_9.2-band_trimmed_vasprun.xml` is cut from the MnTe `9.2-band` file to 2 k-points and 4 bands; a test checks it against the full file.
- The real-file tests write the labels (total, mx, my, mz) and the single eigenvalue block literally, not from the parser's tables.
- Rebuild the fixtures with `python tests/fixtures/orbital_vasprun_builder.py`; a test checks the stored files match the builder.
- Rebuild the trimmed fixture with `python tests/fixtures/orbital_vasprun_builder.py --trim <path to 9.2-band/vasprun.xml>`.
- Full check: `PYTHONPATH=. python -m unittest discover -s tests` runs 108 tests, OK, one skip (76 before).
- Unverified: nonmagnetic and collinear ISPIN = 2 layouts on real files, because no such file is in the data; only synthetic fixtures and pymatgen.
- Unverified: multi-gigabyte files and real .xml.gz files.
- Unverified: the sign and axis convention of mx, my, mz against VASP's own documentation; pymatgen reads them the same way.
- Not done (stage 2 and 3): atom and orbital selection, orbital groups, k-path distances, plots and the full command line.
- Question for Diana: the CrSe "noSOC" runs are noncollinear without SOC; should the plotter treat them like SOC runs (total charge only)?
- Question for Diana: do you have a nonmagnetic or a collinear ISPIN = 2 band run with LORBIT = 11 we may use for checks?
- Question for Diana: the MnTe `9.2-band` folder holds VASPKIT output (BAND.dat, KLABELS); may stage 2 compare its fat bands with them?
- Question for Diana: should the k-point labels come from `<kpoints_labels>` and the line-mode divisions when no labels are given by hand?

## Status for Diana, 26 October 2026

- Main now holds your MAGMOM fix (pull request 1), the rescan and delete actions, and the spin projection step. Main is 1e2eb10.
- Main's history was rewritten twice on 26 October to replace personal email addresses by GitHub no-reply addresses. Please re-clone, or run `git fetch origin` and `git reset --hard origin/main` on a clean checkout.
- Your commit kept its content, message and date; only the author line changed. Please set `git config --global user.email` to your GitHub no-reply address so future commits need no rewrite.
- Rule for everyone: no personal or work email address in any commit, file or message of this repository. Merges into main are done locally by Josh after the tests, never through the GitHub website.
- Request 1, rescan and delete: done, see the section "Rescan and delete, 26 October 2026". Run `altermag rescan CrSe` and `altermag rescan 13-CrSb` once to rebuild the old summary order.
- Request 2, spin projection: done, see the section "Spin projection, 26 October 2026". On your coplanar no-SOC CrSe Wannier run it matches your script's sigma_z to 1.3e-10 on the same k-points.
- Request 3, the orbital-projected band plotter from vasprun.xml: not started. It will be a separate package on its own branch, built in the three stages of your specification, after your go on the scope.
- Please check: the spin order of the raw VASP Hamiltonians is assumed udud (up, down per orbital); your symmetrised file says "permuted uudd->udud" in its header, so it is read as udud too. Tell us if either is wrong.
- Please check: in the index on Josh's laptop the no-SOC Wannier run of CrSe_real had no parent band run, so the spin button stayed hidden until the parent was set in the interface. Say if the importer should have found it.
- Working convention from now on: one branch per task, a dated section at the top of this file, a pull request to main. Josh's session reviews, tests and merges, then writes its own section here.
- Tests on main: `python -m unittest discover -s tests` runs 76 tests, OK, one skip.

## Spin projection, 26 October 2026

- Added `altermag spin-projection`, based on Diana's reconstruct_spin_sz_bands.py, in `altermagagent/spin_projection.py`.
- It reads a Wannier Hamiltonian, builds H(k) on the k-path of the ancestor band run's KPOINTS and computes the averages of sigma_x, sigma_y and sigma_z per band and k-point.
- The k-path, tick labels, Fermi level and reference choice are the same as in `altermag compare-wannier`, so the two figures line up.
- Both commands share one setup function, `wannier_compare.prepare`; `select_reference(..., latest_if_ambiguous=True)` picks the band run with the highest sequence number when several qualify.
- When the reference has ISPIN = 2 without SOC or LNONCOLLINEAR, its spin-up channel gives the k-path and Fermi level, and a note says so.
- This rule applies to every caller without an explicit spin, so the band match in the suggestions no longer fails on such references.
- A present but unreadable wannier90_wsvec.dat falls back to the direct-HR convention with a note, as before the shared setup.
- Inside energy-degenerate groups (spread up to 1e-5 eV) each component is diagonalised with its own operator, so all three panels are well defined.
- The JSON and the caption say "each component diagonalised inside degenerate groups"; the Pauli operators are applied as row swaps, not dense matrices.
- Spinor check: the run's wannier90.win must not say spinors = false (the Wannier90 default when absent), and num_wann must be even.
- When no wannier90.win is found, the result and the caption carry "wannier90.win not found; spinor basis assumed".
- Spin order rule: udud is up, down per orbital; uudd is all up, then all down.
- The order is read from the first header line ("permuted uudd->udud") or from the file name (_udud, _uudd); `--order` overrides it.
- When nothing states the order, udud is used and marked "assumed" in the caption, the JSON and the printed line.
- Caption: "Colours assume a diagonal spin operator in the Wannier basis, order udud (assumed)"; this is an assumption about the basis, not a measured spin.
- For runs with LSORBIT = .TRUE. the caption adds "with SOC these are not the physical spin"; reading .spn spin matrices is not done.
- Output: a PNG (three stacked panels with `--component all`) and a JSON with the order and its source, Fermi energy and source, window means of |sigma|, min, max, residual and caption.
- Figures go to a new folder per call under `~/.altermagagent/figures/<project>`, never into a scanned folder, so two clicks never share a file.
- Run: `altermag spin-projection <project> --run-name <wannier run> [--emin -1] [--emax 1] [--component x|y|z|all] [--order udud|uudd]`.
- Suggestions: step (e) shows "Next: spin projection of <run>/wannier90_hr.dat ..." only after the band comparison has read the Hamiltonian.
- Interface: a "Spin projection" button next to Compare in the band match block opens the figure with the caption and the numbers.
- Validation on CrSe coplanar no-SOC (`coplanar/1-noSOC/04_wannier`, 96 Wannier functions), in a temporary index with parent links set by hand.
- On Wannier90's own 573 band.kpt points, sigma_z agrees with Diana's script to 1.3e-10 at most, and energies to 5e-11 eV, with Fermi 6.3344 eV and udud.
- On the VASP KPOINTS path (640 points) the window -1 to 1 eV holds 8525 states; mean |sigma_z| 0.167, |sigma_x| 0.039, |sigma_y| 0.048.
- With the per-component rotation (files _v3) sigma_z is unchanged; mean |sigma_x| rose by 1.0e-4 and mean |sigma_y| by 1.3e-4.
- The two paths share their corners but not their points, so the per-point comparison was made on Wannier90's points.
- Without a wsvec file the direct-HR convention is used, as in Diana's script, and the JSON says so.
- Re-run after the clean-up gave the same window means, states and 640 points; the files carry a _v2 suffix in the validation folder.
- Tests are in `tests/test_spin_projection.py` and `tests/test_wannier_compare.py`; shared setup is in `tests/_fixtures.py`.
- The test home also redirects the default figure folder, so tests never write under `~/.altermagagent`.
- Full check ran 76 tests, OK, one skip; sigma_z still agrees with Diana's script to 1.3e-10 on her 573 points.
- Open: in the real index the no-SOC Wannier run has no parent link, so its band match has no reference and the suggestion stays hidden until a parent is set.

## Rescan and delete, 26 October 2026

- Added `altermag rescan <project>`: it re-imports the project from its stored root folder under the same name.
- Rescan rebuilds the runs and the summary with the current code, which fixes the stale summary order found on 25 October.
- To fix CrSe and 13-CrSb, run `altermag rescan CrSe` and `altermag rescan 13-CrSb`, or use Rescan in the interface.
- When the stored root folder no longer exists or the summary cannot be built, rescan stops with a clear message; a failed rescan leaves the index as it was.
- User-set parent choices and cached band matches survive a rescan, as they already did on re-import.
- Added `altermag delete <project>`: it removes the project, its runs, band matches and parent choices from the index only.
- Delete never touches folders on disk; cached figures under the index folder are also left in place.
- Interface: each project in the left column has a small menu button (three dots) with Rescan and Delete.
- Delete asks for confirmation; when the selected project is removed, the first remaining project is shown.
- Every rescan and delete message says that no data folder was changed.
- Added five tests in `tests/test_rescan_delete.py` on the demo tree and a temporary index.
- Full check: `python -m unittest discover -s tests` ran 57 tests, OK, one skip.
- Both menu actions were clicked in a headless browser against a temporary index; scanned calculation folders were not changed.

## Structure moment parsing, 25 October 2026

- Fixed `altermagagent/structure.py` to ignore surplus MAGMOM values after all atoms have moments.
- Noncollinear input consumes three numbers per atom; collinear input consumes one.
- Too few values still produce an explanatory note and no arrows.
- Added seven regression tests for surplus values, repeats, continued lines, SAXIS rotation and saved figures with supercell arrows.
- Full check: `.venv/bin/python -m unittest discover -s tests` ran 53 tests with no failures and eight skips.
- Existing structure figures remain cached in the interface.
- Regenerate an affected figure with `altermag structure <project> <run_name>` and refresh the browser after restarting the server.
- Scanned calculation folders and the project index were not changed.

## GUI summary order investigation, 25 October 2026

- The server on port 8080 uses this repository through its editable `.venv` installation.
- The current summary builder puts Structures investigated first.
- The GUI renders `projects.summary_json` in its saved order without rebuilding it.
- CrSe and 13-CrSb retain the older order: Magnetic, SOC, Hamiltonian, Structures, Symmetries, Transport, Status, Last update.
- MnTe-test already stores the current order: Structures, Magnetic, Symmetries, SOC, Hamiltonian, Transport, Status, Last update.
- Local HTTP responses from the same server confirmed the different orders for CrSe and MnTe-test.
- Re-importing each affected project under its existing name regenerates its summary with the current code.
- Refreshing the browser or restarting the server alone does not regenerate stored summaries.
- This was an investigation only; application code, the index and scanned folders were not changed.
- Next: regenerate affected summaries or add an explicit summary cache migration when updating the application. Add option of rescaning the project. Add option of deleting the object from menu.

## Where things stand, 21 October 16:45 EDT (read this first)

- Repository main has everything built today: importer, storyline with evidence tags, diagnosis block, benchmark with negative fixture, band and structure figures, chat with hosted model or Ollama, README, video script, write-up draft, public post text.
- Real data: Diana's MnTe project, four runs (9-scf, 9.2-band, 9.4-band, 9.6-wannier), imported on Josh's laptop. Diagnosis: compensated order supported; non-SOC spin splitting not supported; next calculation named.
- Submission: through Moment only. Josh signed up and created the team. Diana and Poulami: sign up at https://scientificintelligence.org/onboarding?campaign=osi-hackathon-2026 (Globus or Auth0), then join the team with this link: https://scientificintelligence.org/teams/join/r-my-ovgYj0Z5qoJQ8mEIeIs-YdA6CZ-LxH-CxsKNsk Deadline Thursday 22 October 7 PM Central = 8 PM EDT = 3 AM Israel = 5:30 AM India. Internal target 6 PM EDT.
- Deliverables: video up to 4 minutes (docs/video_script.md), public post (docs/public_post.md), one-page write-up with figure (docs/writeup.md, docs/figures/writeup_figure.png), repository link (public on Josh's word before submission).
- Scope is frozen (plan v14 on Josh's laptop). No new features. Tomorrow is corrections, write-up edits, recording, submission.

## Overnight note, 22 October 01:15 EDT (read after the section above)

Three branches were built overnight. None is on main. Merge order for Josh at noon, each only after its check:
- accuracy (commit 8d69b09): importer fixes for parent links (no volume rescale, Hubbard U must match, SOC is a hard requirement) plus tests/accuracy.py on 20 synthetic trees: phase 100 %, status 100 %, parent 98.5 %. Check: Diana's MnTe storyline unchanged after import. Merge first.
- literature (commit 63059c2): section "From the literature, not from these files" under the diagnosis, through Valency with VALENCY_TOKEN; without the token it prints "literature unavailable" and carries on. MnTe shows six papers (Osumi 2023, Lee 2023, Lovesey 2023, and others). Check: run summary with and without the token. Merge second.
- symmetry (commit 23f49fb): magnetic space group, spin group, wannhr_symm operations, all tagged "parsed, pending validation". Check: Diana writes confirm or wrong under "Symmetry validation by Diana". Merge only on her confirm.
Merging: git checkout main && git merge --no-ff accuracy, then run .venv/bin/python -m unittest discover -s tests, then push; repeat per branch. Expect small conflicts in summary.py and HANDOFF.md between literature and symmetry; keep both sections.

Decisions recorded tonight: public post on LinkedIn with the video hosted on YouTube; Valency and marimo are real uses, Adaption is not used.

## Integration candidate, 22 October 08:15 EDT (branch integration, not on main; Josh decides at noon)

- Branch integration = origin/graph + origin/accuracy + origin/literature + origin/symmetry, merged in that order, conflicts resolved by hand.
- Conflicts resolved: importer.py (accuracy and graph parent rules), summary.py (twice), cli.py (twice), app.py, README.md, HANDOFF.md.
- Parent links keep both behaviours: restart links need the same SOC and the same U, no volume rescale (angle tolerance 1 degree).
- When the restart file is gone, the only compatible earlier run is still linked with "restart file not present now", or candidates are listed.
- Runs without a restart tag keep graph's rule: same structure and the same magnetic realization, SOC as a tie-breaker.
- Accuracy (tests/accuracy.py, 20 trees): phase 137/137 100.0 %, status 137/137 100.0 %, parent 135/137 98.5 %.
- Before the scoring fix parent was 130/137 (94.9 %): five "unknown" truths where the restart file was deleted got a "restart file not present now" link.
- In all five the importer named the true parent; the scoring now accepts such a link or candidate list when it holds the true parent (docs/accuracy.md).
- Literature: "From the literature, not from these files" appears once, under the suggestions of the selected branch, in the interface.
- It also follows the Branches list in `altermag summary` and the suggestions in `altermag suggest`; the chat tool search_literature stays.
- Queries now come from the branch's open steps instead of the removed diagnosis; without VALENCY_TOKEN it prints "literature unavailable" and continues.
- Symmetry: step (b) "define MSG" runs spglib on the branch's first converged run (OUTCAR moments, else INCAR MAGMOM) and still writes findsym.in.
- The computed label is tagged "parsed, pending validation"; with a findsym.log the label shows "from findsym.log" and the computed one with agree or differ.
- Interface: a "define MSG" button under the suggestions recomputes and lists the result (checked in headless Chromium on CrSe_real).
- Symmetries point per branch: findsym label if present, else the computed label tagged pending, else "define MSG".
- CrSe_real coplanar-noSOC: findsym.log P-62'm' (BNS 189.225); spglib on 022_scf OUTCAR gives P-62'm' (BNS 189.225), 12 operations: agree.
- CrSe_real coplanar-SOC: computed P-62'm' (BNS 189.225), pending (1-scf OUTCAR).
- CrSe_real collinear-noSOC: computed P6_3'/m'm'c (BNS 194.268), pending (02-noSOC-bands INCAR MAGMOM, no OUTCAR).
- CrSe_real collinear-SOC: computed P6_3'/m'm'c (BNS 194.268), pending (03-SOC-scf INCAR MAGMOM, no OUTCAR).
- MnTe realization-1-SOC: a findsym.log is now in 9-scf and says Cm'c'm (BNS 63.462); spglib on 9-scf OUTCAR gives Cm'c'm (BNS 63.462): agree.
- Demo: realization-1-noSOC P6_3'/m'mc' (BNS 194.269), realization-1-SOC P6_3'/m'm'c (BNS 194.268), both pending.
- Tests: `python -m unittest discover -s tests` runs 30 tests, all pass (one new test for define MSG in tests/test_negative_fixture.py).
- Commands checked on CrSe_real, MnTe and CrSe_demo: import, summary, graph, suggest (every branch), with and without VALENCY_TOKEN; all exit 0.
- `altermag graph` output is identical to branch graph for all three projects, so the merged parent rules change no real-data link.
- Screenshots: ~/development/osi-hackathon/screenshot_integration_CrSe_v1.png and screenshot_integration_MnTe_v1.png.
- They show the branch graph, the suggestions, the define MSG button, the literature section with papers and the symmetry labels.
- Open: Diana confirms the computed labels for coplanar-SOC and the two collinear branches by running findsym on the written findsym.in files.

### Diana, morning order
1. Moment sign-up and join the team.
2. Accept the GitHub invite, clone, install, run the tool on your MnTe folder (altermag import, summary, bands MnTe 9.4-band, gui).
3. Corrections from Diana: every wrong value with the file where the truth is.
4. Fill docs/benchmark_mnte_template.json with your expected values, then run altermag benchmark MnTe docs/benchmark_mnte_template.json.
5. Edit the physics in docs/writeup.md. Read docs/video_script.md and adjust the narration.
6. Answer: is "number of nodes equal to NBANDS" the BSUB slot count? May the trimmed MnTe set be published as the sample?
7. Push on branch diana. Post in Slack.
- Estimate how long the same inventory of your MnTe folder takes you by hand, for the write-up's time comparison.


### Diana, extra step 8: validate the symmetry branch (not on main until you confirm)
- Run: git fetch origin && git checkout symmetry, then altermag import and altermag summary on your MnTe folder.
- The Symmetries point shows a magnetic space group label from the converged OUTCAR moments, tagged "parsed, pending validation".
- MnTe result to check: Cm'c'm (BNS 63.462), type 3, 8 operations, 4 unitary and 4 antiunitary; wannhr_symm mysymmop.dat has 24 operations and its magmom logs keep 4 without time reversal and 0 with, while spglib keeps 4 and 4.
- Write under "Symmetry validation by Diana" in HANDOFF.md: confirm or wrong for the label, the counts, the unitary/antiunitary split, and whether the collinear demo label should be shown. See docs/symmetry_notes.md on that branch.
- Then git checkout diana before you write your other notes.


## Install test by Poulami (posted in Slack 22 October, 20:23 India, 10:53 EDT)
- Computer: Windows (HP laptop), Python 3.14. Proxy workaround that worked: phone hotspot.
- pip install: ok, about 20 minutes on the hotspot. summary: ok. band figure: ok. browser page with graph: ok. figure window and file window: ok. chat: not tried. Errors: none.
- Screenshots in Slack: the summary and branch graph, the run window of bands_nosoc with the 3D structure view, the File window of relax/OUTCAR, the band figure and three structure figures.
- Real-data reproduction by Poulami (22 October, 12:00 EDT): she downloaded Diana's CrSe archive, extracted the coplanar/2-SOC band and scf runs on her Windows laptop, imported them and drew the SOC band figure with the tool. Her figure matches docs/figures/real/CrSe_coplanar_SOC_bands.png: same path labels −M, Γ, M, K, Γ, A, L, H, K, same bands. Second machine, real data, same result.


### Poulami, morning order
1. Moment sign-up and join the team.
2. Install test from your document (poulami-tasks_v2), chat test with Ollama, report block here under "Install test by Poulami", push on branch poulami.
3. If you have old VASP folders: import one and list wrong values under "Corrections from Poulami".

### Rule if someone is unavailable
If no note from a person by 02:00 EDT Thursday, Josh's session takes their tasks.


Read this file first and write it last in every session.

## Done
- Done (main, 22 October, Diana's corrections): summary points now Structures, Magnetic, Symmetries, SOC, Hamiltonian, Transport, Status, Last update.
- Symmetries point shows only the magnetic space group per magnetic realization with its magnetic point group, no BNS numbers or bracket tags; crystal space group stays in Structures.
- Step (b) is done with the label when a magnetic space group exists (findsym.log or computed); the define MSG action shows only without a label.
- Run cards show number, name, date, badges, phase, short convergence for scf runs and the parent line; INCAR tags and changes versus the parent moved into the run dialog.
- 04_wannier forced against 030_bands (test only): no band mapping, 030_bands has 87 bands for 96 Wannier functions; against 031_bands_nbands120: no Fermi energy without a common SCF parent.
- Done (main, 22 October): Wannier band comparison module by Diana Golovanova (altermagagent/wannier_compare.py), `altermag compare-wannier` (alias `wannier-compare`), step (e) band match line and a Compare pill on Wannier runs, cached in the index; docs/wannier_comparison/README.md.
- Results (window -1 to +1 eV, preliminary because INCAR and POSCAR hashes differ): CrSe_real coplanar-SOC raw hr versus 2-bands MAE 0.040 meV, RMSE 0.074 meV, max 0.656 meV.
- The symmetrised CrSe hr (wannier_sym/wannier90_hr.dat_primed_nsymm_6) gives MAE 452.9 meV, RMSE 590.0 meV, max 1686.8 meV, matching the rejection in Diana's 3.1-wannier README.
- MnTe 9.6-wannier raw hr versus 9.4-band: MAE 4.058 meV, RMSE 7.217 meV, max 34.324 meV; versus 9.2-band: MAE 5.726 meV, RMSE 8.985 meV, max 34.324 meV.
- CrSe_real 04_wannier: "no ancestor band run with eigenvalues" because neither it nor 030_bands has an indexed parent; set parents with `altermag parent` to compare.
- Done (main, 22 October): silent screen recording of the CrSe_real walkthrough, 3:58, 12 scenes with caption bands, video/altermagagent_silent_v1.mp4 and video/scenes.json on Josh's laptop (not in git); method in docs/recording.md.
- Done (integration, 22 October): physicist validations in docs/validations.json (six entries from Diana); summary.py reads it (ALTERMAG_VALIDATIONS, else docs/validations.json, plus ALTERMAG_HOME/validations.json) and upgrades matching symmetry and convergence tags; re-import to refresh the Symmetries point.
- Done (integration, 22 October): convergence of every scf run is shown on the run cards, in the run dialog, in `altermag summary` and in the chat data, with a red "unconverged" badge and "N unconverged scf" in the Status point. CrSe_real: 020_scf and 021_scf NOT reached, 022_scf and 023_scf_nbands120 reached; MnTe 9-scf reached. Screenshot: screenshot_integration_CrSe_v2.png on Josh's laptop.
- Done (integration, 22 October): video script v3 on CrSe_real (3:15 planned), write-up and public post on branches and six steps, figure docs/figures/writeup_figure_v2.png. Open: the 020/021 non-convergence is in the index only, not shown in the interface.
- Done (branch graph, 22 October morning, not on main): Diana's CrSe test requests 1 to 9.
  Runs are grouped into branches by structure (StructureMatcher, no volume rescale, 1 % lattice tolerance), magnetic realization and SOC.
  The realization comes from INCAR MAGMOM (or converged OUTCAR moments): collinear, coplanar, non-coplanar or non-magnetic, pattern up to a global rotation within 0.1 μB.
  CrSe_real: coplanar-noSOC (7 runs), coplanar-SOC (3), collinear-noSOC (1), collinear-SOC (2); MnTe: realization-1-SOC (4); demo: realization-1-noSOC (3), realization-1-SOC (4).
  The summary has eight points: the spin splitting point is gone (the scan stays in `altermag bands`); SOC, Hamiltonian and Symmetries are one short line.
  The diagnosis block is replaced by per-branch "Suggestions for the next DFT calculations": steps (a) to (f) with done, partial, missing or not applicable, and the next calculation.
  Interface: a git-like SVG graph of the branches above the runs; clicking a node or lane opens the branch header (3D view), its run cards and its suggestions.
  New commands: `altermag graph <project>` and `altermag suggest <project> <branch>`; missing MSG writes `~/.altermagagent/msg/<project>/<branch>/findsym.in`.
  The parsers branch is merged into graph and wired: MSG label from findsym.log, wannier90.win key info, INCAR changes versus the parent on run cards.
  Parent links: when ICHARG=11 or ISTART>0 and the restart file is gone, the only compatible earlier run is linked with "restart file not present now".
  Tests: 7 negative fixture tests (step (d) missing without bands_nosoc) and the parser tests pass.
- Done (branch accuracy): measured accuracy on 20 synthetic trees, phase 100 %, status 100 %, parent 98.5 % (was 78.8 %), see docs/accuracy.md.
- Fix (accuracy): parent matching no longer rescales volume (ltol 1 %), so a strained cell is another structure; the Structures point is unchanged.
- Fix (accuracy): ICHARG=11 and ISTART links need the same LDAU, LDAUTYPE, LDAUL, LDAUU and LDAUJ; the note says "same U".
- Fix (accuracy): SOC is a hard requirement for ICHARG=11 and ISTART links; without a usable file the note says "parent unknown (restart file missing)".
- Known limit (accuracy): a file copied with cp -p keeps the parent's time and can reorder runs; left as is.
- Done (Josh's session, branch literature, not merged): a literature line kept apart from the file evidence.
- New `altermagagent/literature.py`: a small MCP client for Valency (urllib only, 10 s timeout) and `search_literature(query, limit=3)`.
- It returns title, first author, year, source, url and citation count, or an empty list with "literature unavailable: <reason>"; it never raises.
- `altermag summary` prints "From the literature, not from these files" after the diagnosis; `altermag literature <project>` prints only that section.
- The interface shows the same section under the diagnosis with links, filled in a worker thread; the chat has a `search_literature` tool.
- The chat prompt says literature answers are labelled "from the literature", cite title and year, and are never merged with file evidence.
- Results are cached per search for 24 hours in `~/.altermagagent/literature_cache.json`; the token comes from `VALENCY_TOKEN` only.
- Tested on MnTe and CrSe_demo with a real token; without the token the summary runs and prints "literature unavailable: VALENCY_TOKEN is not set".
- A bad token gives "literature unavailable: the Valency token was rejected (HTTP 401)"; `python -m unittest discover -s tests`: 10 tests pass.
- Done: importer wall time is measured and stored in the project row; MnTe (hackaton_exa) scanned 4 runs and 175 files in 0.86 s, CrSe demo 7 runs and 61 files in 0.77 s; `altermag import` and `summary` print it.
- Done: docs/video_script.md is version 2, built around "can these existing calculations support the claim we want to make?"; 2:35 planned, 3:00 with the two optional scenes (symmetry, literature) if merged by noon.

As of 21 October, 14:30 EDT:

- First running version is on `main`.
- Importer, nine summary points, gaps list, run cards with file links, three-panel interface and command line all work.
- Tested on the synthetic CrSe tree and on Diana's real MnTe folder.
- MnTe result: three SOC runs (9-scf, 9.2-band, 9.6-wannier).
- MnTe result: Mn moments ±2.25, net moment 0.00, space group P6_3/mmc.
- MnTe result: the Wannier Hamiltonian is present and the missing transport step is listed as a gap.
- The interface opens a project directly with `?project=NAME` in the address.

As of 21 October, evening EDT (Josh's session):

- Band figure tool: `altermag bands <project> <run_name> [--emin] [--emax] [--style paper]` in `altermagagent/bands.py`.
- Output: PNG, SVG, interactive HTML and a result JSON in `~/.altermagagent/figures/<project>/`, never in the run folder.
- Interface: the Plot pill and the run window of band runs draw the figure on first open; energy inputs and Redraw are below it.
- Agent: `plot_bands` is a tool the model can call; tested with a stand-in client, not yet with a real API key.
- MnTe 9.2-band: KPOINTS is VASPKIT line mode with labels Γ, Y, C₀, Σ₀, Γ, S; 100 k-points, 60 bands, SOC, one colour.
- CrSe demo bands_nosoc: two colours; largest spin splitting 0.42 eV, band 4, k-point 3.

As of 21 October, late evening EDT (Josh's session):

- Diana's procedure is in `docs/dft_run_procedure.md`; the rules taken from it are in `docs/detection_rules_from_procedure.md`.
- Each run now has a phase (relax, scf, band, wannier, symmetrisation, projection, transport, other), shown on the run card.
- Sub-folders without INCAR (for example `wannier_sym`, `AHC`) are read as extra steps of their run.
- Each run has procedure checks (EDIFF, vasp_ncl for SOC, amn/mmn/eig, num_bands, num_wann, slots); warnings print in `altermag summary`.
- MnTe phases: 9-scf = scf, 9.2-band = band, 9.6-wannier = wannier + symmetrisation, transport.
- MnTe Hamiltonian and Symmetries points: symmetrised Hamiltonian present, wannhr_symm logs found, mysymmop.dat lists 24 operations.
- MnTe Transport point: calculated, 9.6-wannier/AHC/ahc.dat, setting task fullahc; the transport gap is gone.
- MnTe checks: EDIFF is 1E-5 in all three runs; 9.6-wannier run.wexac.sh has 40 slots, NBANDS is 60.

As of 21 October, night EDT (Josh's session):

- Done: marimo notebook `notebooks/altermagagent_walkthrough.py` with summary, diagnosis, runs, band sliders, structure and benchmark.
- Tested with `marimo run` and headless Chromium on MnTe; screenshot `docs/figures/marimo_walkthrough_v1.png`.
- The chat system prompt has a "Workflow knowledge" section built from the rules document.

As of 21 October, afternoon EDT (Josh's session):

- The chat works without an API key through a local Ollama server, found at `OLLAMA_HOST` (default `http://localhost:11434`).
- Order: `ANTHROPIC_API_KEY` set uses Anthropic unchanged; else Ollama if it answers; else the notice, which now points to the README.
- Ollama model comes from `ALTERMAG_OLLAMA_MODEL` (default `llama3.1:8b`); plot_bands works as a tool, at most 3 rounds.
- A model without tool support gets plain chat, and the answer says to draw figures with the Redraw button.
- Every call is logged with provider, model and token counts in `~/.altermagagent/agent_calls.log`.
- Tested on MnTe with llama3.2:3b on Josh's laptop: both questions answered, and the band figure for -1 to 1 eV was drawn.
- The 3b model also called plot_bands for a plain storyline question; a larger model (llama3.1:8b or qwen) should do better.

As of 21 October, night EDT (Josh's session):

- Structure figure: `altermag structure <project> <run_name> [--supercell 1 1 1]` in `altermagagent/structure.py`.
- It reads POSCAR (or CONTCAR) and INCAR MAGMOM through pymatgen and draws atoms, cell edges and moment arrows.
- Output: PNG, SVG, JSON and a 3Dmol.js HTML view in `~/.altermagagent/figures/<project>/<run>_structure.*`; the HTML needs internet.
- Interface: a Structure pill on every run with POSCAR or CONTCAR; the scf file pill is now called "Structure file".
- MnTe 9-scf (SOC): Mn1 (3.98, 2.30, 0), Mn2 (-3.98, -2.30, 0), Te 0, net 0; arrows opposite and in-plane.
- CrSe scf_nosoc (collinear): Cr1 +3, Cr2 -3 along z, Se 0, net 0; arrows opposite.
- The arrows are the INCAR starting moments, not the OUTCAR moments; the notes say so.

As of 21 October, night EDT (Josh's session, install test and docs):

- Fresh-machine install test passed: clean clone, new python3.14 venv, README install and demo commands, separate index through `ALTERMAG_HOME`.
- No warnings or errors; no packaging fix was needed.
- Pip install timing: about 3 minutes with an empty pip cache, 37 seconds with a warm cache; `pip install -e .` took 2 seconds.
- Python 3.12 check with `pip download` (no install): every wheel exists for Linux x86_64 and Mac arm64.
- On Windows the seven pinned packages have 3.12 wheels; uvloop is skipped there by its own platform marker.
- README rewritten for newcomers: install per system, five-command demo, commands, `?project=`, chat, index and figures, known limits, screenshots.
- Video script and recording checklist in `docs/video_script.md`.

As of 21 October, night EDT (Josh's session, outside review adopted):

- Every summary point and every check has an evidence state: output found, parsed, validated or insufficient evidence.
- The state is a small tag in the interface and in brackets in `altermag summary`, with a note saying what was checked.
- Magnetic is validated only when the net moment is below the declared tolerance 0.01 μB per cell; the tolerance is printed.
- Transport on MnTe: ahc.dat output found and parsed (1000 rows); not validated, no convergence or unit check exists.
- A diagnosis block replaces any material label: "Supported by these files", "Not supported by these files", "Next calculation".
- The block ends with "Literature classification is not read from these files."; it never says a material is an altermagnet.
- MnTe: supported = compensated collinear Mn ±2.25 order (9-scf/OUTCAR) and hexagonal P6_3/mmc (9-scf/CONTCAR).
- MnTe: not supported = non-relativistic spin splitting, because the folder has no collinear non-SOC run.
- MnTe next: collinear non-SOC SCF (ENCUT 450, MAGMOM 4.60 -4.60 0 0), then a band run on the 9.2-band k-path.
- Splitting is worded as the largest up-down difference among the sampled k-points and bands in the window, with the k-point count.
- K-point indices are now 1-based in both the summary and the band figure tool.
- Parents: a run can be "parent unknown" or "candidates: a, b" with the reason; the SOC setting breaks ties.
- A "change parent" selector on each run card, and `altermag parent <project> <run> <name|unknown|auto>`, store a choice in the index only.
- The choices live in table parent_choices of `~/.altermagagent/index.sqlite` and survive a re-import.
- Procedure checks are labelled "local procedure check (Diana's notes), pending her confirmation", not physics errors.
- The structure dialog caption says "Input magnetic configuration from INCAR MAGMOM".
- Benchmark: `altermag benchmark <project> <file.json>`, `docs/benchmark.md`, `demo/benchmark_crse.json`, `docs/benchmark_mnte_template.json`.
- Negative fixture: `python demo/make_demo.py --without-nosoc-bands`; diagnosis says not supported, splitting row prints UNKNOWN.
- Results: CrSe demo 7 pass; negative fixture 6 pass, 1 unknown; MnTe template 8 unknown until Diana fills it.
- `python -m unittest discover -s tests`: 4 tests pass.
- Headless browser check on MnTe, CrSe_demo and CrSe_negative: tags, diagnosis block, parent selector and caption all shown.

As of 21 October, late night EDT (Josh's session, write-up):

- Write-up draft (one page, 582 words) in `docs/writeup.md`, with the figure `docs/figures/writeup_figure.png`.
- The figure puts the MnTe 9.4-band band structure beside the diagnosis block, rendered as text from `altermag summary MnTe`.
- Public post draft (115 words) in `docs/public_post.md`; repository and video links are placeholders.

As of 21 October, near midnight EDT (Josh's session, branch symmetry only):

- New `altermagagent/symmetry.py`: magnetic space group (spglib), spin group (spinspg, collinear non-SOC only), wannhr_symm operations.
- Every label is tagged "parsed, pending validation"; the diagnosis adds one Supported line with the magnetic group label only.
- MnTe: Cm'c'm (BNS 63.462), type 3, 8 operations, 4 unitary, 4 antiunitary; wannhr_symm 24 operations, differs (logs keep 0 antiunitary).
- CrSe demo: P6_3'/m'mc' (BNS 194.269), type 3, 24 operations, 12 and 12; spin group G 6/mmm, H -3m.
- spinspg 0.3.2 added to requirements.txt; run `.venv/bin/pip install -r requirements.txt` after checkout, then re-import.
- Tests: 9 pass. Headless browser check on MnTe and CrSe_demo: the Symmetries point and the diagnosis line render.

## Next, in order

00. Josh: review branch graph (merges origin/parsers), run the tests, then merge into main; it conflicts with the literature and symmetry branches in summary.py (the diagnosis block is gone).
00a. Diana: check the CrSe branch names and the six suggestion steps per branch, and run findsym on the written findsym.in of coplanar-SOC.

0. Diana: fill the expected values in `docs/benchmark_mnte_template.json` and confirm the procedure checks wording.
0b. Diana: edit the physics in docs/writeup.md on branch diana.
0a. Diana: run the collinear non-SOC SCF and band run listed under "Next calculation" for MnTe, if she agrees.

1. Diana: confirm that "#BSUB -n" is what she means by "number of nodes" and that the AHC folder counts as the transport step.
2. Chat test on real data once `ANTHROPIC_API_KEY` is set, including the paper-style band figure request.
2a. Chat test with llama3.1:8b through Ollama, to check answer quality against the 3b test model.
3. Diana: compare the 9.2-band figure in `~/.altermagagent/figures/MnTe/` with BAND.jpg in the run folder.
4. Agent tool `plot_structure` in agent.py (TOOLS entry plus a branch in run_tool), and arrows from OUTCAR moments as an option.
5. Diana: read `docs/video_script.md` and change any wording; Josh records the screen with the checklist in it.
6. Day two items after that.

## Blocked / waiting

- `ANTHROPIC_API_KEY` is not set yet.
- GitHub user names of Diana and Poulami are needed for the collaborator invites.
- Diana: may the trimmed MnTe set be published?
- Answered by Diana, 21 Oct: the space group and moments are correct. The files in `9.6-wannier/wannier_sym` come from wannhr_symm, the Wannier90 symmetrisation tool, which also writes symmetry analysis files. Day two: parse its logs into the Symmetries point.
- Poulami: GitHub name, Python version, and whether she has VASP.
- This year's deadline time.

## Diana, your morning (from 06:00 Israel)

1. Pull `main`.
2. Run the tool on your full MnTe folder and on any other project folder you have.
3. For every wrong value, write "point, shown value, correct value, file where the truth is" below under "Corrections from Diana".
4. Answer here: are the Mn moments of 2.25 and the P6_3/mmc label right?
5. Answer here: does 9.2-band have a labelled line-mode KPOINTS file?
6. Write the physics wording you want for the nine points under "Wording from Diana".
7. Write the two-minute narration under "Narration".

### Symmetry validation by Diana

Confirmed by Diana on 22 October, recorded in docs/validations.json (status "validated by physicist"):
- CrSe_real collinear-noSOC and collinear-SOC: magnetic space group P6_3'/m'm'c (BNS 194.268) is correct.
- CrSe_real coplanar-noSOC: P-62'm' (BNS 189.225) agrees with her findsym.log.
- MnTe realization-1-SOC: Cm'c'm (BNS 63.462) agrees with her own 9-scf/findsym.log.
- CrSe_real coplanar/1-noSOC/020_scf and 021_scf are known unconverged trials.
- CrSe_real: the three collinear runs are a started, uncompleted project.
- CrSe_real: the umbrella folder has no calculations.
- `altermag summary` now tags these labels and the 020/021 convergence lines "validated by physicist (Diana, 22 Oct)" when the computed value matches, else "differs from validated value".
- Open: coplanar-SOC P-62'm' is still "parsed, pending validation".

Branch `symmetry`, not merged into main. Labels and questions are in `docs/symmetry_notes.md`.
Write "confirm" or "wrong, correct value, file" after each line.

- MnTe magnetic space group Cm'c'm (BNS 63.462, type 3), moments from 9-scf/OUTCAR:
- MnTe operation counts: 24 in mysymmop.dat (P6_3/mmc), 8 in the magnetic group:
- MnTe unitary and antiunitary split: spglib 4 and 4, wannhr_symm logs 4 and 0:
- CrSe demo spin group: G 6/mmm, H -3m, collinear, no translation or inversion swaps the spins:

### Corrections from Diana

### Wording from Diana

### Narration

## Poulami, your morning (from 09:00 India)

Note for Poulami: chat works without a key if Ollama is installed, see README.

1. Task B: install test following the brief in Slack; paste the report block below under "Install test by Poulami".
2. Task C: if you have VASP, run the MnTe sample on your cluster and post the Drive link here.
3. Task D: run the tool on your own DFT folders and list wrong values below under "Corrections from Poulami".
4. Task E: write README notes below under "README notes".

### Install test by Poulami

### Corrections from Poulami

### README notes

## Rule if someone is unavailable

If there is no note from a person by 02:00 EDT Thursday, Josh's session takes over their tasks.

## Report format for Slack

```
What works: ...
What failed: ...
What is needed: ...
```
