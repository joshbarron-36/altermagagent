# Handoff note

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
