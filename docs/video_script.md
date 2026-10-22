# Video script, version 3

One question runs through the whole video: "Can these existing calculations support the claim we want to make, and what is the next calculation?"
The files never answer "is it an altermagnet"; they answer what is done, what is missing, and what to run next.

The recording uses Diana's real CrSe project (13 runs, four branches); MnTe gets one line only.
Planned length: 3 minutes 25 seconds; target 3 minutes, hard limit 4 minutes.
Diana reads the narration; each paragraph is under 40 words and starts at the time shown.

## Narration

**0:00 to 0:15. The problem and the question**

My CrSe project has 13 VASP runs from two months of work, and nowhere is its history written down.
One question matters: can these existing calculations support the claim we want to make, and what is the next calculation?

**0:15 to 0:30. Import, read-only**

I point AltermagAgent at the folder.
It only reads; it never writes inside it.
It scanned 13 runs and 1127 files in under four seconds.
My MnTe folder, with four runs, is scanned in seconds.

**0:30 to 0:45. The summary**

Eight points sum up the project, each with its evidence file and a tag.
Magnetic lists two realizations: coplanar in ten runs, net moment zero in OUTCAR, and collinear in three runs, from INCAR only.

**0:45 to 1:00. The branch graph**

The runs are grouped into branches and drawn like a git history.
There are four branches, split by magnetic realization and by spin-orbit coupling, or SOC.
I click the coplanar SOC node.

**1:00 to 1:15. The branch's runs**

The branch lists its three runs in order, each with INCAR shown as changes versus its parent.
The first SCF reads a WAVECAR, but no earlier run has the same SOC setting, so its parent is unknown.

**1:15 to 1:40. The suggestions in six steps**

Below the runs, the suggestions follow my six working steps, from structure to transport, each marked done, partial or missing.
On the coplanar SOC branch all six are done: symmetry gives P-62'm', and step e shows the Wannier band match.

On the collinear non-SOC branch, the next calculation is a new SCF, with the exact INCAR tags to set.

**1:40 to 2:00. Symmetries and convergence**

The Symmetries row gives one magnetic space group per magnetic realization, with its magnetic point group.
On the coplanar non-SOC branch, step b reads P-62'm' from my own findsym.log, and the computed label agrees.

The tool also reads every OSZICAR, and the convergence is visible on the run cards.
Runs 020 and 021 carry a red "unconverged" badge: they stopped above EDIFF; step c cites run 023, which converged.

**2:00 to 2:20. Band figure and structure**

The SOC band run is redrawn in a window from minus one to one electronvolt, with the path labels from my KPOINTS file.
The structure view shows the input moments, and I can rotate it.

**2:20 to 2:35. Local procedure checks and the Wannier card**

Each run is compared with my own procedure notes, such as job slots against NBANDS.
These are local checks, pending my confirmation, not physics errors.
The Wannier card lists the energy windows and projections.

**2:35 to 2:45. Wannier Hamiltonian versus VASP bands: mean error 0.04 meV**

My comparison script overlays the raw Hamiltonian of 3.1-wannier on the VASP bands of 2-bands.
Within 1 eV of the Fermi level the mean error is 0.04 meV, the RMSE 0.07 meV and the maximum 0.66 meV.

**2:45 to 3:00. The literature section**

Below the suggestions, a separate heading lists published papers from a literature search.
They come from the literature, not from these files, and they never change the suggestions.

**3:00 to 3:15. The chat**

I ask the chat: "Can you prepare the next step of band structure calculations for the collinear SOC structure, based on the similar calculations for the coplanar SOC branch?"
It prepares the next calculation as a preview, from the coplanar sibling branch: it compares KPOINTS and INCAR, gives the order of steps and what to copy, and never submits.

**3:15 to 3:25. Close**

AltermagAgent is free and open, under Apache 2.0.
On 20 synthetic trees with known answers it scored 100 percent on phase and status, and 98.5 percent on parents.
Run it on your own folders.

## Recording checklist for Josh

### How the recording is made

- The screen recording is made by script on Josh's laptop, browser at 1920 by 1080, zoom 100 percent, browser window only.
- The interface runs at http://localhost:8080/?project=CrSe_real; the terminal shots use a terminal at the same size.
- Diana records the narration as a voice memo; the audio is merged with the screen recording afterwards.
- So each scene must hold on screen at least as long as its paragraph read at a calm pace; when in doubt, hold 2 seconds longer.
- The times above are the plan; after merging, the scene lengths follow Diana's audio, not the other way round.

### Before recording

1. Use one index for every step: `export ALTERMAG_HOME=<recording index>` in each terminal, or leave it unset for `~/.altermagagent`.
2. Run `source ~/.altermagagent/env` in the terminal that starts the interface, so the literature section finds papers.
3. Set the chat key in that terminal: `export ANTHROPIC_API_KEY=sk-ant-...`; without a key use Ollama with `llama3.1:8b`, never the 3b model.
4. Import MnTe once before recording, so it is in the project list.
5. Start `altermag gui`; close other tabs; hide bookmarks and notifications.
6. Open the coplanar SOC branch once, so the literature search is cached and appears without "Searching...".
7. Rehearse the chat question once, so the answer time is known.

### Click sequence per scene

- 0:00 The question: hold on the folder tree of `~/development/osi-hackathon/data/diana_CrSe/CrSe-test` in a file browser (coplanar and collinear folders open).
- 0:15 Import: in the terminal run `altermag import ~/development/osi-hackathon/data/diana_CrSe/CrSe-test --name CrSe_real`.
- 0:15 continued: hold on the line "Scanned 13 runs and 1127 files in ... seconds; the folder was only read." (3.19 s in the last run).
- 0:30 Summary: switch to the browser, reload http://localhost:8080/?project=CrSe_real and hold on the eight points; point the cursor at the bold labels of Structures, Magnetic and Symmetries.
- 0:45 Branch graph: scroll to "Branches" and hold on the four lanes; then click a node in the coplanar-SOC lane.
- 1:00 Runs: scroll to the coplanar-SOC run cards; hold on 2-bands (its INCAR line versus parent 1-scf), then on 1-scf with "parent unknown" and its reason.
- 1:15 Suggestions: on the coplanar-SOC branch, hold on step (b) with its label, the step (e) band match line and the "Recommended next calculation" line.
- 1:40 Symmetries: hold on the Symmetries row of the summary, then click a node in the coplanar-noSOC lane and hold on step (b) "from findsym.log".
- 1:40 continued: scroll to the run cards and hold on the red "unconverged" badges and "NOT converged" lines of 020_scf and 021_scf.
- 2:00 Band figure: click the coplanar-SOC lane again, click the "Plot" pill of 2-bands, set the window to -1 and 1, click Redraw, hold on the figure.
- 2:00 continued: close the band dialog, scroll to the branch header 3D view, drag to rotate it for a few seconds (the 3D view needs internet).
- 2:20 Procedure checks: click the 3.1-wannier run card, scroll to "Local procedure checks (Diana's notes), pending her confirmation", hold.
- 2:20 continued: close the run dialog and hold on the 3.1-wannier card line "wannier90.win: num_wann 96, bands 108, outer window ... projections".
- 2:35 Wannier comparison: on the coplanar-SOC branch, click the Compare pill of 3.1-wannier and hold on the figure with its error box.
- 2:45 Literature: scroll below the suggestions to "From the literature, not from these files" and hold on the note and the papers it found.
- 3:00 Chat: type the next-step question above in the chat, click Send, cut the wait, and hold 18 s on the first screen of the answer with its "Preview only" line.
- 3:15 Close: in the terminal run `python tests/accuracy.py` and hold on the three lines phase 137 / 137, status 137 / 137, parent 135 / 137.
- 3:15 continued: end on a title card with the repository link, "Apache 2.0" and `altermag import <folder> --name <project>`.

### Figure files to show

- Band figure, written when Redraw is clicked: `$ALTERMAG_HOME/figures/CrSe_real/coplanar_2-SOC_2-bands_bands_screen_-1_1.png` (with `.svg`, `.html`, `.json`).
- Structure figure, written when a structure view is drawn: `$ALTERMAG_HOME/figures/CrSe_real/<run>_structure.png`, for example `coplanar_1-noSOC_020_scf_structure.png`.
- Write-up figure, for the title card or thumbnail if needed: `docs/figures/writeup_figure_v2.png` (summary and branch graph of CrSe_real).
- Wannier comparison figure: `$ALTERMAG_HOME/figures/CrSe_real/coplanar_2-SOC_3.1-wannier_wannier_compare_<hash>/overlay.png`.
