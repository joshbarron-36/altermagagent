# Measured accuracy of phase, status and parent detection

## What was generated

`demo/make_variants.py` writes seeded synthetic CrSe project trees, each with a `tree_NN.truth.json` file next to it.
`tests/accuracy.py` generates 20 trees (seeds 1 to 20), imports each into a temporary index and compares the importer with the truth.
Run it with `python tests/accuracy.py`; `python -m unittest discover -s tests` runs it too through `tests/test_accuracy.py`.

Each tree draws its own story from the seed:

- Folder names are picked at random from habits such as `01_relax`, `opt`, `scf`, `static`, `scf_U4`, `U3`, `bands_path`, `soc`, `soc/bands`, `wan90`, `old_try` and `test`.
- One structure (hexagonal CrSe) or two: half of the trees add a strained cell (a +3 %, c -1.5 %) in a sub-folder such as `strain_3pct/` or `a+3%/`.
- Optional Hubbard U runs share the structure of the plain scf run and differ only in `LDAUU` (3, 4 or 5 eV).
- Status: completed, failed with a ZBRENT or EDDDAV message in `vasp.log`, `vasp.out` or `slurm-NNNNN.out`, running, or stopped without any marker.
- Running means files changed in the last 20 minutes and no completion marker; stopped runs have only inputs or a cut-off OUTCAR.
- About 40 % of runs have no `vasprun.xml`, and some completed runs have only `vasprun.xml` and no OUTCAR.
- Parent links: band runs by `ICHARG = 11` with a copied CHGCAR, Wannier runs by `ISTART = 1` with a copied WAVECAR, all others by file time only.
- A relax after a failed `old_try` relax starts from its CONTCAR without a restart tag, as the ZBRENT message suggests.
- Half of the copies keep the parent's file time (`cp -p`); about 15 % of the CHGCAR and WAVECAR files are deleted afterwards to save disk.
- All file times are set with `os.utime` in one consistent story: inputs at the run start, outputs at its end, a few hours to three days apart.
- INCAR spelling varies: `ICHARG=11   # read CHGCAR`, `icharg = 11 ; LORBIT = 11`, `LSORBIT = T`, `IBRION = 2` left in a static run with `NSW = 0`.

The truth for each run is its phase, its status and its parent.
The parent truth is one of three: a run name when the files decide it, "candidates" with the true parent among them, or "unknown".
The truth uses what a careful reader can see in the files: the restart tag, the file the tag needs, the structure, the SOC setting and `LDAUU`.
Without a restart tag the truth is the single earlier run on the same structure, or "candidates" when there are several.
A parent counts as correct when the names agree, when both say "unknown", or when both say "candidates" and the importer list holds the true parent.

## Final numbers

| Field | Correct / total | Percentage |
|---|---|---|
| Phase | 137 / 137 | 100.0 % |
| Status | 137 / 137 | 100.0 % |
| Parent | 135 / 137 | 98.5 % |

The first measurement gave parent 108 / 137 (78.8 %) with 29 mismatches; three importer fixes brought it to 135 / 137.

## Importer fixes made on this branch

- Strained cells: parent matching now uses `StructureMatcher(scale=False, ltol=0.01, angle_tol=1)`, so a 2 % strain is another structure (18 mismatches before).
- Hubbard U: an `ICHARG = 11` or `ISTART` link needs the same `LDAU` flag and the same `LDAUTYPE`, `LDAUL`, `LDAUU` and `LDAUJ` values; the note then says "same U" (7 before).
- SOC: for `ICHARG = 11` and `ISTART` links the SOC setting is a hard requirement, not a tie-breaker (2 before).
- When compatible runs exist but none kept the restart file, the note reads "parent unknown (restart file missing)" instead of naming an incompatible run.

## Which structure tolerance applies where

- Parent links (`assign_parents` in `altermagagent/importer.py`): no volume rescaling, lattice lengths within 1 %, angles within 1 degree.
- The Structures summary point (`altermagagent/summary.py`): unchanged, pymatgen `StructureMatcher` defaults, which rescale volume and allow 20 % length differences.
- So a strained cell counts as the same distinct structure in the summary, but not as a parent source for a run on the unstrained cell.

Runs without a restart tag keep the old rule: same structure, SOC as a tie-breaker, U not compared, because a fresh scf_U run legitimately starts from a relax without U.

The generator was corrected once during this work: later runs on the strained cell now have the first strained run as their true parent.
Before, the truth said "unknown" for them, which contradicted the stated truth rule for runs without a restart tag.

## After the merge with branch graph (integration, 22 October)

Branch graph links a restart run whose CHGCAR or WAVECAR is gone to the only compatible earlier run, with the note "restart file not present now".
The truth says "unknown" for these runs, so the plain rule scored 130 / 137 (94.9 %); in all five extra cases the importer named the true parent.
The scoring now accepts such a link, or such a candidate list, when it holds the true parent; the result is again 135 / 137 (98.5 %).

## Mismatch list

Both remaining mismatches are parent mismatches with the same cause, a run order taken from a copied file's time.
A CHGCAR copied with `cp -p` keeps the parent's time, so the band run sorts before a later run and the importer names one parent where the files allow two.
The named parent is the true one in both cases, but the decision rests on the copy time, not on the restart tag.
This is left as a known limit: the importer takes a run's start as the oldest file time in its folder.

| Tree | Run | Expected | Got |
|---|---|---|---|
| tree_07 | bands | candidates including 02_scf | 02_scf |
| tree_11 | band | candidates including static | static |

## What the synthetic set does not cover

The trees follow one person's tidy story with one material, one code and folder names drawn from a short list, so 100 % on phase and status is an upper bound.
Real folders hold habits the generator does not write: runs restarted inside the same folder and old error logs left behind by earlier attempts.
Files moved between machines with new times, and symbolic links, are not generated either.
Other codes (Quantum ESPRESSO, ABINIT, CP2K), non-VASP Wannier or transport steps, and band paths given as explicit k-point lists are not generated at all.
Copied files that keep their original time (`cp -p`, `rsync -a`) can still reorder runs, as the mismatch list shows.
Status markers are limited to the ones the importer knows, so a real failure with another message would be missed both here and in practice.
