# Instructions for every Claude Code session in this repository

These rules apply to every session: Diana's, Poulami's and Josh's.

## What the project is

AltermagAgent reads trees of VASP calculation folders without ever changing them.
It rebuilds the project storyline and shows it in a three-panel browser interface with a chat about the project.
It is built at the OSI hackathon on 21 and 22 October 2026, as a free and open tool for the community.

## Layout of the code

- `altermagagent/importer.py`: walks a VASP folder tree read-only and writes each run into the index.
- `altermagagent/summary.py`: builds the nine summary points and the gaps list from the indexed runs.
- `altermagagent/files.py`: read-only file helpers for excerpts, folder listings and the file links on each run card.
- `altermagagent/agent.py`: the project chat, which answers only from the indexed summary and runs.
- `altermagagent/app.py`: the three-panel browser interface (projects, storyline, chat), built with NiceGUI.
- `altermagagent/cli.py`: the `altermag` command with `import`, `summary` and `gui`.
- `altermagagent/db.py`: the SQLite index under `~/.altermagagent`, never inside a scanned folder.
- `altermagagent/benchmark.py`: compares a project with a benchmark JSON file, see `docs/benchmark.md`.
- `tests/test_negative_fixture.py`: the negative fixture test, run with `python -m unittest discover -s tests`.
- `demo/make_demo.py`: writes a synthetic CrSe project tree for demos, so no VASP is needed.

## Rules

- Never write anything inside a scanned data folder.
- Never open CHGCAR, WAVECAR or PROCAR; show only their size and date.
- The index lives in `~/.altermagagent/index.sqlite`.
- The agent answers only from indexed data and names the run and file behind every fact.
- No API key in the code; read `ANTHROPIC_API_KEY` from the environment.
- Every function has a one-line header comment directly above it.
- No sentence in any document is longer than 200 characters.
- No attribution lines of any kind in code, documents or commit messages.

## Team rules

- One branch per person, named by first name: `diana`, `poulami`, `josh`.
- Nobody edits another person's files.
- Only Josh's session merges into `main`, and only after a test.
- `HANDOFF.md` is the only place for "done, next, blocked".
- Read `HANDOFF.md` first and write it last in every session.
- Report in the Slack channel in the short format "what works, what failed, what is needed".

## How to run

```bash
cd altermagagent
python3.14 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/pip install -e .
```

```bash
altermag import <root> --name <project>   # read-only scan of <root> into the index
altermag summary <project>                # print the nine summary points, the gaps and the runs in order
altermag benchmark <project> <file.json>  # compare with expected values: pass, fail or unknown per row
altermag gui                              # start the interface on http://localhost:8080
```

```bash
python demo/make_demo.py                  # writes the synthetic tree demo/CrSe_demo
python demo/make_demo.py --without-nosoc-bands   # negative fixture without the bands_nosoc run
```

## Real data

Diana's MnTe project lives outside the repository on Josh's laptop at `~/development/osi-hackathon/data/diana_MnTe/hackaton_exa`.
On other machines, import your own project folder instead.

## Day plan

The two-day plan is on Josh's laptop as `altermagagent-plan_v14.html`.
Day one is the folder reader and the storyline.
Scope is frozen for submission: see plan v14. Deferred after the hackathon: symmetry groups, planning of new runs, SSH access to clusters.
