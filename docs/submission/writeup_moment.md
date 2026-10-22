AltermagAgent reads old VASP folder trees without changing them and says which claims the files support, which they do not, and what to run next. On Diana Golovanova's real CrSe project it found two magnetic realizations in four branches, two unconverged SCF runs, and a magnetic space group that matches her own findsym.log.

## Approach

**The problem.** Computational materials projects live in scattered VASP folders whose history is written down nowhere. Months later, nobody can say which run supports which conclusion, or what is missing. Workflow tools record provenance only for work done inside them.

**What we built:**

- **Branch graph:** runs grouped by structure, magnetic realization and SOC, drawn like a git history; each parent link gives a reason, candidates or "unknown".
- **Evidence tags:** every summary point names its files and carries a tag: output found, parsed, validated or insufficient evidence.
- **Six-step suggestions** per branch, from structure to transport, each done, partial or missing, then the next calculation.
- **Define MSG:** spglib computes the magnetic space group from the moments and compares it with the user's findsym.log, or writes findsym.in.
- **Wannier comparison** by Diana Golovanova: Wannier bands against VASP bands, with mean and maximum error.
- **Literature kept apart:** papers found through the Valency MCP platform appear under "From the literature, not from these files" and never change a suggestion.
- **Chat** that answers only from the indexed files and cites them, and a **marimo** notebook that walks through the same analysis.

[IMAGE 1]
*The real CrSe project: eight tagged summary points and 13 runs in four branches.*

## Results

**CrSe** (13 runs, 1127 files):

- Two magnetic realizations: coplanar in ten runs (net 0.00 μB per cell, OUTCAR) and collinear in three (INCAR MAGMOM only).
- Four branches: coplanar non-SOC (7 runs), coplanar SOC (3), collinear non-SOC (1) and collinear SOC (2).
- On coplanar non-SOC, spglib gives P-62'm' (BNS 189.225), the same label as Diana's findsym.log.
- OSZICAR shows that 020_scf and 021_scf stopped at 200 and 160 steps above EDIFF; Diana confirmed both as known unconverged trials.
- On coplanar SOC, the Wannier bands match the VASP bands with a mean absolute error of 0.040 meV within 1 eV of the Fermi level (preliminary).

[IMAGE 2]
*CrSe coplanar SOC: Wannier-interpolated bands over the VASP bands.*

**MnTe** (4 runs, 175 files):

- Computed Cm'c'm (BNS 63.462) from the 9-scf OUTCAR moments, confirmed by Diana's own findsym.log.
- Anomalous Hall output (9.6-wannier/AHC/ahc.dat) was found and parsed, not validated.

[IMAGE 3]
*MnTe with SOC: band structure drawn from the 9.4-band run.*

**Measurements:**

- The 13 CrSe runs are scanned read-only in about four seconds on a laptop.
- On 20 synthetic CrSe trees with known truth (137 runs): phase 100 %, status 100 %, parent link 98.5 % (docs/accuracy.md).
- Poulami Sadhukhan installed it on a second machine, a Windows laptop, in about twenty minutes over a mobile hotspot, with no errors.

**What it refuses to claim:**

- It never states that a material is an altermagnet on the basis of the files.
- A computed magnetic space group stays "pending validation" until a findsym.log or the physicist confirms it.
- Transport output is reported as found and parsed, never as validated.
- Published papers never count as evidence.

## Next steps

- Spin group analysis next to the magnetic space group.
- SSH access to computing clusters, so the tool can read running projects where they live.
- Planning of new calculations: from the suggested next step to a prepared, checked input folder.

## Reproduce in five commands

```bash
python3.14 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e .
python demo/make_demo.py
altermag import demo/CrSe_demo --name CrSe
altermag summary CrSe
```

The demo tree is synthetic, so no VASP is needed; `altermag gui` opens the interface.

## Team and credits

- **Diana Golovanova**, Weizmann Institute of Science: idea, data, physics, procedure rules and the Wannier comparison module.
- **Josh Barron:** design and software build (directed with AI coding agents), integration of Diana's module and validations, video and write-up.
- **Poulami Sadhukhan**, CSIR National Chemical Laboratory: independent installation.

The CrSe and MnTe data were computed by Diana Golovanova with VASP (PAW, GGA-PBE) and Wannier90 through VASP2WANNIER90.
Symmetry analysis uses spglib, spinspg and FINDSYM; full references for these and the other codes are in the README.

With thanks to two hackathon sponsors: [Valency](https://valency.io), whose MCP platform powers the built-in literature search, and [marimo](https://marimo.io), which runs the walkthrough notebook.

## Links

- Repository: https://github.com/joshbarron-36/altermagagent
- Video: in the repository under docs/video
- Licence: Apache 2.0
