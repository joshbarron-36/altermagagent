# Benchmark of expected values

This page lists the values that AltermagAgent must read from two project trees.
Each row names the value, the source file, the units, the tolerance and what the comparison tests.
A benchmark file is JSON; the command compares it with an imported project.

```bash
altermag benchmark <project> <benchmark.json>
```

Each row prints PASS, FAIL or UNKNOWN.
UNKNOWN means that the expected value is not filled in, or that the value cannot be read from the files.
The command exits with code 1 when a row fails.

## Synthetic CrSe demo tree

Write the tree with `python demo/make_demo.py`, import it, then run `altermag benchmark CrSe_demo demo/benchmark_crse.json`.
The values come from the synthetic files that `demo/make_demo.py` writes.

| Value | Expected | Source file | Units | Tolerance | What the comparison tests |
|---|---|---|---|---|---|
| Net moment | 0.00 | scf_nosoc/OUTCAR | μB per cell | 0.01 | The summed site moments cancel |
| Cr site moment | 3.12 | scf_nosoc/OUTCAR | μB | 0.05 | The OUTCAR magnetization table is read per site |
| Space group number | 194 | scf_nosoc/CONTCAR | none | 0 | spglib (symprec 0.01) finds P6_3/mmc |
| Lattice a | 3.684 | scf_nosoc/CONTCAR | Å | 0.001 | The relaxed lattice is read |
| Lattice c | 6.019 | scf_nosoc/CONTCAR | Å | 0.001 | The relaxed lattice is read |
| Spin splitting | 0.42 | bands_nosoc/vasprun.xml | eV | 0.01 | Largest up-down difference among the 6 sampled k-points and bands within ±1.5 eV of E_F |
| SCF energy | −36.9027 | scf_nosoc/OSZICAR | eV | 0.0001 | The last F= value of OSZICAR is read |

The splitting value is a largest up-down difference among sampled k-points along the k-path.
It is not a proof of splitting over the full Brillouin zone.

## Diana's MnTe project

The template is `docs/benchmark_mnte_template.json`.
Diana fills every `expected` value; until then every row prints UNKNOWN.
The field `tool_value_2026_10_21` holds what the tool read on 21 October, for comparison only.

| Value | Tool value on 21 Oct | Source file | Units | Tolerance | What the comparison tests |
|---|---|---|---|---|---|
| Net moment | 0.001 | 9-scf/OUTCAR | μB per cell | 0.01 | The summed site moment vectors cancel |
| Mn site moment | 2.246 | 9-scf/OUTCAR | μB | 0.05 | The OUTCAR magnetization table is read per site |
| Space group number | 194 | 9-scf/CONTCAR | none | 0 | spglib (symprec 0.01) finds P6_3/mmc |
| Lattice a | 4.1475 | 9-scf/CONTCAR | Å | 0.001 | The structure file is read |
| Lattice c | 6.71 | 9-scf/CONTCAR | Å | 0.001 | The structure file is read |
| Symmetry operations | 24 | 9.6-wannier/wannier_sym/mysymmop.dat | operations | 0 | The wannhr_symm log is read; not a check of the magnetic state |
| Spin splitting | none | a collinear non-SOC band run, not in the folder | eV | 0.01 | UNKNOWN until that run exists |
| SCF energy | −22.667769 | 9-scf/OSZICAR | eV | 0.0001 | The last F= value of OSZICAR is read |

Diana confirmed the Mn moments and the space group on 21 October; the other values await her check.

## Negative fixture

`python demo/make_demo.py --without-nosoc-bands` writes `demo/CrSe_demo_without_nosoc_bands`, the same tree without the bands_nosoc run.
On that tree no collinear non-SOC run has eigenvalues in vasprun.xml.
The diagnosis must list non-relativistic spin splitting under "Not supported by these files".
The benchmark row "spin splitting" must print UNKNOWN; every other row of `demo/benchmark_crse.json` must pass.

```bash
python demo/make_demo.py --without-nosoc-bands
altermag import demo/CrSe_demo_without_nosoc_bands --name CrSe_negative
altermag benchmark CrSe_negative demo/benchmark_crse.json
```

The test `tests/test_negative_fixture.py` checks both trees in a temporary folder with a temporary index.

```bash
python -m unittest discover -s tests
```

## Quantities a row can name

| quantity | Extra fields | Value read |
|---|---|---|
| net_moment | run (optional) | Size of the summed site moments in OUTCAR |
| site_moment | run, element | Largest absolute site moment of that element in OUTCAR |
| spacegroup_number | run | spglib space group number of the structure file |
| lattice_a, lattice_c | run | Lattice lengths of the structure file |
| splitting_max | none | Largest up-down difference of the spin splitting point |
| symmops | none | Operations listed in mysymmop.dat |
| total_energy | run | Last energy from OSZICAR or OUTCAR |
