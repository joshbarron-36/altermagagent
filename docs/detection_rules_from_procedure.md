Source: derived from Diana Golovanova's working notes (docs/dft_run_procedure.md), shared 21 October 2026. Private to the team until she says otherwise.

# Detection rules from Diana's DFT procedure

Each rule names a phase, the files or INCAR tags that mark it, the checks her procedure states, and what the summary says.
Rules marked "implemented" are in `altermagagent/importer.py` (field `phase`, `steps`, `checks`) and `altermagagent/summary.py`.
Rules marked "documented only" are not yet machine-checked.

## How a run is found

- R0. A run is a folder with INCAR plus POSCAR or CONTCAR (implemented).
- R0a. Sub-folders without INCAR inside a run are tool folders (for example `wannier_sym`, `AHC`); their files are attached to the run as extra steps (implemented, up to three levels deep).
- R0b. A tool folder that is not inside any run becomes its own run with phase symmetrisation, transport or projection (implemented).
- Phase values: relax, scf, band, wannier, symmetrisation, projection, transport, other.
- The run card shows the phase plus extra steps, for example "wannier + symmetrisation, transport".

## Phase rules, in the order they are tested

### R1. WANNIER90 (phase `wannier`, implemented)

- Marked by: `LWANNIER90 = .TRUE.` in INCAR, or `wannier90.win` in the folder, or `NUM_WANN` in INCAR.
- Her INCAR for this step: `ISTART = 0, ICHARG = 11, ISYM = -1, LWANNIER90 = .TRUE., LWAVE = .FALSE., LCHARGE = .FALSE.`, KPOINTS regular mesh.
- Files she names: wannier90.win, wannier90.amn, wannier90.mmn, wannier90.eig, wannier90_hr.dat, wannier90.wout, run.wexac.sh, submit_wannier.sh, repair_amn.py.
- Checks (implemented):
  - wannier90.amn, wannier90.mmn and wannier90.eig exist ("Check that the .mmn .amn and .eig files have been created").
  - NUM_WANN is set in INCAR.
  - wannier90.win `num_bands` equals INCAR NBANDS.
  - wannier90.win `num_wann` equals INCAR NUM_WANN.
  - wannier90.win has `write_hr = .true.`.
  - `#BSUB -n` (or `#SBATCH -n`) in the job script (run.wexac.sh or another *.sh) equals NBANDS ("set the number of nodes equal to NBANDS").
  - `#BSUB -n` in submit_wannier.sh equals NUM_WANN.
  - When symmetrisation follows, the VASP path in the job script contains 5.4.4 (the group's recommendation); reported as info only.
- Checks (documented only): `write_xyz = .true.`; disentanglement window keys dis_win_max, dis_win_min, dis_froz_max, dis_froz_min; projections in uudd spin order; `use_bloch_phase` commented.
- Summary when it fires: Hamiltonian point "generated, RUN/wannier90_hr.dat, DATE" when wannier90_hr.dat exists; "prepared" when only wannier90.win exists.

### R2. BAND (phase `band`, implemented)

- Marked by: `ICHARG = 11` plus any of KLABELS, BAND.dat, REFORMATED_BAND.dat or REFORMATTED_BAND.dat (vaspkit option 211), KLINES.dat, BAND_GAP, or a line-mode KPOINTS.
- Note: her text says REFORMATED_BAND.dat; vaspkit actually writes REFORMATTED_BAND.dat (seen in MnTe 9.2-band). Both spellings are accepted.
- Inputs she names: CHGCAR and CHG copied from SCF, KPOINTS along characteristic points, KLABELS.
- Outputs and scripts: BAND.dat, REFORMATED_BAND.dat, KLABELS from vaspkit 211; band_plot.m on the PC.
- Checks (implemented, info when missing): KLABELS present; BAND.dat present.
- Checks (documented only): CHGCAR in the folder came from the SCF parent; the last comment in KLABELS was deleted before plotting.
- Summary when it fires: the run counts as a bands run for the band figure tool and the SOC bands gap.

### R3. Other non-self-consistent step (phase `other`, implemented)

- Marked by: `ICHARG = 11` on a regular mesh with none of the band markers.
- In her procedure this is the Fermi surface step (folder such as `5-fermi`), with `ISTART = 0, ICHARG = 11, LWAVE = .F., LCHARG = .F., LREAL = .F., NWRITE = 3`.
- Program: pyprocar `fermi2D` on PROCAR, INCAR, OUTCAR, KPOINTS, POSCAR copied to the PC.
- Summary when it fires: phase "other" with the reason "ICHARG = 11 on a regular mesh".

### R4. Other recognised steps (phase `other`, implemented)

- wt.in: WannierTools input for surface states (SlabSS_calc, SlabArc_calc); outputs read with `gnuplot surfdos_l.gnu`.
- trace.txt: vasp2trace output for band representations (folder such as `2_vasp2trace`).
- band.conf, band.yaml, phonopy_disp.yaml, POSCAR-001 and following: phonopy phonon modes.
- findsym.in: symmetr input for the linear response symmetry tensor; recorded as an extra step, it does not change the phase.

### R5. RELAX (phase `relax`, implemented)

- Marked by: NSW > 0 and IBRION not -1.
- Her values: `ISIF = 1` (3, or 4 for 2D), `PREC = Normal`, `IBRION = 2`, `NSW = 200`, `LORBIT = 11`, `EDIFFG = -1E-2`.
- Second cycle: `ISIF = 2, IBRION = 2, POTIM = 0.5, ADDGRID = .T., EDIFFG = -1E-4`; on error, `ISIF = 2, IBRION = 1` from CONTCAR.
- Checks (documented only): forces after POSITION in OUTCAR are zero; KPOINT_num*LENGTH_a in [15..30]; vacuum no more than 20 Å for slabs.

### R6. SCF (phase `scf`, implemented)

- Marked by: everything else with an INCAR, normally `ICHARG = 2, ISTART = 0, ISIF = 2, IBRION = -1, NSW = 0, LCHARG = .TRUE.`.
- Outputs she names: CHGCAR, CHG, WAVECAR (WAVECAR only with `LWAVE = .TRUE.` and fixed NBANDS).
- Inputs: POSCAR, POTCAR (element order as in POSCAR), INCAR, KPOINTS, run.wexac.sh.
- Checks (implemented): EDIFF at most 1E-6 ("EDIFF = 1E-6 is enough for good convergence"), shown as a warning otherwise.
- Checks (documented only):
  - ENCUT equals the largest ENMAX in POTCAR times 1.3.
  - NELMIN 4 to 10, and the dE in OUTCAR approaches 10^-1.
  - NELMDL: if dE at the start of the SCF did not reach 10^-5, NELMDL should be enlarged (negative, magnitude above 2).
  - IALGO = 38 for big systems.
  - Restart from a non-converged run uses ICHARG = 0, ISTART = 1.
  - KPOINT_num*LENGTH_a in [30..50].
  - Non-magnetic SOC: LSORBIT = T, ISPIN = 2, MAGMOM = 3*NATOMS zeros.

### R7. VASP binary (all VASP phases, implemented)

- Marked by: the `vasp_std`, `vasp_ncl` or `vasp_gam` path on an uncommented line of run.wexac.sh.
- Check: SOC (LSORBIT = .TRUE.) or LNONCOLLINEAR needs vasp_ncl; otherwise std. Mismatch is a warning.

### R8. Symmetrisation with wannhr_symm (phase `symmetrisation`, implemented)

- Marked by any of: protmat.log, mysymmop.dat, magmom_*.log (magmom_primed.log, magmom_unprimed.log), new_rpt_block.dat.
- Also marked by: wannier90_hr.dat_*nsymm* (for example wannier90_hr.dat_unprimed_nsymm_4) and the two symmhr_addrptblock scripts.
- Inputs that go with them but do not mark the step alone: poscar.in, wann.in, locaxis.in, magmom.in (magnetic case only), wannier90_hr.dat, run.wexac.sh.
- Non-magnetic command: `python2 symmhr_addrptblock.py` in a py27 environment.
- Magnetic commands: `python2 symmhr_addrptblock_magnetic_conjugate_when_primed.py 1.0 wannier90_hr.dat`.
- Then the same script with `-1.0 wannier90_hr.dat_unprimed_nsymm_4`.
- Extra data read (implemented): the number of operations listed in mysymmop.dat.
- Checks (documented only): Wannier basis of wannier90_hr.dat is in uudd order, all orbitals of one atom before the next; magmom.in exists when MAGMOM is not zero; bands checked after symmetrisation.
- Summary when it fires: Hamiltonian and Symmetries points say "symmetrised Hamiltonian present, wannhr_symm logs found (FOLDER: logs, mysymmop.dat lists N operations)".
- Gap when a Hamiltonian exists without these files: "Hamiltonian not symmetrised: no wannhr_symm logs found".

### R9. Orbital projections (phase `projection`, implemented)

- A group projection method A: parse_vaspxml.ini, SELECTED_ATOMS_LIST (vaspkit 214), `python /path/to/tools/parse_vaspxml_v7.py` in the band folder.
- Its outputs are PBAND_atom-… and DOS_atom-…, plotted with plot_pband_pdos.m.
- Method B: vaspkit 213 and projband.py in the band folder.
- Marked by: PBAND_*, DOS_atom*, SELECTED_ATOMS_LIST, parse_vaspxml.ini.
- In a band folder the run keeps phase `band` and shows "+ projection"; a separate folder gets phase `projection`.
- Check (documented only): parse_vaspxml_v7.py is not correct for split SPIN_UP and SPIN_DOWN.

### R10. Transport, SHC and AHC (phase `transport`, implemented)

- Her SHC procedure inputs: wannier90_hr.dat, wpos (Wannier centres in Cartesian coordinates from wannier_name_center.txt, via calc_wann.m), setting, 1-submit.sh, 2-submit.sh, 3-submit.sh.
- Outputs: files in the `out` and `result` folders; MnTe holds ahc.dat with scripts 2_rmerge.sh, 3_get_shc.sh, 4_get_ahc.sh.
- Marked as calculated by: ahc*.dat, shc*.dat, anything under a `result` folder, and the postw90 or WannierBerri outputs (*.wpout, *-ahc*, *_kubo*, *boltzwann*, *-berry*, *wberri*).
- Marked as prepared by: setting plus 1-submit.sh without outputs.
- Extra data read (implemented): the task line of `setting` (ahc, shc, fullahc, fullshc).
- Checks (documented only): wpos agrees with the centres in wannier90.wout; Wannier spin order is uuuddd; icompiler/2021.3 and iimkl/2023b modules loaded.
- Summary when it fires: Transport point "calculated, N output files (RUN/AHC/ahc.dat); setting task: fullahc in RUN/AHC".
- Berry curvature: her document names no separate Berry curvature output; AHC is the Berry curvature route here. The postw90 `-berry` outputs are still matched.

### R11. Steps documented only (not detected yet)

- Wilson loop: POSCAR, wannier90_hr.dat, wannier_name_center.txt, band_1D.m, Wilson_loop_hr_new.m (run on the PC).
- Symtopo: `symtopo geninp` and `symtopo caltopo`, folder such as `SOC/6-symtopo`.
- Charge density waves: get_phmode.py, phmode_q*b* files, plot_phonon_mode_nondegenerated.m, ph_modes.mat, conv_vectors.ipynb, *.vesta.
- Phonons with phonopy: 1_sub_vaspscf.sh, sample_input folder, plot_phonon_projection.m.
- Quantum Espresso and EPW: pw.x, ph.x, q2r.x, matdyn.x, EPW/bin/pp.py, bands.x, kmesh.pl, plotband.x, pb.save, epw2.in, scf_band.kpt.
- GGA+U: LDAU, LDAUTYPE = 2, LDAUL, LDAUU, LDAUJ, LDAUPRINT = 1, LMAXMIX = 6, LASPH = T, LORBIT = 11; AMIX from 0.2 on bad convergence.
- Slab: vacuum in POSCAR, VESTA transform, KPATH perpendicular to the slab.
- Parallelisation: IBZKPT gives the number of irreducible k-points; KPAR and NCORE split the cores.

## External programs she names

| Program | Used for |
|---|---|
| VASP (vasp_std, vasp_ncl; 6.3.2 and 5.4.4) | SCF, band, Wannier, relaxation, Fermi surface and phonon force runs |
| vaspkit | 211 band data (BAND.dat, REFORMATED_BAND.dat, KLABELS), 213 projected bands, 214 SELECTED_ATOMS_LIST |
| Wannier90 (wannier90.x, 3.1.0) | Maximally localised Wannier functions and wannier90_hr.dat |
| repair_amn.py | Fixing spacing in amn and mmn files after a segmentation fault |
| WannierTools (wannhr_symm utility, wt.x) | Hamiltonian symmetrisation; surface states and Fermi arcs |
| symmhr_addrptblock.py and symmhr_addrptblock_magnetic_conjugate_when_primed.py (python2, py27) | The wannhr_symm scripts, non-magnetic and magnetic |
| MATLAB scripts band_plot.m, plot_pband_pdos.m, calc_wann.m, band_1D.m, Wilson_loop_hr_new.m, plot_phonon_projection.m, plot_phonon_mode_nondegenerated.m | Plots and helper calculations on the PC |
| parse_vaspxml_v7.py (method A tools) | Orbital projections from vasprun.xml |
| projband.py (method B) | Orbital projections after vaspkit 213 |
| pyprocar | 2D Fermi surface cuts (fermi2D) and band unfolding |
| SHC/AHC code (1-submit.sh, 2-submit.sh, 3-submit.sh, Fortran built with ifort) | Spin and anomalous Hall conductivity from wannier90_hr.dat |
| phonopy | Phonon supercells and bands; space group and primitive cell (`--symmetry`) |
| get_phmode.py, conv_vectors.ipynb | Charge density wave displaced structures and their VESTA view |
| VESTA | Structure views and slab construction |
| vasp2trace and the Bilbao server | Band representations (trace.txt) |
| symtopo | Topology from DFT results |
| symmetr (linear-response-symmetry, conda env for_symmetr) | Linear response symmetry tensors from findsym.in |
| Quantum Espresso (pw.x, ph.x, q2r.x, matdyn.x, bands.x, plotband.x) and EPW | Phonons, electron-phonon coupling and Wannier interpolation |
| kmesh.pl (Wannier90 utility) | k-mesh lists for EPW |
| gnuplot | surfdos_l.gnu plots from WannierTools |
| LSF on the cluster (bsub, bjobs, bhist, bkill, bmod, bpeek, bstop, bresume, bswitch, blimits) | Job submission and control |
