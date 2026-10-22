# Symmetry labels: what they mean and what Diana must confirm

Every label below is shown as "parsed, pending validation" until a physicist confirms it.
The code is in `altermagagent/symmetry.py`; the labels appear in the Symmetries point and in one diagnosis line.

## What each label means

- Space group: the crystal symmetry from the atom positions only, moments ignored (spglib, symprec 0.01 Å).
- Magnetic space group (BNS label): the operations that keep both the atoms and the converged OUTCAR moments, in the BNS setting.
- OG label: the same magnetic group written in the Opechowski-Guccione setting; for type 3 groups it equals the BNS label.
- uni number: spglib's serial number for the magnetic space group (1 to 1651).
- Type: 1 no time reversal at all, 2 grey group, 3 time reversal only with rotations or mirrors, 4 time reversal with a translation.
- Unitary operation: a space operation that keeps the moments without time reversal.
- Antiunitary operation: a space operation that keeps the moments only when time reversal is attached.
- Operation words: identity, inversion, mirror, glide, n-fold rotation, n-fold screw, rotoinversion, read from the rotation and translation.
- Moment input: x, y, z blocks of the last OUTCAR magnetization give axial vectors (SOC runs); one block gives scalars (collinear runs).
- Moment tolerance: 0.05 μB; site moments that differ by less are treated as equal.
- Collinear scalar moments: spglib does not rotate a scalar moment, so the collinear label carries no spin direction.
- wannhr_symm operations: the 24 operations listed in `wannier_sym/mysymmop.dat`, with the "MagSymop ? True" flags of the two magmom logs.
- wannhr_symm comparison: "consistent" when the operation set and the unitary and antiunitary subsets equal spglib's; else "differs" with the reason.
- Spin group (collinear non-SOC runs only): spinspg operations where the spin part may differ from the space part.
- Spin group label: spinspg gives operations, not a symbol, so the label is the point group G of all operations and the subgroup H that keeps the spin.
- Spin-only group: collinear means every rotation of the spins about the common axis is a symmetry.
- Opposite-spin relation: whether the operations that swap spin up and spin down include a pure translation or inversion, or only other operations.

## Results on Josh's laptop, 21 October

- MnTe, moments from 9-scf/OUTCAR (Mn 3.89, 2.25, 0 and its opposite): magnetic space group Cm'c'm, BNS 63.462, type 3.
- MnTe: 8 operations in the 4-atom cell, 4 unitary (identity, inversion, 2-fold screw, mirror), 4 antiunitary.
- MnTe: wannhr_symm lists 24 operations, the same set as spglib's P6_3/mmc (194).
- MnTe: the wannhr_symm logs keep 4 operations without time reversal (the same 4 as spglib) and 0 with time reversal.
- MnTe: spglib keeps 4 antiunitary operations (2-fold rotations, mirror and glide with time reversal), so the comparison prints "differs".
- CrSe demo scf_nosoc (synthetic, Cr +3.12 and −3.12): magnetic space group P6_3'/m'mc', BNS 194.269, type 3, 24 operations, 12 and 12.
- CrSe demo spin group: G 6/mmm, H −3m, spin-only group collinear, 12 of 24 operations flip the spin.
- CrSe demo: no pure translation or inversion swaps the two spin sublattices.

## What Diana must confirm

1. The BNS label Cm'c'm (63.462) for MnTe with the moments as in 9-scf/OUTCAR.
2. The operation count: 24 in mysymmop.dat for P6_3/mmc, and 8 for the magnetic group in the 4-atom cell.
3. The unitary versus antiunitary split: spglib 4 and 4, wannhr_symm logs 4 and 0; which one is right, and why the primed log keeps none.
4. The spin group for the CrSe demo: G 6/mmm with H −3m, and the opposite-spin relation stated above.
5. Whether a collinear magnetic space group label without a spin direction (CrSe demo) should be shown at all.
