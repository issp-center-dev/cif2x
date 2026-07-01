# SrVO3 RESPACK sample

`cif2x -t respack input.yaml SrVO3.cif` emits `scf.in`, `nscf.in`, and `input.in`.
Provide `SrVO3.cif`, `pp.csv`, and the pseudopotentials under `./pseudo`.
`N_wannier = 3` (t2g) and the energy windows are physics choices in
`respack.in_tmpl`. Run order: scf → nscf → `qe2respack` → `respack` (the
`respack-wannier-py` port; SCDM, `N_initial_guess = 0`).
