# SrVO3 RESPACK sample

`cif2x -t respack input.yaml SrVO3.cif` emits `scf.in`, `nscf.in`, and
`input.in`. The structure (`SrVO3.cif`), the pseudopotential mapping
(`pp.csv`), and the cutoffs (`cutoff.csv`: ecutwfc 80 Ry / ecutrho 320 Ry —
ONCV UPF headers carry no usable cutoff metadata) are included; only the
pseudopotential files themselves are not bundled. Download the SG15 ONCV set
(PBE, v1.0) and rename to the `<element>.<name>.UPF` convention:

```bash
mkdir -p pseudo && cd pseudo
for el in Sr V O; do
  curl -LO "http://www.quantum-simulation.org/potentials/sg15_oncv/upf/${el}_ONCV_PBE-1.0.upf"
  mv ${el}_ONCV_PBE-1.0.upf ${el}.ONCV_PBE-1.0.UPF
done
```

`N_wannier = 3` (t2g) and the energy windows are physics choices in
`respack.in_tmpl`. Run order: scf → nscf → `qe2respack` → `respack` (the
`respack-wannier-py` port; SCDM, `N_initial_guess = 0`). See the manual's
"RESPACK workflow" tutorial (`docs/tutorial/cif2x/respack/`) for a walkthrough
with the generated outputs included.
