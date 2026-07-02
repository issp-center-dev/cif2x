# RESPACK tutorial fixtures (SrVO3)

Reproducible input set for the "RESPACK workflow" section of the cif2x
tutorial. `cif2x -t respack input.yaml SrVO3.cif` emits `scf/scf.in`,
`nscf/nscf.in`, and `input.in` (all committed here for the manual to include).

Pseudopotentials are **not** bundled. Download the SG15 ONCV set (PBE, v1.0)
and rename the files to the `<element>.<name>.UPF` convention cif2x expects:

```bash
mkdir -p pseudo && cd pseudo
for el in Sr V O; do
  curl -LO "http://www.quantum-simulation.org/potentials/sg15_oncv/upf/${el}_ONCV_PBE-1.0.upf"
  mv ${el}_ONCV_PBE-1.0.upf ${el}.ONCV_PBE-1.0.UPF
done
```

`pp.csv` maps the elements to these files; `cutoff.csv` pins
`ecutwfc = 80 Ry` / `ecutrho = 320 Ry` (the ONCV UPF headers carry no usable
cutoff metadata, so the values are supplied via `optional.cutoff_file`).
`N_wannier = 3` (V t2g) and the energy window in `respack.in_tmpl` are physics
choices. Run order: `pw.x` (scf) → `pw.x` (nscf) → `qe2respack` → RESPACK.
