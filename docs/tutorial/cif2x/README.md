# cif2x tutorial files (Quantum ESPRESSO, Co3SnS2)

Input set for the Quantum ESPRESSO walk-through in the user guide
("Tutorial" chapter). Run from this directory:

```bash
cif2x -t qe input.yaml Co3SnS2_nosym.cif
```

The scf input is written to `scf/scf.in`; `scf/` holds the expected result.

Files:

- `input.yaml`: the input parameter file (template + `content`, pseudopotential
  index in `optional`).
- `scf.in_tmpl`: the namelist template merged into the generated input.
- `Co3SnS2_nosym.cif`: the crystal structure.
- `pseudo/pp_psl_pbe_rrkjus.csv`: the element-to-pseudopotential index. The
  UPF files themselves are **not** bundled; download the pslibrary PBE RRKJUS
  set into `pseudo/` (file names as listed in the CSV), or add
  `cutoff_file: ./pseudo/cutoff.csv` to `optional` with a CSV produced by
  `cif2x-pp-cutoff` so that the cutoffs can be resolved without the files.
- `respack/`: the RESPACK workflow example (SrVO3), see its own README.
