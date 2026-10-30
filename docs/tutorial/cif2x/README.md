# cif2x tutorial files (Quantum ESPRESSO, Co3SnS2)

Input set for the Quantum ESPRESSO walk-through in the user guide
("Tutorial" chapter). Run from this directory:

```bash
cif2x -t QE input.yaml Co3SnS2_nosym.cif
```

The scf input is written to `scf/scf.in`. The committed `scf/scf.in` is a
reference result produced with the pslibrary UPF files in place (the cutoff
values in it come from those files).

Files:

- `input.yaml`: the input parameter file (template + `content`, pseudopotential
  index in `optional`).
- `scf.in_tmpl`: the namelist template merged into the generated input.
- `Co3SnS2_nosym.cif`: the crystal structure.
- `pseudo/pp_psl_pbe_rrkjus.csv`: the element-to-pseudopotential index. The
  UPF files themselves are **not** bundled: download the pslibrary PBE RRKJUS
  set into `pseudo/` with the file names listed in the CSV. The cutoffs are
  then read from the UPF headers. Once the files are there you can also build
  a cutoff index with `cif2x-pp-cutoff -o pseudo/cutoff.csv pseudo/*.UPF` and
  add `cutoff_file: ./pseudo/cutoff.csv` to `optional`, so that later runs do
  not need to open the UPF files.
- `respack/`: the RESPACK workflow example (SrVO3), see its own README.
