# Changelog

## 1.2.0 (unreleased)

### Changed behaviour

- The `qe-tools` dependency is gone; Quantum ESPRESSO templates are read with
  f90nml. Indexed keys such as `starting_magnetization(2)` keep their index.
- The `phonopy` dependency is gone as well: nothing in cif2x or getcif
  imported it.
- A task-level `optional` block now overrides the global `optional` section
  (it was the other way round).
- `k_resolution` derives the k mesh from the reciprocal lattice vector
  lengths, which changes the mesh for non-orthogonal cells.
- The `nbnd` estimate is no longer doubled for collinear (`nspin = 2`) runs.
- VASP: the POTCAR must be resolvable; a missing or unreadable POTCAR is an
  error. Set `content.potcar.skip: true` to omit it deliberately. Only `INCAR`
  and `KPOINTS` are read from `template_dir`.
- `--api-key-file` / `api_key_file` accept any file name; an explicitly given
  file that does not exist is an error.
- `${range(...)}` with a float step rounds its values and never includes the
  end value; a `${...}` inside a longer string is an error.

### New

- Magnetic moments of mCIF structure data are carried into the inputs:
  `nspin`/`noncolin` and `starting_magnetization(i)` (plus `angle1`/`angle2`)
  for Quantum ESPRESSO, `MAGMOM` with `ISPIN`/`LNONCOLLINEAR` for VASP, and
  up/down electrons with `scf.SpinPolarization` for OpenMX. Antiferromagnetic
  cells survive `use_primitive` and `use_ibrav`.
- `pp_file` takes an optional `zval` column (valence charge); otherwise the
  UPF header is read (v1 and v2).
- OpenMX: `optional.vps_variant` chooses the hard/soft pseudopotential of
  Fe, Co, Ni, Cu, Zn.
- `cif2x-pp-cutoff` (formerly `src/utils/pp_cutoff.py`) is installed as a
  command and builds the `cutoff_file` CSV from UPF files.
- A "Parameters for RESPACK" section in the manual, and corrected VASP
  `kpoints` keys, QE card options, appendix examples, and structure keys.

### Fixed

- `cif2x -h`, `cif2x --version` and the getcif equivalents return at once:
  the generators and pymatgen are imported on first use instead of at
  start-up (issue #1; help took about 4 s before).
- Every user-caused failure (missing files, disordered structures, missing
  pseudopotential rows, bad options) ends with a one-line error and exit
  status 1 instead of a traceback; write failures are reported.
- A bare `ATOMIC_POSITIONS:` / `K_POINTS:` in `content` no longer crashes.
- `getcif --version` printed `develop` because of a wrong distribution name.
- The Bi "precise" basis in the OpenMX table had a typo.
- pyakaikkr's progress output no longer leaks into `--dry-run` output.
- Floating-point noise no longer appears in OpenMX coordinates nor splits
  magnetic species.

### Development

- Tests on Python 3.11 to 3.14; `ruff check src tests` in CI; the shipped
  samples are regenerated and compared with their references in the test
  suite (the AkaiKKR samples only where pyakaikkr is installed, which CI
  does not do).

## 1.1.0 (2024-09-14)

Earlier releases are described in the user guide.
