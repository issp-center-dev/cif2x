# cif2x RESPACK target (`-t respack`) — design

Date: 2026-06-30

## Goal

Add a `respack` target so `cif2x -t respack input.yaml structure.cif` emits the
whole RESPACK workflow in one run: the Quantum ESPRESSO inputs RESPACK consumes
(`scf` / `nscf` / `bands`, via the existing QE generator) **and** the RESPACK
control input `input.in` (a new generator). The structure-derived fields of
`input.in` are auto-filled (number of Wannier functions, the high-symmetry
k-path, the interpolation grid); physics parameters come from a template /
`content`. Initial guesses use **SCDM** (`N_initial_guess = 0`), so no manual
Gaussian projection block is emitted.

Scope: the QE target stays unchanged; `respack` reuses it for QE-mode tasks.

## Background (verified)

- Targets are dispatched in `main.py` via `normalize_target` + a `Struct2X`
  class; `input_validator.TARGETS` holds per-target aliases/required/allowed
  keys. Each `cif2x` run uses ONE target, but iterates a `tasks:` list.
- `Struct2QE` already computes the Wannier counts: `_find_nbnd_info`
  (`struct2qe.py:179-204`) sets `self.num_wann`, `self.nexclude`,
  `self.num_bands` from `pp_file` `orbitals`/`nexclude` (loaded by
  `_set_pseudo_info`). These are exactly `input.in`'s `N_wannier` and the QE
  `nbnd` for the nscf step.
- `cards._generate_band_path` builds a high-symmetry path from pymatgen
  `HighSymmKpath` (added for `bands`); `input.in`'s `&param_interpolation`
  block needs the same high-symmetry points (without per-segment divisions).
- RESPACK `input.in` (reference: `samples/*/input.in` in the user's
  `respack-wannier-py`) is Fortran namelists + trailing non-namelist blocks:
  `&param_chiqw` (Ecut_for_eps, flg_cRPA), `&param_wannier` (N_wannier,
  Lower/Upper_energy_window, N_initial_guess) [+ optional Gaussian initial-guess
  lines], `&param_interpolation` (dense, N_sym_points) + a high-symmetry k-coords
  block, `&param_visualization`, `&param_calc_int`. Per the reference parser
  (`input_parser.py:367-372`), `N_initial_guess = 0` means SCDM-only and the
  Gaussian block is omitted.
- `f90nml` is already a dependency; `Struct2QE` parses QE templates via
  `QEInputGeneral` (qe_tools), which is QE-specific — NOT reused for RESPACK.

## Decisions

- `-t respack` routes **per task by `mode`**: QE modes
  (`scf`/`nscf`/`bands`/`relax`/`vc-relax`/`dos`) → `Struct2QE`; `mode: respack`
  → new `Struct2RESPACK`; anything else → `InputValidationError`.
- `input.in` is built **template + content merge** (like QE): the template holds
  the `&param_*` namelists (physics); cif2x fills the structure-derived blanks
  and appends the auto-generated k-path block.
- **SCDM** initial guess: set `N_initial_guess = 0`, emit no Gaussian block (v1).
- Wannier-count computation is **factored into a shared helper** so QE and
  RESPACK produce identical `num_wann`/`nexclude`.

## Components

### 1. Shared Wannier-count helper

Extract the `num_wann`/`nexclude` computation currently inside
`Struct2QE._find_nbnd_info` into a reusable function (new
`src/cif2x/wannier.py`, or a function in `cif2x/qe/tools.py`), e.g.
`wannier_counts(species, pp_list, *, spinor_factor) -> {num_wann, nexclude, num_bands, nbnd}`.
`Struct2QE._find_nbnd_info` calls it (behavior unchanged — covered by existing
QE tests); `Struct2RESPACK` calls it for `N_wannier`.

### 2. `Struct2RESPACK` (`src/cif2x/struct2respack.py`)

`__init__(params, struct)` then `write_input(output_file, output_dir, dry_run=False)`,
matching the other generators (and supporting `--dry-run` via `utils.dryrun_emit`).

Build `input.in`:
1. Load the RESPACK template (`params["template"]`) namelists with `f90nml`;
   deep-merge `params["content"]` on top (reuse `deepupdate`).
2. Resolve pseudopotential / orbital info (reuse `Struct2QE._set_pseudo_info`
   logic or share it) and compute `num_wann` via the shared helper.
3. Auto-fill:
   - `&param_wannier`: `N_wannier = num_wann`; `N_initial_guess = 0` (SCDM);
     `Lower_energy_window` / `Upper_energy_window` from content (physics).
   - `&param_interpolation`: `dense` from content (or a k-resolution default);
     `N_sym_points = <number of high-symmetry points>`.
   - `&param_chiqw` / `&param_calc_int` / `&param_visualization`: passthrough
     from template/content.
4. Render: write the `&param_*` namelists in RESPACK order, and **after
   `&param_interpolation`'s `/`, append the high-symmetry k-coords block** — one
   line per point `kx ky kz  ! label`, count = `N_sym_points`, from
   `HighSymmKpath`. (Cross-check the exact RESPACK `m_rdinput` trailing-block
   format/order against the reference `respack-wannier-py` `input_parser.py`
   during implementation.)

### 3. Dispatch + validation

- `main._run`: when `target == "respack"`, select the generator per task:
  `mode in _QE_MODES → Struct2QE`, `mode == "respack" → Struct2RESPACK`, else
  `InputValidationError("task N: unknown mode '<m>' for target 'respack'")`.
  Other targets keep the single-generator loop. Factor the per-task body
  (params build, optional injection, write_input) so both paths share it.
- `input_validator.TARGETS["respack"]`: `aliases = ("respack",)`;
  per-task allowed keys = QE allowed ∪ respack-specific (`template`, `content`,
  `output_file`, `output_dir`, `mode`, `optional`); required = `mode`,
  `output_file`. Mode legality (QE-mode vs `respack`) is checked at dispatch.

### 4. QE side (no new code)

The `scf`/`nscf`/`bands` tasks are ordinary `Struct2QE` tasks. For RESPACK the
nscf task uses `K_POINTS` `crystal` (uniform explicit mesh, already supported)
and `nbnd` is auto-filled from `num_bands` (already supported). Configuration
lives in the sample `input.yaml`.

## Error handling

- `respack` task with an unrecognized `mode` → `InputValidationError` (clean exit).
- Missing `pp_file`/`orbitals` when `num_wann` is needed → the existing
  `Struct2QE` error path (reused), surfaced clearly.
- Template parse failure → `InputValidationError` with the file name.

## Testing (network-free where possible)

- **Shared helper:** `wannier_counts` returns the expected `num_wann`/`nexclude`
  for a small species list + a stub pp table; existing QE `nbnd` tests still pass
  (behavior unchanged).
- **Struct2RESPACK:** with a small structure (`use_ibrav=False`), a minimal
  RESPACK template, and a stub pp table, assert the rendered `input.in` has
  `N_wannier` = computed num_wann, `N_initial_guess = 0`, no Gaussian block, a
  `&param_interpolation` k-coords block whose count matches the high-symmetry
  points, and that template/content physics values pass through. Use
  `pytest.importorskip("pymatgen")`.
- **Dispatch:** `-t respack` routes a `mode: scf` task to `Struct2QE` and a
  `mode: respack` task to `Struct2RESPACK` (monkeypatch both, assert each
  `write_input` is called for the right task); an unknown mode raises.
- **Validation:** respack target accepts the sample tasks; an unknown top-level
  or task key is rejected.

## Documentation + sample

- Command reference (en/ja): add `respack` to the targets list.
- A new section / tutorial example showing the full chain `input.yaml`
  (`tasks:` = scf, nscf, bands QE tasks + a `mode: respack` task) plus a RESPACK
  `input.in` template, and noting SCDM (`N_initial_guess = 0`) and that energy
  windows / Ecut / flg_cRPA are user physics.
- A `sample/cif2x/respack/<case>/` directory mirroring the existing samples
  (structure file + input.yaml + RESPACK template), if a runnable example is
  wanted (decide in the plan; may be docs-only inline like the sweep examples).

## Out of scope (v1)

- Manual Gaussian initial-guess projections (SCDM only).
- Spin-polarized split `input.in.up` / `input.in.dn` (LSDA) — note as a
  follow-up.
- `qe2respack` invocation (a runtime conversion step the user runs; cif2x only
  emits inputs).
- Auto-deriving energy windows / cRPA parameters (physics; template/content).
- Cross-tool orchestration beyond per-task routing within one `-t respack` run.

## Open items to resolve in the plan (via cross-check)

- Exact RESPACK `input.in` trailing-block format/order (k-coords after
  `&param_interpolation`; whether labels are allowed inline) — verify against the
  reference `respack-wannier-py` `input_parser.py` / `m_rdinput`.
- Whether to reuse `Struct2QE._set_pseudo_info` directly (composition) or extract
  a shared pp loader alongside the `wannier_counts` helper.
