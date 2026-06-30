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

- **Target tool = `respack-wannier-py` (the Python port)**, not Fortran RESPACK.
  v1's SCDM/`N_initial_guess = 0` (no Gaussian block) is supported by the port
  per `input_parser.py:367-372`; the Fortran `calc_wannier` would reject it. This
  assumption is stated in the docs.
- `-t respack` routes **per task by `mode`**: the RESPACK workflow QE modes
  **`scf` / `nscf` / `bands`** → `Struct2QE`; `mode: respack` → new
  `Struct2RESPACK`; any other mode (incl. `relax`/`vc-relax`/`dos`) →
  `InputValidationError` (respack's surface is intentionally narrowed to the
  workflow it consumes).
- `input.in` is built **template + content merge** (like QE): the template holds
  the `&param_*` namelists (physics); cif2x fills the structure-derived blanks
  and appends the auto-generated multi-segment k-path block.
- **SCDM** initial guess: cif2x always sets `N_initial_guess = 0` and emits no
  Gaussian block (v1). If the template/content supplies a nonzero
  `N_initial_guess` or a Gaussian projection block, that is **rejected** with
  `InputValidationError` (v1 is SCDM-only; no silent override).
- **Template contract (v1): namelist-only.** A RESPACK template contains only
  `&param_*` namelists (parsed by `f90nml`). Any non-blank, non-comment line
  after the last namelist's `/` is rejected with `InputValidationError` — cif2x
  generates ALL non-namelist blocks (the k-coords block) itself, so a template
  must not carry a Gaussian or k-path block.
- **Spin scope (v1):** `spinor_factor` comes from a **dedicated respack task
  field `spinor`** (bool, default `false` → factor 1; `true` → factor 2), NOT
  from a QE `&system` block (which must not appear in `input.in`). A spin-polarised
  `nspin = 2` (LSDA) respack configuration is **rejected** with
  `InputValidationError` (the split `input.in.up`/`.dn` is out of scope); the docs
  require `spinor` to match the QE nscf task.
- Wannier-count computation is **factored into a neutral shared helper**
  (`src/cif2x/wannier.py`, not under `cif2x/qe/`) so QE and RESPACK produce
  identical `num_wann`/`nexclude`.

## RESPACK `input.in` format (resolved via cross-check of the reference parser)

From `respack-wannier-py` `input_parser.py:108-230` (`parse_band_path`,
`_read_band_path_coords`) and the `samples/*/input.in` corpus:

- `&param_interpolation` supports **multiple disjoint paths**. `N_sym_points` is
  a **list** — one positive integer per path (number of high-symmetry points on
  that path, each `>= 2`), with trailing-zero padding (Fortran fixed array up to
  10). A scalar means a single path.
- The k-coords block follows the `&param_interpolation` `/` (default
  `reading_sk_format = 0`): for each path, for each high-symmetry point, one line
  of three fractional (crystal) coordinates `kx ky kz`; an inline `! label`
  comment is allowed and ignored by the parser. cif2x emits format 0.
- `dense` is three integers (scalar → `(n,n,n)`, or a 3-list); `Ndiv` defaults to
  40 if omitted (left to template/content).
- `&param_wannier`: `N_wannier`, `Lower_energy_window`, `Upper_energy_window`,
  `N_initial_guess` (0 ⇒ SCDM, no Gaussian block).

Example (single path, srvo3-style):

```
&param_interpolation
dense = 8, 8, 8
N_sym_points=5/
0.50 0.50 0.50  ! R
0.00 0.00 0.00  ! G
0.50 0.00 0.00  ! X
0.50 0.50 0.00  ! M
0.00 0.00 0.00  ! G
```

For a structure whose `HighSymmKpath` path has two disjoint segments
(e.g. `[[Γ,X,M,Γ,R,X],[M,R]]`), emit `N_sym_points=6,2/` then the 6 coords of
the first segment followed by the 2 of the second — no interpolation bridges the
break.

## Components

### 1. Shared Wannier-count helper

Extract the `num_wann`/`nexclude` computation currently inside
`Struct2QE._find_nbnd_info` into a reusable function in a **neutral module**
`src/cif2x/wannier.py` (NOT `cif2x/qe/`, so RESPACK does not depend on the QE
namespace):

```
wannier_counts(composition: Mapping[str, int], pp_list, *, spinor_factor: int)
    -> {"num_wann", "nexclude", "num_bands", "nbnd"}
```

- `composition` = per-element atom counts (e.g. `{"Sr":1,"V":1,"O":3}`); `pp_list`
  is the pp table keyed by element symbol (the existing `Struct2QE.pp_list`).
- Formula (verbatim from current `_find_nbnd_info`,
  `struct2qe.py:188-203`): `nexclude = Σ_element count·pp[el].nexclude · sf`;
  `num_wann = Σ_element count·(Σ orbital_degeneracy) · sf`;
  `num_bands = num_wann*2 + (num_wann//2 + (num_wann//2)%2)`;
  `nbnd = nexclude + num_bands`; where `sf = spinor_factor`.
- A shared **pp loader** (thin extraction of `Struct2QE._set_pseudo_info`) returns
  this `pp_list`; both QE and RESPACK use it. Missing `orbitals`/`nexclude`
  columns raise `InputValidationError`.

`Struct2QE._find_nbnd_info` calls `wannier_counts` (behavior unchanged — covered
by existing QE tests); `Struct2RESPACK` calls it for `N_wannier` with
`spinor_factor` from the task's `spinor` field.

### 2. `Struct2RESPACK` (`src/cif2x/struct2respack.py`)

`__init__(params, struct)` then `write_input(output_file, output_dir, dry_run=False)`,
matching the other generators (and supporting `--dry-run` via `utils.dryrun_emit`).

Build `input.in`:
1. Load the RESPACK template (`params["template"]`, namelist-only) with `f90nml`;
   reject any non-blank/non-comment line after the last namelist
   (`InputValidationError`). Deep-merge `params["content"]` on top (reuse
   `deepupdate`). Namelist keys are case-insensitive (f90nml lowercases them);
   internal handling is lowercase and the rendered output uses lowercase keys
   (RESPACK reads keys case-insensitively, per the reference parser's `.lower()`).
2. Resolve pseudopotential / orbital info via the shared pp loader, reading
   `pp_file`/`pseudo_dir` from the **shared `optional:` block** (same source as QE
   tasks; injected per-task by `deepupdate`). Determine `spinor_factor` from the
   task's `spinor` field (default 1; `true` ⇒ 2); reject `nspin = 2`
   (`InputValidationError`). Compute `num_wann` via `wannier_counts`.
3. Auto-fill:
   - `&param_wannier`: `N_wannier = num_wann`; force `N_initial_guess = 0` (SCDM)
     and reject a template/content nonzero `N_initial_guess` or Gaussian block
     with `InputValidationError`. `Lower_energy_window`/`Upper_energy_window` come
     from template/content and are **validated** (both present, numeric,
     `lower < upper`), else `InputValidationError`.
   - `&param_interpolation`: `dense` from template/content, validated as three
     ints (scalar → `(n,n,n)` or a 3-list); no invented physics default. The
     k-path comes from `HighSymmKpath(struct.structure)`: emit
     `N_sym_points = [len(seg) for seg in kpath["path"]]` (one entry per disjoint
     segment, each `>= 2`) and the segment-by-segment coords block.
   - `&param_chiqw` / `&param_calc_int` / `&param_visualization`: passthrough from
     template/content.
4. Render (format resolved above): write the `&param_*` namelists in RESPACK
   order; after `&param_interpolation`'s `/`, append the multi-segment k-coords
   block (format 0): per segment, one `kx ky kz  ! label` line per high-symmetry
   point (fractional coords). Surface pymatgen's "path may be incorrect" warning
   through the logger (as the QE bands path does) when the cell is non-standard.

### 3. Dispatch + validation

- `main._run`: when `target == "respack"`, select the generator per task:
  `mode in {"scf","nscf","bands"} → Struct2QE`, `mode == "respack" →
  Struct2RESPACK`, else `InputValidationError("task N: unknown mode '<m>' for
  target 'respack'")` (relax/vc-relax/dos are NOT accepted under respack). Other
  targets keep the single-generator loop. Factor the per-task body (params build,
  optional injection, write_input) so both paths share it.
- `input_validator` validates respack tasks **per mode**: `scf`/`nscf`/`bands`
  tasks must satisfy the existing QE task schema (so QE required/allowed keys still
  apply), and a `mode: respack` task uses the RESPACK schema (allowed:
  `mode`, `template`, `content`, `output_file`, `output_dir`, `optional`,
  `spinor`; required: `mode`, `output_file`, `template`). Any other mode →
  `InputValidationError`.

### 4. QE side (no new code, but documented requirements)

The `scf`/`nscf`/`bands` tasks are ordinary `Struct2QE` tasks. For RESPACK the
nscf task uses `K_POINTS` `crystal` (uniform explicit mesh, already supported)
and `nbnd` is auto-filled from `num_bands` (already supported). **`qe2respack`
requires the nscf run to set `nosym = .true.` and `noinv = .true.`** — the sample
`input.yaml` sets these in the nscf task's `content.system`. When the respack
target sees a `mode: nscf` task whose `content.system` lacks `nosym`/`noinv`
(or sets them false), it emits a `logger.warning` (not an error — the user may
have reasons), and the docs call out the requirement. (Auto-injection is a
possible follow-up.)

## Error handling

- `respack` task with an unrecognized/unsupported `mode` → `InputValidationError`.
- Template with non-blank/non-comment content after the last namelist →
  `InputValidationError` (v1 templates are namelist-only).
- Missing `pp_file`/`orbitals` columns when `num_wann` is needed →
  `InputValidationError` (the shared pp loader/helper raises it).
- `N_initial_guess != 0` or a Gaussian block in template/content →
  `InputValidationError` (v1 is SCDM-only).
- `nspin = 2` in a respack task → `InputValidationError` (LSDA out of scope).
- Missing/invalid energy windows (`Lower`/`Upper` absent, non-numeric, or
  `lower >= upper`) → `InputValidationError`.
- `dense` not three ints → `InputValidationError`.
- Template parse failure → `InputValidationError` with the file name.

## Testing (network-free where possible)

- **Shared helper:** `wannier_counts` returns the expected
  `num_wann`/`nexclude`/`num_bands`/`nbnd` for a small species list + a stub pp
  table, for both `spinor_factor=1` and `2`; existing QE `nbnd` tests still pass
  (behavior unchanged).
- **Struct2RESPACK render (golden-style):** with a small structure
  (`use_primitive=True`, `use_ibrav=False`), a minimal RESPACK template, and a
  stub pp table, render `input.in` and **parse it back** (f90nml for the
  namelists + a small reader for the k-coords block) to assert: `N_wannier` =
  computed num_wann; `N_initial_guess = 0` and no Gaussian lines; `N_sym_points`
  equals the per-segment lengths from `HighSymmKpath` and the coords block has
  that many lines per segment; `dense` and the windows pass through. Include a
  multi-segment structure (e.g. simple cubic) so the segment list is exercised.
  Use `pytest.importorskip("pymatgen")`.
- **Spin:** `spinor: true` doubles `num_wann`/`N_wannier` vs the default; assert
  both.
- **Template contract:** a template with a trailing non-namelist line is rejected.
- **Validation/error paths:** nonzero `N_initial_guess`, missing window,
  `lower >= upper`, bad `dense`, `nspin = 2` each raise `InputValidationError`.
- **Dispatch:** `-t respack` routes a `mode: scf` task to `Struct2QE` and a
  `mode: respack` task to `Struct2RESPACK` (monkeypatch both, assert each
  `write_input` is called for the right task); `mode: relax` (and unknown) raise.
- **Validation:** respack target accepts the sample tasks; unknown top-level/task
  keys rejected.

## Documentation + sample

- Command reference (en/ja): add `respack` to the targets list.
- A new section / tutorial example showing the full chain `input.yaml`
  (`tasks:` = scf, nscf, bands QE tasks + a `mode: respack` task) plus a RESPACK
  `input.in` template. Call out: SCDM (`N_initial_guess = 0`); energy windows /
  Ecut / flg_cRPA are user physics; the nscf task must set
  `nosym = .true.`/`noinv = .true.` (qe2respack); the documented run order
  scf → nscf → `qe2respack` → `respack_wannier`; the target is the
  `respack-wannier-py` port (SCDM); and `use_primitive: true` for a correct
  k-path.
- A `sample/cif2x/respack/<case>/` directory mirroring the existing samples
  (structure file + input.yaml + RESPACK template), if a runnable example is
  wanted (decide in the plan; may be docs-only inline like the sweep examples).

## Risks (documented)

- **Reciprocal-basis match:** `HighSymmKpath` coordinates are in the standardized
  primitive reciprocal basis, so the RESPACK path is only consistent with a
  primitive cell. respack targets `use_primitive: true` (and `use_ibrav: false`),
  same as the QE bands path; pymatgen's "path may be incorrect" warning is
  surfaced for non-standard cells.
- **Fortran-RESPACK incompatibility:** SCDM (`N_initial_guess = 0`, no Gaussian
  block) is supported by `respack-wannier-py`; the Fortran `calc_wannier` would
  reject it. Stated in the docs.
- **Directory/prefix alignment:** cif2x emits inputs only; the user must keep the
  QE `prefix`/`outdir` and the `qe2respack`/`respack_wannier` working directories
  consistent. The sample + docs show a consistent layout.

## Out of scope (v1)

- Manual Gaussian initial-guess projections (SCDM only).
- Spin-polarized split `input.in.up` / `input.in.dn` (LSDA) — follow-up.
- `qe2respack` invocation and any execution orchestration (cif2x only emits
  inputs; the run order scf→nscf→qe2respack→respack is documented, not executed).
- Auto-deriving energy windows / cRPA parameters (physics; template/content).
- Auto-injecting nscf `nosym`/`noinv` (documented in the sample for v1).

## Open items to resolve in the plan

- The exact mechanics of extracting the shared pp loader from
  `Struct2QE._set_pseudo_info` (which attributes it must keep populating on
  `Struct2QE` for backward compatibility) — QE behavior must stay unchanged
  (covered by existing QE tests).
- Whether `HighSymmKpath` segment labels (e.g. `\Gamma`) are emitted verbatim as
  inline comments or normalized to ASCII (comments are ignored by RESPACK, so
  either is acceptable; pick one in the plan).
