# RESPACK Tutorial (SrVO3) — Design

**Date:** 2026-07-01
**Status:** approved (pending written-spec review)

## Goal

Add a fully reproducible RESPACK tutorial to the cif2x manual (Japanese and
English), using SrVO3 as the example. A reader should be able to download the
required pseudopotentials, run `cif2x -t respack`, and obtain exactly the
`scf.in` / `nscf.in` / `input.in` files shown in the tutorial. The downstream
RESPACK run (QE SCF/NSCF → `qe2respack` → RESPACK) is presented as the "next
step" with the exact commands only — no computed physics results (e.g. cRPA U
values) are shown, because those are outside cif2x's scope and not CI-verifiable.

## Scope

**In scope**
- A new tutorial section "RESPACK ワークフロー" / "RESPACK workflow" appended to
  the existing cif2x tutorial pages (`docs/{ja,en}/source/cif2x/tutorial/index.rst`).
- A reproducible fixture set under `docs/tutorial/cif2x/respack/`.
- Committed cif2x-generated outputs (`scf/scf.in`, `nscf/nscf.in`, `input.in`),
  shown via `literalinclude`.
- Pseudopotential download instructions with URLs (ONCV / norm-conserving).

**Out of scope**
- Bundling pseudopotential UPF files in the repository (download instructions
  only — see "Pseudopotentials" below).
- Actual QE/`qe2respack`/RESPACK execution results. Only the run commands are
  given; the user verifies the physics on their own machine.
- Any change to cif2x source code. This is a documentation-only change. If
  generating the fixtures surfaces a real cif2x bug, that is handled separately.

## Pseudopotentials

- **Library:** ONCV (Optimized Norm-Conserving Vanderbilt), PBE. Norm-conserving
  pseudopotentials are the natural choice for RESPACK/cRPA because RESPACK
  consumes the Kohn–Sham wavefunctions directly.
- **Source (primary):** PseudoDojo — `http://www.pseudo-dojo.org/` (PBE,
  standard accuracy, UPF format; Sr, V, O all available with recommended
  cutoffs). Alternative: SG15 ONCV — `http://www.quantum-simulation.org/potentials/sg15_oncv/`.
- **Not bundled.** The repository ships download instructions and a `pp.csv`
  that names the exact UPF files. The reader downloads the same pseudopotential
  versions to reproduce the outputs.
- **Reproducibility note:** the committed `scf.in`/`nscf.in` carry cutoff values
  that depend on the exact UPF version. The tutorial pins the pseudopotential
  set (element, source, version, filename) precisely so the reader reproduces
  the same numbers. During implementation the UPFs are downloaded locally, cif2x
  is run to generate the outputs, and the exact filenames/cutoffs are recorded
  in `pp.csv` and (if the ONCV UPF header lacks suggested cutoffs) an explicit
  cutoff column / `content.system.ecutwfc` value.

## Fixture set (`docs/tutorial/cif2x/respack/`)

| File | Responsibility |
|------|----------------|
| `SrVO3.cif` | Standardized primitive cell of cubic-perovskite SrVO3 (Pm-3m, a ≈ 3.8425 Å; Sr (0,0,0), V (½,½,½), O at the three face centers). Single formula unit → primitive = conventional. |
| `input.yaml` | 3-task workflow (scf / nscf / respack). Based on `sample/cif2x/respack/SrVO3/input.yaml`, with `optional.pseudo_dir`/`pp_file` pointing at the ONCV set. `structure.use_primitive: true`, `use_ibrav: false`. |
| `respack.in_tmpl` | Namelist-only RESPACK template: `&param_chiqw` (cRPA on), `&param_wannier` (`N_wannier = 3` for the V-t2g manifold + energy windows), `&param_interpolation` (`dense`), `&param_calc_int`. |
| `pp.csv` | `element,pseudopotential,nexclude,orbitals` rows for Sr, V, O naming the ONCV UPF files. |
| `scf/scf.in` | Committed cif2x output for the scf task. |
| `nscf/nscf.in` | Committed cif2x output for the nscf task. |
| `input.in` | Committed cif2x output for the respack task (namelists + auto-generated k-path block). |
| `README.md` | Short pointer mirroring the `sample/` README, plus the pseudopotential download note. |

The existing `sample/cif2x/respack/SrVO3/` sample stays as-is (a minimal
pointer); the tutorial fixtures are the reproducible, output-carrying copy under
`docs/tutorial/`, consistent with how the QE tutorial uses `docs/tutorial/cif2x/`.

## Tutorial section content

Appended to `tutorial/index.rst` (both languages), same depth/idiom as the
existing QE and parameter-sweep sections. Subsections:

1. **Workflow overview.** Where cif2x sits in the RESPACK pipeline: QE SCF/NSCF
   → `qe2respack` → RESPACK. cif2x's job ends at generating the input files.
   Run order: scf → nscf → `qe2respack` → respack.
2. **Input parameter file.** `literalinclude` of `input.yaml`. Explain the
   3-task structure; why `structure.use_primitive: true` / `use_ibrav: false` is
   required (pymatgen high-symmetry k-path assumes the standardized primitive
   cell); why the nscf task sets `system.nosym: true` / `noinv: true` (required
   by `qe2respack`).
3. **RESPACK template.** `literalinclude` of `respack.in_tmpl`. Split what cif2x
   auto-fills (`&param_interpolation` k-path coordinates + `N_sym_points`; SCDM
   forced via `N_initial_guess = 0`) from what the user supplies as physics
   (`N_wannier`, `Lower/Upper_energy_window`, `dense`, cRPA settings). Note the
   template is namelist-only.
4. **Preparing pseudopotentials.** ONCV download URLs; how to build `pp.csv`;
   `optional.pseudo_dir` vs QE's runtime `control.pseudo_dir` (reuse the
   existing tutorial's explanation).
5. **Generating the inputs.** `cif2x -t respack input.yaml SrVO3.cif`;
   `literalinclude` of `scf/scf.in`, `nscf/nscf.in`, `input.in`. Point out the
   k-path block appended after `&param_interpolation`'s terminator in `input.in`.
6. **Next step (commands only).** The downstream run: `pw.x` scf → `pw.x` nscf →
   `qe2respack` → `RESPACK` (calc_wannier / calc_chiqw / calc_w3d / calc_j3d as
   appropriate). Commands only; no result output. Link to the RESPACK docs.

## Verification

- `cd docs/ja && make html` and `cd docs/en && make html` build without new
  warnings (all `literalinclude` paths resolve).
- The committed `input.in` k-path block matches live `cif2x -t respack` output
  for `SrVO3.cif` (regenerate and diff — must be identical).
- `scf.in` / `nscf.in` are the actual cif2x output for the pinned ONCV set.
- Full existing test suite still passes (no source change expected; run as a
  guard).

## Risks / open points

- **ONCV cutoffs:** if the PseudoDojo UPF headers don't expose cutoffs in the
  form cif2x parses, the tutorial sets `ecutwfc`/`ecutrho` explicitly in
  `input.yaml` (documented) rather than relying on auto-fill. Resolved at
  implementation when the UPFs are downloaded.
- **SrVO3 primitive cell:** must be given already standardized so the written
  k-path is consistent with the QE cell. The bundled `SrVO3.cif` is authored as
  the primitive cubic cell to satisfy this.
