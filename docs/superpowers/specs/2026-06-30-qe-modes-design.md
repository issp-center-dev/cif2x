# cif2x QE mode expansion (relax / vc-relax / bands; dos via docs) — design

Date: 2026-06-30

## Goal

Extend the cif2x Quantum ESPRESSO generator beyond `scf`/`nscf` so that
`relax`, `vc-relax`, and `bands` produce correct pw.x inputs with
auto-generated structure cards. `dos` is covered by a documented task chain
(no new code). The headline new capability is a band-path `K_POINTS {crystal_b}`
generator built from pymatgen's `HighSymmKpath`.

Scope: the QE target only (`src/cif2x/qe/`).

## Background (verified against the code)

- `calc_mode.create_modeproc(mode, qe)` routes `scf`/`nscf` → `QEmode_pw`, all
  others → `QEmode_generic` (template passthrough, no structure cards).
- `QEmode_pw.card_table` maps `CELL_PARAMETERS`, `ATOMIC_SPECIES`,
  `ATOMIC_POSITIONS`, `K_POINTS` to generators. `update_namelist` sets
  `calculation = self.qe.mode` and updates struct/cutoff/nspin/nbnd/control
  blocks generically (not scf-specific).
- `generate_k_points(qe, params)` supports `option` ∈ {`gamma`, `automatic`,
  `crystal`} (uniform meshes); unknown option currently returns `{}` (drops the
  card). The `crystal` integer is a per-point weight, NOT a segment length.
- `content.py` card renderer: a scalar `data` item renders as one line; a list
  item renders as its elements joined by two spaces; floats as `%.6f`.
- `HighSymmKpath(structure).kpath` → `{"kpoints": {label: frac_coords},
  "path": [[label, ...], ...]}` (segments may be disjoint; labels are LaTeX-ish,
  e.g. `"\\Gamma"`).

## Decisions

- **Routing:** add `relax`, `vc-relax`, `bands` to the `QEmode_pw` list in
  `create_modeproc`. No new namelist logic — `calculation` is set to the mode
  string automatically; `&ions`/`&cell` come from the user's template/content.
- **bands K_POINTS:** new option `crystal_b` (only — `tpiba_b` is deferred,
  YAGNI). Built from `HighSymmKpath`. Numeric-only rows (no inline labels).
- **Raw passthrough escape hatch:** `generate_k_points` no longer drops unknown
  options or explicit user data — a verbatim `K_POINTS` block survives
  `QEmode_pw`.
- **bands guard:** when `qe.mode == "bands"`, a non-band K_POINTS option
  (e.g. the default `automatic`) is rejected with a clear error, instead of
  silently emitting wrong physics.
- **dos:** no code. pw.x has no `calculation='dos'`; the DOS pw step is `nscf`
  (already supported) and dos.x input (`&DOS`, no cards) works via
  `QEmode_generic`. Covered by a sample task chain + docs.

## Components

### `src/cif2x/qe/calc_mode.py`

`create_modeproc`: route list becomes
`["scf", "nscf", "relax", "vc-relax", "bands"]` → `QEmode_pw`; else
`QEmode_generic`. No other change. (`_update_nbnd_info` still auto-fills `nbnd`
only when the template leaves it empty — documented caveat for bands, see
Limitations.)

### `src/cif2x/qe/cards.py` — `generate_k_points`

1. **Raw passthrough first:** if `option` has no generator (unknown) OR the
   caller supplied explicit `data` for a `_b` option, return the card unchanged
   (`{'key': 'K_POINTS', 'option': option, 'data': <user data>}`) instead of
   `{}`. This preserves a hand-written `crystal_b` block from a template.
2. **`crystal_b` branch** (new):
   - Resolve the path: if `params["path"]` (a list of label sequences) is given,
     resolve each label against `HighSymmKpath(qe.struct.structure).kpath["kpoints"]`;
     a missing label raises a clear error listing the available labels. Otherwise
     use the auto `kpath["path"]`.
   - `line_npoints` (int, default 20): the QE per-line integer = number of points
     generated from this k-point to the next.
   - Emit, for each segment, every high-symmetry point as a row
     `[kx, ky, kz, line_npoints]`; **the terminal point of each disjoint segment
     (and the final point overall) carries `nk = 0`** so bands.x does not
     interpolate across a discontinuity. (The exact crystal_b break convention is
     cross-checked against QE `INPUT_PW` during implementation.)
   - `data = [[N]] + rows`, where `N` is the total number of emitted rows.
   - `option` is `crystal_b`.
3. **bands guard:** the caller (QEmode_pw, knowing `qe.mode`) rejects a
   `bands` task whose resolved K_POINTS option is not `crystal_b`/`tpiba_b`.
   (Implementation detail: the guard lives where `qe.mode` is in scope — either a
   check in `QEmode_pw.update_cards` or passing `qe.mode` into the K_POINTS
   generator. Decide in the plan; behaviour is "clear error, no automatic mesh
   for bands".)

### dos / samples

- `sample/cif2x/qe/<dos-case>/input.yaml`: a `tasks:` chain scf → nscf (dense
  mesh) → dos (`mode: dos`, `&DOS` namelist in template/content).

## Error handling

- bands without a band K_POINTS option → clear error (no silent automatic mesh).
- unknown high-symmetry label in a `path` override → error listing available
  labels.
- Unknown/explicit K_POINTS data → passed through (not an error).

## Testing (network-free, deterministic)

`tests/test_qe_calc_mode.py` / `tests/test_qe_cards.py`:

- **Routing:** `create_modeproc("relax", qe)`, `"vc-relax"`, `"bands"` each
  return a `QEmode_pw` instance; `"dos"` returns `QEmode_generic`.
- **relax/vc-relax:** on a small `use_ibrav=False` structure, `update_namelist`
  sets `calculation` to the mode, and `update_cards` yields CELL_PARAMETERS,
  ATOMIC_SPECIES, ATOMIC_POSITIONS, and an `automatic` K_POINTS card.
- **bands crystal_b:** fixed small structure; assert the card has
  `option == "crystal_b"`, a leading count row, each row `[kx,ky,kz,nk]` with
  `nk == line_npoints` for interior points and `nk == 0` at each disjoint-segment
  terminal (the highest-risk assertion); compare coords at rounded precision.
- **bands guard:** a `bands` task with `automatic`/no option → error.
- **raw passthrough:** an unknown option, and an explicit-`data` `crystal_b`
  block, both survive `generate_k_points` unchanged.
- **path override:** a bogus label → error whose message contains the available
  label set.

## Documentation

- QE tutorial / basic-usage (en + ja) and samples: show `bands`
  (`option: crystal_b`, `line_npoints`), `relax`/`vc-relax`, and the dos chain.
- State caveats: for `bands`, the template must set `nbnd` explicitly if
  conduction bands are wanted (the generator only fills `nbnd` when left empty,
  to a valence-oriented default); `dos` is a chained workflow, not a pw.x
  calculation.

## Out of scope / limitations (v1)

- `tpiba_b` (only `crystal_b`).
- Inline high-symmetry labels in the crystal_b card (numeric only).
- `use_ibrav: true` combined with `bands` (bands targets `use_ibrav=False`;
  ibrav+bands is a documented limitation pending a separate verification of the
  `system`/CELL_PARAMETERS path).
- Multi-executable orchestration for DOS (cif2x emits inputs; the user runs the
  chain).
- Per-segment `line_npoints` (single global value for v1; the param name leaves
  room to extend later).
