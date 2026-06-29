# cif2x input.yaml validation — design

Date: 2026-06-30

## Goal

Replace cif2x's ad-hoc, traceback-heavy input handling (raw `KeyError` /
`RuntimeError` / library-internal `sys.exit`) with up-front schema validation of
`input.yaml` that produces clear, user-facing error messages. Validate the
target name and the required/known keys at startup, before any heavy work
(YAML parse of structure, pymatgen `Structure` construction).

Scope: **cif2x only**. `getcif` already validates its own input
(`QueryMaterialsProject` checks `properties`/`fields`), so it is out of scope.

## Decisions

- **Hand-rolled lightweight validator** (no `jsonschema`/`pydantic` dependency),
  so error messages are fully controlled.
- **Strict** unknown-key rejection — but only at the two *structural* levels
  (top-level keys, and per-task keys). The free-form blocks `optional:`,
  `structure:`, and `content:` are NOT key-checked (only type-checked as dicts),
  because they legitimately carry open-ended keys
  (`pseudo_dir`, `pp_file`, `cutoff_file`, `data_path`, `pseudo_map`, `tol_deg`, …).
- Validate the **raw task dict** (as written in YAML), *before* `deepupdate`
  injects the top-level `optional` block into each task — otherwise the injected
  nested `optional` key would be misclassified.
- **Preserve current required-key semantics** exactly (see table); only the
  messages and failure mode change.
- Output-file collision detection is **out of scope** (YAGNI): sweeps expand via
  `inflate` into subdirectories, so duplicate detection is error-prone and a hard
  error could break deliberate overwrite workflows.

## Per-target rules (data-driven)

A single table drives both validation and dispatch:

| canonical          | aliases                       | required task keys   | allowed task keys                                          |
|--------------------|-------------------------------|----------------------|-----------------------------------------------------------|
| `quantum_espresso` | qe, espresso, quantum_espresso| `mode`, `output_file`| mode, template, content, output_file, output_dir          |
| `vasp`             | vasp                          | (none)               | template, template_dir, content, output_file, output_dir  |
| `openmx`           | openmx                        | `output_file`        | mode, template, content, output_file, output_dir, precision|
| `akaikkr`          | akaikkr                       | `output_file`        | mode, template, content, output_file, output_dir, workdir |

`optional` is also accepted at task level (allowed everywhere) since per-task
optional overrides are conceptually supported by `deepupdate`.

This matches `main.py` today: QE requires `mode` + `output_file`; openmx/akaikkr
require only `output_file`; vasp requires neither. (Note: `Struct2QE` itself
treats a missing `mode` as `"generic"`, so the CLI contract is deliberately
stricter than the lower-level class — preserved as-is.)

## Components

### `src/cif2x/input_validator.py` (new, pure/testable)

- `class InputValidationError(Exception)` — carries a user-facing message.
- `TARGETS` — the data-driven table above (aliases, required/allowed keys; the
  `Struct2X` class is wired in `main.py` to avoid import cycles, or referenced
  here by name).
- `normalize_target(target: str) -> str` — canonicalize aliases
  (case-insensitive); raise `InputValidationError` listing valid targets on
  unknown input.
- `validate_input(info_dict, target_canonical) -> None` — raise
  `InputValidationError` on the first violation:
  1. `info_dict` is a non-None dict.
  2. unknown top-level keys (allowed: `structure`, `optional`, `tasks`) → error.
  3. `structure`/`optional` present → must be dicts (contents not key-checked).
  4. `tasks` present, is a non-empty list; each item is a dict.
  5. per task (raw dict): required keys present; unknown keys (vs per-target
     allowed set + `optional`) → error; `output_file`/`output_dir` if present are
     strings.

### `src/cif2x/main.py` (modified)

Reordered flow with friendly failure:

```
normalize_target(args.target)            # fail fast on bad target, no file I/O
load YAML (friendly message on parse err)
validate_input(info_dict, target)
struct = Cif2Struct(cif_file, info_dict.get("structure", {}))
for task in tasks: params = {**task}; deepupdate(params, {"optional": optional})
    TARGETS[target].cls(params, struct).write_input(output_file, output_dir, dry_run)
```

- `main()` wraps the body in `try/except InputValidationError` → `logger.error`
  + `sys.exit(1)`, **no traceback**.
- The four near-identical target branches collapse into one loop keyed off the
  `TARGETS` table.
- Inline `mode`/`output_file` checks are removed (now centralized in the
  validator).
- Structure-level validation (e.g. `cell_shape` required for molecular `.xyz`
  inputs) stays in `Cif2Struct`; not duplicated here.

## Error message examples

```
Error: unsupported target 'qe2'. Choose from: quantum_espresso (qe, espresso), vasp, openmx, akaikkr.
Error: task 1: 'output_file' is required for target 'quantum_espresso'.
Error: task 2: unknown key 'output_fil' for target 'vasp'. Allowed: template, template_dir, content, output_file, output_dir.
Error: 'tasks' is empty; nothing to generate.
```

## Testing

`tests/test_input_validator.py` (pytest):

- `normalize_target`: each alias canonicalizes; unknown target → `InputValidationError`
  with the candidate list in the message.
- `validate_input` raises for: missing `mode` (QE), missing `output_file`,
  unknown top-level key, unknown task key, `tasks` not a list, `tasks`
  empty/absent, empty (None) input, task entry not a dict.
- `validate_input` passes when `optional`/`structure`/`content` carry free-form
  keys (no false positives).
- Regression: every `sample/cif2x/<target>/*/input.yaml` validates successfully
  under its directory's target (collected via glob).

No linter/test runner is configured in the repo; tests run with `pytest` (added
to dev-dependencies).

## Out of scope

- getcif validation (already validated).
- Output-file collision detection.
- Strict validation of `structure:`/`optional:`/`content:` block contents.
- GitHub Actions CI to run pytest (separate improvement item).
