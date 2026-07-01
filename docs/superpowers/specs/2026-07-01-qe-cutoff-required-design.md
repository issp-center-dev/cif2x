# QE cutoff: reject instead of silently defaulting to 0.0 — Design

**Date:** 2026-07-01
**Status:** approved (pending written-spec review)

## Problem

When cif2x generates a QE input and the plane-wave cutoffs (`ecutwfc`/`ecutrho`)
are left empty in `content`, it resolves them by, in order: the cutoff CSV
(`optional.cutoff_file`), then the UPF header / "Suggested" lines. If none of
these yields a value, `Struct2QE._find_elem_cutoff` silently substitutes `0.0`:

```python
if ecutwfc is None:
    ecutwfc = 0.0
if ecutrho is None:
    ecutrho = 0.0
```

and `_find_elem_cutoff_from_file` only logs `"cutoff information not found"`
(the original `raise` is commented out). The result is a generated `scf.in`
with `ecutwfc = 0.0` — a silently broken input that QE will reject downstream
with a confusing error, far from the real cause.

This bites norm-conserving (ONCV/PseudoDojo) pseudopotentials in particular:
their UPFs frequently do not carry `wfc_cutoff`/`rho_cutoff`, so the cutoff must
come from the CSV or an explicit `content` value — and if the user forgets,
today they get `0.0` with no hard failure.

## Goal

Fail loudly at generation time when a required cutoff cannot be resolved,
instead of emitting `ecutwfc = 0.0`. The error must name the offending
element(s) and tell the user how to fix it.

## Behavior

`_find_elem_cutoff(ename)` raises `InputValidationError` when, after consulting
the cutoff table and the UPF file, either `ecutwfc` or `ecutrho` is still
`None`. Message names the element and both remedies, e.g.:

> cutoff information not found for element 'Sr'. Provide it via
> `optional.cutoff_file` or set `content.system.ecutwfc`/`ecutrho` explicitly.

Unchanged (all still succeed exactly as before):
- Explicit `content.system.ecutwfc`/`ecutrho`: `_find_cutoff_info` is only
  called when those keys are empty, so explicit values never hit this path.
- Cutoff CSV hit (`_find_elem_cutoff_from_table`).
- UPF header / "Suggested" line hit (`_find_elem_cutoff_from_file`).

Only the previously-silent "nothing found → 0.0" case changes, and it changes
from a broken file to a clear error. No existing sample or test relies on the
`0.0` fallback (verified: every QE sample supplies `cutoff_file` or uses PSL
UPFs that carry header cutoffs).

## Changes

- `src/cif2x/struct2qe.py`
  - `_find_elem_cutoff`: replace the `0.0` substitution with a raised
    `InputValidationError` naming the element when either value is `None`.
  - `_find_elem_cutoff_from_file`: drop the dead commented-out `raise` /
    `logger.error("cutoff information not found")` at the end, since the caller
    now enforces the requirement. Keep the file-open / parse `logger.error`
    diagnostics that point at a specific file.
  - Import `InputValidationError` (from `cif2x.input_validator`) if not already.
- `tests/` — new test module (or add to an existing QE cutoff test file):
  - cutoff missing from both CSV and UPF → `InputValidationError` naming the
    element.
  - cutoff present in CSV → resolves, no error (regression guard).
  - cutoff present in UPF header → resolves, no error (regression guard).

## Verification

- New tests pass; full existing suite still green (no sample/test depended on
  the `0.0` fallback).
- Message names the element and both remedies.

## Out of scope

- The RESPACK tutorial (separate spec, `2026-07-01-respack-tutorial-design.md`)
  builds on this: it supplies ONCV cutoffs explicitly and relies on this clear
  error if a reader omits them.
- Any change to the cutoff-lookup order or the CSV/UPF parsing itself.
