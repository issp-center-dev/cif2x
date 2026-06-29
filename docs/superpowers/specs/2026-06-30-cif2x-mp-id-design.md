# cif2x `--mp-id` (fetch structure from Materials Project) — design

Date: 2026-06-30

## Goal

Let `cif2x` build inputs directly from a Materials Project material id, turning
the current two-step "`getcif` to download a CIF, then `cif2x` to generate
inputs" into one command:

```
cif2x -t qe --mp-id mp-149 input.yaml
```

The positional CIF file becomes optional; exactly one of {CIF file, `--mp-id`}
must be supplied.

Scope: `cif2x` only (one material id per run). The fetched structure is written
to a transient temporary CIF and read through the existing `Cif2Struct`, so the
result is identical to running `getcif` then `cif2x`.

## Decisions

- **Temp-CIF round-trip** (not passing a `Structure` object): fetch the MP
  structure, write it to a temp CIF applying `symprec`, then hand the path to the
  unmodified `Cif2Struct`. This is faithful to the existing `getcif` to `cif2x`
  flow and keeps `Cif2Struct` unchanged.
- **Do not persist** the downloaded CIF (temporary directory, auto-cleaned).
- **Reuse getcif's API-key + symprec semantics** via a shared helper (single
  source of truth); `symprec` is a CLI option defaulting to 0.1, and
  `symprec == 0` disables symmetry refinement (mirrors getcif).
- All MP fetch failures are converted to `InputValidationError` so `main()`
  exits cleanly (status 1, no traceback) — they all occur before `Cif2Struct`
  is constructed.
- `mp_api` is imported lazily so `cif2x` keeps working (and its tests keep
  importing) when only local CIFs are used.

## Components

### `src/getcif/mp.py` (new, shared)

Pure-ish MP helpers, usable by both tools. Top-level imports avoid `mp_api`
(lazy inside `fetch_structure`).

- `resolve_api_key(api_key_file="materials_project.key") -> str | None`
  The key-resolution logic currently inline in `_setup_dbinfo`, extracted
  **verbatim** so getcif behavior does not drift: when `api_key_file` ends in
  `.key` and exists, take the first line that does not start with `#` after
  `strip()`; otherwise `None` (so `MPRester` falls back to env var / pymatgen
  settings). This preserves the existing edge case where a leading blank line
  yields `""` (treated as "no key" → fallback); changing that is out of scope.
- `_import_mprester() -> type` — module-level helper that does the lazy
  `from mp_api.client import MPRester` and returns the class. This is the
  **patchable seam**: tests monkeypatch `getcif.mp._import_mprester` to inject a
  fake rester, so `mp_api` is never imported (or required) in unit tests.
- `fetch_structure(material_id, *, api_key=None) -> Structure`
  `MPRester = _import_mprester()`; inside
  `with MPRester(api_key=api_key, mute_progress_bars=True) as mpr:` call
  `mpr.materials.summary.search(material_ids=[material_id], fields=["structure"])`.
  Return `docs[0].structure`; raise `LookupError` if `docs` is empty. (Returns
  the MP final stored structure — not a conventionalized cell.)

`getcif/main.py` `_setup_dbinfo` is refactored to
`self.api_key = resolve_api_key(api_key_file)` (single source of truth). The
existing eager `mp_api` import in `getcif/main.py` is unchanged.

### `src/cif2x/mp_source.py` (new)

- `fetch_to_cif(material_id, dest_path, *, symprec=0.1, api_key_file="materials_project.key") -> None`
  - `api_key = resolve_api_key(api_key_file)`
  - The whole fetch-and-write body is wrapped so **every** failure (auth,
    network, and `structure.to(...)` serialization) becomes `InputValidationError`
    (exit 1, no traceback), not just the fetch call:
    - `LookupError` (empty docs) → `"material '<id>' not found in the Materials Project."`
    - any other exception → `"failed to fetch '<id>' from the Materials Project: <e>"`
  - on success: `structure.to(dest_path, fmt="cif", symprec=(symprec or None))`

Imports `InputValidationError` from `cif2x.input_validator` (same package).

### `src/cif2x/main.py` (modified)

CLI:
- `cif_file` positional becomes `nargs="?"` (optional); `input_file` stays
  required.
- add `--mp-id`, `--symprec` (type float, default 0.1), `--api-key-file`
  (default `materials_project.key`). `--symprec` / `--api-key-file` are used only
  when `--mp-id` is set; with a positional CIF they are ignored (no warning).
- Because `cif_file` is now optional, supplying neither source is no longer an
  argparse "required argument" error but the `InputValidationError` below
  (`provide a CIF file or --mp-id.`); tests assert this new message.

`_run` flow:
```
target = normalize_target(args.target)
load + validate YAML
struct_params = info_dict.get("structure") or {}
_require_one_source(args)            # exactly one of cif_file / mp_id
if args.mp_id:
    with tempfile.TemporaryDirectory() as td:
        cif_path = Path(td) / "structure.cif"
        fetch_to_cif(args.mp_id, cif_path, symprec=args.symprec,
                     api_key_file=args.api_key_file)
        struct = Cif2Struct(str(cif_path), struct_params)
else:
    struct = Cif2Struct(args.cif_file, struct_params)
# dispatch loop unchanged
```

`_require_one_source(args)` raises `InputValidationError`:
- both set → `"provide either a CIF file or --mp-id, not both."`
- neither → `"provide a CIF file or --mp-id."`

`Cif2Struct` is built inside the `TemporaryDirectory` block; it reads the file
synchronously in `__init__` and keeps the `Structure` in memory, so the temp dir
can be removed afterward.

## Error handling

| Situation | Message (exit 1, no traceback) |
|-----------|--------------------------------|
| both sources | `provide either a CIF file or --mp-id, not both.` |
| no source | `provide a CIF file or --mp-id.` |
| unknown/deprecated id (empty docs) | `material 'mp-x' not found in the Materials Project.` |
| auth / network / other MP error | `failed to fetch 'mp-x' from the Materials Project: <e>` |

A missing API key is not pre-rejected: `resolve_api_key` returns `None` and
`MPRester(api_key=None)` falls back to env/pymatgen settings (parity with
getcif). A genuine 401 surfaces through the generic fetch-failure message.

## Testing (network-free)

`tests/test_mp_source.py` / `tests/test_main_cli.py`:
- `resolve_api_key`: first non-`#` line returned; missing file → `None`;
  non-`.key` filename → `None`.
- `fetch_structure`: monkeypatch `getcif.mp._import_mprester` to return a fake
  `MPRester` class (a context manager whose `materials.summary.search` returns a
  one-element list); assert it is called with `material_ids=[id]`,
  `fields=["structure"]`; empty list → `LookupError`. (mp_api is never imported.)
- `fetch_to_cif`: monkeypatch `fetch_structure` to return a small pymatgen
  `Structure`; assert a CIF is written at `dest_path`; `LookupError` and generic
  exceptions both become `InputValidationError` with the id in the message.
- CLI mutual-exclusion: `cif_file` + `--mp-id` → `InputValidationError`; neither
  → `InputValidationError`.
- E2E: monkeypatch `main_mod.fetch_to_cif` (writes a known CIF) and
  `Cif2Struct`/the target class; assert the `--mp-id` path reaches
  `write_input`.

## Documentation

- Add a `--mp-id` example to the cif2x tutorial / basic-usage, alongside the
  existing positional-CIF form.
- Document the **known limitations**, both shared with the existing
  getcif → cif2x flow:
  - CIF serialization can drop site properties such as `magmom`.
  - MP disordered/partial-occupancy structures are flagged `is_composite` and
    rejected by the QE generator.
  - `symprec == 0` means no symmetry refinement.

## Out of scope

- Multiple material ids in one run (cif2x builds exactly one structure).
- Persisting the downloaded CIF (`getcif` already does that).
- Passing a `Structure` object straight into `Cif2Struct` (would change its
  interface and skip symprec).
- A `database:` block in cif2x's `input.yaml` (API key comes from
  `--api-key-file` / env).
