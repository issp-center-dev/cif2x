# cif2x parameter-sweep examples for VASP / OpenMX / AkaiKKR — design

Date: 2026-06-30

## Goal

The `${...}` parameter-sweep syntax already works for every cif2x target, but the
tutorial only demonstrates it for Quantum ESPRESSO. Add per-target sweep examples
to the tutorial (English + Japanese) and a regression test that the sweep
expansion works through each target's content class — so the documented
behaviour is real and stays working.

Scope: documentation + tests only. No production-code changes (the sweep already
works for all four targets).

## Background (verified)

- `inflate(content)` lives in `src/cif2x/utils.py` and is shared: every generator
  calls it — `struct2qe.py:63`, `struct2vasp.py:74`, `struct2openmx.py:142`,
  `struct2akaikkr.py:63`. It returns a list of `(dirkey, content)`; with no
  `${...}` it returns `[("", content)]`.
- `inflate` is content-agnostic: it uses `content.serialize()` /
  `type(content).deserialize(...)`. Each target ships a small content class with
  those methods (QE `cif2x/qe/content.py`; VASP/OpenMX/AkaiKKR each define a
  `Content` in their `struct2*.py`).
- The directory key is derived from the swept value (`_to_string`): scalars →
  `str(value)` (e.g. `400`), lists → joined with `x` (e.g. `4x4x4`).
- The `${...}` syntax is already documented generically in
  `docs/{en,ja}/source/cif2x/filespec/index.rst`; the tutorial
  (`docs/{en,ja}/source/cif2x/tutorial/index.rst`, "Specifying parameter sets")
  shows only a QE `K_POINTS.grid` example.
- `tests/test_inflate.py` already covers the generic `utils.inflate` (list,
  range, no-placeholder, cartesian product) via a `DictContent` stand-in and the
  QE `Content`. It does NOT exercise the VASP/OpenMX/AkaiKKR content classes.
- `bzqlty` is a valid AkaiKKR content key (`akaikkr/read_input.py:199`,
  `make_input.py:140`) — a natural Brillouin-zone-quality convergence knob.

## Decisions

- Examples live inline in the existing tutorial sweep section (no new sample
  directories — those need pseudopotentials/structure files to be runnable).
- One example per target, each a convergence-style sweep mirroring the QE one:
  - VASP: `incar.ENCUT` → `${[400, 600, 800]}`
  - OpenMX: `scf.energycutoff` → `${[150, 200, 250]}`
  - AkaiKKR: `bzqlty` → `${[12, 16, 20]}`
- Regression test extends `tests/test_inflate.py` by running `inflate` through
  each target's real `Content` class (not just the generic stand-in).

## Components

### Documentation (`docs/en` and `docs/ja` `cif2x/tutorial/index.rst`)

In the "Specifying parameter sets" / corresponding JA section, after the existing
QE example, add a sentence that the sweep applies to all targets (the output
sub-directory name comes from the swept value), then three small `code-block:: yaml`
examples:

```yaml
# VASP — cutoff convergence
content:
  incar:
    ENCUT: ${ [400, 600, 800] }      # -> 400/, 600/, 800/
```
```yaml
# OpenMX — energy-cutoff convergence
content:
  scf.energycutoff: ${ [150, 200, 250] }   # -> 150/, 200/, 250/
```
```yaml
# AkaiKKR — Brillouin-zone-quality convergence
content:
  bzqlty: ${ [12, 16, 20] }          # -> 12/, 16/, 20/
```

EN and JA stay consistent. (`filespec` is unchanged — the syntax is already
generic there.)

### Regression test (`tests/test_inflate.py`)

Add cases that build each target's content class and run the shared `inflate`:

- VASP: `from cif2x.struct2vasp import Content as VaspContent`;
  `VaspContent(incar={"ENCUT": "${[400, 600]}"}, kpoints={}, poscar={}, potcar={})`
  → `inflate` yields 2 entries; the two `incar["ENCUT"]` values are `400` and
  `600`; dir keys are `"400"`/`"600"`.
- AkaiKKR: `from cif2x.struct2akaikkr import Content as KkrContent`;
  `KkrContent(go="go", bzqlty="${[12, 16]}")` → 2 entries; `bzqlty` values
  `12`/`16`.
- OpenMX: `from cif2x.struct2openmx import Content as OpenmxContent` if it is
  constructible from a plain mapping (`Content.from_dict({"scf.energycutoff":
  "${[150, 200]}"})`); assert 2 entries with the substituted values. If its
  constructor needs target-specific objects that make a unit test brittle, omit
  the OpenMX case and rely on VASP + AkaiKKR + the existing QE/DictContent
  coverage (note this in the test/commit). 

Use `cif2x.utils.inflate` (the shared function the generators call). Imports that
pull heavy optional deps should use `pytest.importorskip` consistent with the
existing QE tests (e.g. skip if `qe_tools`/`pymatgen` import side effects block
loading a target module).

## Testing / verification

- `pytest tests/test_inflate.py` — existing + new cases pass.
- Manual: render one target's swept content (e.g. run `inflate` on the VASP
  content above) and confirm the directory keys are `400`/`600`/`800`.
- `docutils` lint on the two edited tutorial files (no severe/error).

## Out of scope

- New runnable sample directories (would require pseudopotentials/structures).
- Any change to the sweep mechanism itself or to `filespec` syntax docs.
- Per-target sweep support in code (already works).
