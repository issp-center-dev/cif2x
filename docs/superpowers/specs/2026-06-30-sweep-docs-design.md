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

- There are **two near-identical `inflate` implementations**: VASP, OpenMX, and
  AkaiKKR use `cif2x.utils.inflate` (each does `from cif2x.utils import *` and
  calls `inflate(self.content)` — `struct2vasp.py:74`, `struct2openmx.py:142`,
  `struct2akaikkr.py:63`); QE uses its own `cif2x.qe.content.inflate`
  (`struct2qe.py:15,63`). Both return a list of `(dirkey, content)` and, with no
  `${...}`, return `[("", content)]`. This work tests the three non-QE targets
  through `cif2x.utils.inflate`; QE's path is already covered by existing
  `tests/test_inflate.py` cases against `cif2x.qe.content.inflate`. (Unifying the
  two implementations is out of scope; see Risks.)
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
QE example, add the exact bridging text below, then a note about where the swept
key sits in `content` for each target (VASP groups parameters under
`incar`/`kpoints`/`poscar`/`potcar`; OpenMX and AkaiKKR use flat `content` keys),
then three small `code-block:: yaml` examples.

Exact text (EN):

> The same ``${...}`` syntax applies to every target, not only Quantum ESPRESSO;
> the output sub-directory name is derived from the swept value. The swept key
> sits wherever the parameter lives in ``content`` — for VASP that is under
> ``incar`` (or ``kpoints``/``poscar``/``potcar``); for OpenMX and AkaiKKR the
> ``content`` keys are flat. For example:

Exact text (JA):

> 同じ ``${...}`` 構文は Quantum ESPRESSO だけでなく全ターゲットで利用でき、
> 出力サブディレクトリ名はスイープ値から生成されます。スイープ対象のキーは
> ``content`` 内でそのパラメータが定義される階層に記述します。VASP では ``incar``
> (または ``kpoints``/``poscar``/``potcar``)の下、OpenMX と AkaiKKR では
> ``content`` 直下のフラットなキーです。例:

Each example's ``content:`` block sits under an entry of ``tasks:``, exactly like
the existing QE example earlier in the section. Exact one-liners to append:

> EN: "As in the earlier example, each ``content:`` block below belongs to a
> ``tasks:`` entry."
> JA: 「以下の各 ``content:`` ブロックは、前述の例と同様に ``tasks:`` の各エントリ
> 配下に置かれます。」

Both the spaced (``${ [400, 600, 800] }``) and compact (``${[400,600,800]}``) forms
are accepted — whitespace inside ``${...}`` is ignored (verified) — so the tutorial
may use the spaced form (matching the existing QE example) while the tests use the
compact form. The unquoted form parses as a plain string under the safe YAML
loader (verified for ``grid``/``ENCUT``/``bzqlty``), so no surrounding quotes are
needed — keep the examples unquoted, consistent with the existing QE example.

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

Add **one required case per non-QE target**, each building that target's real
content class and running `cif2x.utils.inflate` (the function those targets call).
The content classes are plain dict-wrappers with **no init-time validation and no
filesystem access** (verified), so minimal construction is safe. Exact
constructor signatures (verified):

- VASP — `cif2x.struct2vasp.Content(incar, kpoints, poscar, potcar, **kwargs)`:
  `Content(incar={"ENCUT": "${[400, 600]}"}, kpoints={}, poscar={}, potcar={})`
  → 2 entries; the `(dirkey, incar["ENCUT"])` pairs are exactly
  `("400", 400)` and `("600", 600)`.
- AkaiKKR — `cif2x.struct2akaikkr.Content(**kwargs)`:
  `Content(go="go", bzqlty="${[12, 16]}")` → 2 entries; `(dirkey, bzqlty)` pairs
  `("12", 12)` and `("16", 16)`.
- OpenMX — `cif2x.struct2openmx.Content` (a case-insensitive dict subclass) built
  via `Content.from_dict({"scf.energycutoff": "${[150, 200]}"})` →
  2 entries; `(dirkey, scf.energycutoff)` pairs `("150", 150)` and `("200", 200)`.
  This case is **required** (the class is constructible from a plain mapping, so
  there is no reason to omit it).

Field accessors (verified): VASP/AkaiKKR `Content.__getitem__` reads the wrapped
dict (`content["incar"]["ENCUT"]`, `content["bzqlty"]`); OpenMX `Content` is a
dict subclass (`content["scf.energycutoff"]`).

Test requirements:
- Each target case does `mod = pytest.importorskip("cif2x.struct2vasp")` (and
  `..struct2openmx`, `..struct2akaikkr`) **inside the test body**, then uses
  `mod.Content`. Guarding on the target module — not a specific library name —
  skips only that case when ANY of its transitive optional deps (pymatgen, monty,
  pandas, …) is missing, without affecting the generic/QE cases. (Tests run under
  the existing `tests/conftest.py`, which puts `src/` on `sys.path`.)
- Assert the `(dirkey, value)` pairs **together** (sorted list of tuples), not the
  keys and values separately, so a mismatch between the generated sub-directory
  name and the substituted value is caught.
- Also add the paired `(dirkey, value)` assertion to the existing QE inflate
  case(s) (currently value-only), so the separate `cif2x.qe.content.inflate` path
  is validated the same way. Keep the generic (`DictContent`) cases unchanged.
  (Verified: both `cif2x.utils.inflate` and `cif2x.qe.content.inflate` return the
  same `list[(dirkey, content)]` shape, so the paired-assertion helper applies to
  both.)

## Testing / verification

- `pytest tests/test_inflate.py` — existing + new cases pass.
- Manual: render one target's swept content (e.g. run `inflate` on the VASP
  content above) and confirm the directory keys are `400`/`600`/`800`.
- Best-effort reST lint on the two edited tutorial files with the exact command
  `python3 -m docutils <file> /dev/null` (run only if `docutils` is importable),
  expecting no `severe`/`error` lines. CI is expected to run in an environment
  where all three target modules import, so the new per-case `importorskip`
  guards skip nothing there (they only protect minimal local environments, and a
  target module's import error would surface in CI rather than be masked).

## Risks

- QE keeps a separate `cif2x.qe.content.inflate`, so sweep behaviour could
  diverge from the other three targets even with the new tests green. Unifying
  the two implementations is deliberately out of scope here; the divergence risk
  is noted for a future cleanup.
- Tests build real target content classes; a missing optional dependency
  (e.g. pymatgen) is handled by per-case `pytest.importorskip` so it skips only
  that case. (Verified: the content classes do no init-time validation or
  filesystem I/O, so construction itself is side-effect-free.)

## Out of scope

- New runnable sample directories (would require pseudopotentials/structures).
- Any change to the sweep mechanism itself or to `filespec` syntax docs.
- Per-target sweep support in code (already works).
