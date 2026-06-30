# cif2x parameter-sweep examples (VASP/OpenMX/AkaiKKR) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Document the `${...}` parameter sweep for VASP, OpenMX, and AkaiKKR in the tutorial (EN+JA) and add a regression test that the sweep expands correctly through each target's content class.

**Architecture:** Documentation + tests only — the sweep already works for all targets. Tests exercise the shared `cif2x.utils.inflate` (used by VASP/OpenMX/AkaiKKR) through each target's real `Content` class; the existing QE cases (separate `cif2x.qe.content.inflate`) gain paired `(dirkey, value)` assertions.

**Tech Stack:** Python, pytest, Sphinx reStructuredText. Tests run from repo root (`tests/conftest.py` puts `src/` on `sys.path`).

---

## Background facts (verified by running them)

- `cif2x.utils.inflate(content)` returns `list[(dirkey, content)]`; with no `${...}` it returns `[("", content)]`. Used by VASP/OpenMX/AkaiKKR. QE uses the near-identical `cif2x.qe.content.inflate`; both return the same shape.
- Verified end-to-end results (these exact tuples are what the tests assert):
  - VASP `Content(incar={"ENCUT": "${[400, 600]}"}, kpoints={}, poscar={}, potcar={})` → `[("400", 400), ("600", 600)]` via `c["incar"]["ENCUT"]`.
  - AkaiKKR `Content(go="go", bzqlty="${[12, 16]}")` → `[("12", 12), ("16", 16)]` via `c["bzqlty"]`.
  - OpenMX `Content.from_dict({"scf.energycutoff": "${[150, 200]}"})` → `[("150", 150), ("200", 200)]` via `c["scf.energycutoff"]`.
- The content classes do no init-time validation or filesystem I/O; minimal construction is safe.
- Unquoted `${ [...] }` parses as a plain string under the safe YAML loader (verified for `grid`/`ENCUT`/`bzqlty`); whitespace inside `${...}` is ignored. Keep tutorial examples unquoted (matches the existing QE example).
- `tests/test_inflate.py` already imports `from cif2x.utils import inflate as utils_inflate` and `from cif2x.qe.content import Content, inflate as qe_inflate`; it does NOT yet import `pytest`. Its QE cases assert values only.
- EN tutorial: the QE sweep example ends with the sentence beginning "When ``K_POINTS`` is given as above, ..." in `docs/en/source/cif2x/tutorial/index.rst` (section "Specifying parameter sets"). JA tutorial: section "パラメータセットを指定する", the QE `grid:` example at line ~76 followed by an explanatory sentence, in `docs/ja/source/cif2x/tutorial/index.rst`.

## File structure

- **Modify** `tests/test_inflate.py` — add `import pytest`; add 3 target cases; add paired assertions to the 2 existing QE cases.
- **Modify** `docs/en/source/cif2x/tutorial/index.rst` — per-target examples after the QE example.
- **Modify** `docs/ja/source/cif2x/tutorial/index.rst` — same in Japanese.

---

## Task 1: Regression test — sweep through each target's content class

**Files:**
- Modify: `tests/test_inflate.py`

- [ ] **Step 1: Add `import pytest` and the three target cases (failing)**

At the top of `tests/test_inflate.py`, add `import pytest` (after `import json`). Then append:

```python
def test_vasp_content_inflate_expands_sweep():
    mod = pytest.importorskip("cif2x.struct2vasp")
    content = mod.Content(incar={"ENCUT": "${[400, 600]}"},
                          kpoints={}, poscar={}, potcar={})
    results = utils_inflate(content)
    assert sorted((k, c["incar"]["ENCUT"]) for k, c in results) == \
        [("400", 400), ("600", 600)]


def test_akaikkr_content_inflate_expands_sweep():
    mod = pytest.importorskip("cif2x.struct2akaikkr")
    content = mod.Content(go="go", bzqlty="${[12, 16]}")
    results = utils_inflate(content)
    assert sorted((k, c["bzqlty"]) for k, c in results) == \
        [("12", 12), ("16", 16)]


def test_openmx_content_inflate_expands_sweep():
    mod = pytest.importorskip("cif2x.struct2openmx")
    content = mod.Content.from_dict({"scf.energycutoff": "${[150, 200]}"})
    results = utils_inflate(content)
    assert sorted((k, c["scf.energycutoff"]) for k, c in results) == \
        [("150", 150), ("200", 200)]
```

- [ ] **Step 2: Run to verify current state**

Run: `pytest tests/test_inflate.py -q`
Expected: the three new tests are collected and run (they should PASS immediately because the sweep already works — this test locks that behavior). If any FAILS, inspect whether the target's `Content` API changed; do not weaken the assertion without checking the actual `(dirkey, value)` output.

- [ ] **Step 3: Strengthen the existing QE cases with paired assertions**

In `tests/test_inflate.py`, in `test_qe_inflate_expands_list_placeholder`, after the existing
`assert sorted(cc.namelist["system"]["ecutwfc"] for _, cc in results) == [30, 40]`, add:

```python
    assert sorted((k, cc.namelist["system"]["ecutwfc"]) for k, cc in results) == \
        [("30", 30), ("40", 40)]
```

In `test_qe_inflate_expands_range_placeholder`, after its existing value assertion
(`... == [30, 40]`), add the same paired assertion:

```python
    assert sorted((k, cc.namelist["system"]["ecutwfc"]) for k, cc in results) == \
        [("30", 30), ("40", 40)]
```

- [ ] **Step 4: Run the test file**

Run: `pytest tests/test_inflate.py -q`
Expected: PASS (existing generic + QE cases with the new paired assertions + the 3 new target cases).

- [ ] **Step 5: Commit**

```bash
git add tests/test_inflate.py
git commit -m "test: cover parameter sweeps for VASP/OpenMX/AkaiKKR content classes"
```
End the commit body with a blank line then:
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

---

## Task 2: Tutorial documentation (EN + JA)

**Files:**
- Modify: `docs/en/source/cif2x/tutorial/index.rst`
- Modify: `docs/ja/source/cif2x/tutorial/index.rst`

- [ ] **Step 1: EN — add per-target examples after the QE sweep example**

In `docs/en/source/cif2x/tutorial/index.rst`, immediately AFTER the existing
sentence that begins "When ``K_POINTS`` is given as above, ... ``12x12x12/``,
respectively." (the end of the QE sweep example), insert:

```rst

The same ``${...}`` syntax applies to every target, not only Quantum ESPRESSO;
the output sub-directory name is derived from the swept value. As in the example
above, each ``content:`` block below belongs to a ``tasks:`` entry. The swept key
sits wherever the parameter lives in ``content`` — for VASP under ``incar`` (or
``kpoints``/``poscar``/``potcar``); for OpenMX and AkaiKKR the ``content`` keys are
flat.

VASP — cutoff convergence:

.. code-block:: yaml

   content:
     incar:
       ENCUT: ${ [400, 600, 800] }

generates the input files in ``400/``, ``600/``, ``800/``.

OpenMX — energy-cutoff convergence:

.. code-block:: yaml

   content:
     scf.energycutoff: ${ [150, 200, 250] }

generates the input files in ``150/``, ``200/``, ``250/``.

AkaiKKR — Brillouin-zone-quality convergence:

.. code-block:: yaml

   content:
     bzqlty: ${ [12, 16, 20] }

generates the input files in ``12/``, ``16/``, ``20/``.
```

- [ ] **Step 2: JA — add the same after the JA QE sweep example**

In `docs/ja/source/cif2x/tutorial/index.rst`, immediately AFTER the explanatory
sentence that follows the QE `grid: ${ [ [4,4,4], ... ] }` example in the
"パラメータセットを指定する" section, insert:

```rst

同じ ``${...}`` 構文は Quantum ESPRESSO だけでなく全ターゲットで利用でき、出力
サブディレクトリ名はスイープ値から生成されます。以下の各 ``content:`` ブロックは、
前述の例と同様に ``tasks:`` の各エントリ配下に置かれます。スイープ対象のキーは
``content`` 内でそのパラメータが定義される階層に記述します。VASP では ``incar``
(または ``kpoints``/``poscar``/``potcar``)の下、OpenMX と AkaiKKR では ``content``
直下のフラットなキーです。

VASP(カットオフ収束):

.. code-block:: yaml

   content:
     incar:
       ENCUT: ${ [400, 600, 800] }

により ``400/``, ``600/``, ``800/`` に入力ファイルが生成されます。

OpenMX(エネルギーカットオフ収束):

.. code-block:: yaml

   content:
     scf.energycutoff: ${ [150, 200, 250] }

により ``150/``, ``200/``, ``250/`` に入力ファイルが生成されます。

AkaiKKR(ブリルアンゾーン品質の収束):

.. code-block:: yaml

   content:
     bzqlty: ${ [12, 16, 20] }

により ``12/``, ``16/``, ``20/`` に入力ファイルが生成されます。
```

- [ ] **Step 3: Lint (best-effort) and commit**

Run (only if docutils is importable):
```bash
python3 -c "import docutils" 2>/dev/null && for f in docs/en/source/cif2x/tutorial/index.rst docs/ja/source/cif2x/tutorial/index.rst; do python3 -m docutils "$f" /dev/null 2>&1 | grep -iE "severe|error" || echo "$f ok"; done || echo "docutils unavailable"
```
Expected: `... ok` for both (no severe/error lines).

```bash
git add docs/en/source/cif2x/tutorial/index.rst docs/ja/source/cif2x/tutorial/index.rst
git commit -m "docs: add VASP/OpenMX/AkaiKKR parameter-sweep examples to the tutorial"
```
End the commit body with a blank line then:
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

---

## Task 3: Full-suite verification

**Files:** none (verification only)

- [ ] **Step 1: Run the entire test suite**

Run: `pytest -q`
Expected: all tests pass, including the new and strengthened `test_inflate.py` cases.

- [ ] **Step 2: Final commit (only if fixups were needed)**

```bash
git add -A && git commit -m "test: verify parameter-sweep coverage" || echo "nothing to commit"
```

---

## Self-review notes

- **Spec coverage:** per-target tutorial examples EN+JA with the tasks-context + YAML-structure note and unquoted syntax (Task 2); required regression cases for VASP/OpenMX/AkaiKKR via `cif2x.utils.inflate` with module-level `pytest.importorskip` and paired `(dirkey, value)` assertions, plus paired assertions added to the existing QE cases (Task 1); full-suite verification (Task 3). All spec sections map to a task.
- **Verified values:** every asserted tuple was confirmed by running `inflate` against the real content classes before writing the plan.
- **No production-code change:** the sweep already works for all targets; QE keeps its separate `inflate` (documented divergence risk, out of scope).
- **Name/accessor consistency:** `utils_inflate`/`qe_inflate`, `Content(...)`/`Content.from_dict(...)`, and accessors `c["incar"]["ENCUT"]`/`c["bzqlty"]`/`c["scf.energycutoff"]` match the verified behavior across tasks.
