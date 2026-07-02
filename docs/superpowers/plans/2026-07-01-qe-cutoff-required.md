# QE cutoff required (reject instead of 0.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** When cif2x cannot resolve a QE plane-wave cutoff for an element, raise a clear `InputValidationError` instead of silently generating `ecutwfc = 0.0`.

**Architecture:** `Struct2QE._find_elem_cutoff` is the single per-element resolver that chains the cutoff CSV then the UPF file. Its final `None → 0.0` substitution becomes a raised error naming the element; the dead commented-out `raise`/`logger.error` tail in `_find_elem_cutoff_from_file` is removed since the caller now enforces the requirement.

**Tech Stack:** Python, pytest, pandas, BeautifulSoup (bs4), pymatgen.

---

### Task 1: Reject unresolved cutoffs

**Files:**
- Modify: `src/cif2x/struct2qe.py` (add import; `_find_elem_cutoff` ~lines 87-101; `_find_elem_cutoff_from_file` tail ~lines 154-157)
- Test: `tests/test_qe_pp_filename.py` (append)

Context — the current resolver (`src/cif2x/struct2qe.py`):

```python
    def _find_elem_cutoff(self, ename):
        ecutwfc, ecutrho = None, None

        if ecutwfc is None or ecutrho is None:
            ecutwfc, ecutrho = self._find_elem_cutoff_from_table(ename)

        if ecutwfc is None or ecutrho is None:
            ecutwfc, ecutrho = self._find_elem_cutoff_from_file(ename)

        if ecutwfc is None:
            ecutwfc = 0.0
        if ecutrho is None:
            ecutrho = 0.0

        return ecutwfc, ecutrho
```

and the dead tail of `_find_elem_cutoff_from_file`:

```python
        if ecutwfc is None or ecutrho is None:
            # raise ValueError("cutoff information not found")
            logger.error("cutoff information not found")
        return ecutwfc, ecutrho
```

`_find_elem_cutoff_from_table` returns `(None, None)` when `self.cutoff_list is None`; `_find_elem_cutoff_from_file` returns `(None, None)` when `self.pp_list is None`. Both `_find_elem_cutoff_from_table` and `_pp_filename` read `self.is_soc`. `InputValidationError` (a plain `Exception` subclass) lives in `cif2x.input_validator`, which imports nothing from `cif2x` — no circular-import risk. It is not yet imported in `struct2qe.py`.

- [ ] **Step 1: Write the failing test + regression guards**

Append to `tests/test_qe_pp_filename.py`. Add this import near the top (after the existing `from cif2x.struct2qe import Struct2QE`):

```python
from cif2x.input_validator import InputValidationError
```

Then append these tests:

```python
def test_find_elem_cutoff_missing_raises():
    # Neither the cutoff CSV nor a UPF file is available: the resolver must
    # reject with a clear error naming the element, not fall back to 0.0.
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.cutoff_list = None
    qe.pp_list = None
    with pytest.raises(InputValidationError, match="Sr"):
        qe._find_elem_cutoff("Sr")


def test_find_elem_cutoff_resolves_from_table():
    # Regression guard: a CSV hit resolves without error.
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.pp_list = pd.DataFrame({"pseudopotential": ["pbe-spn-rrkjus_psl.1.0.0"]},
                              index=["Fe"])
    qe.cutoff_list = pd.DataFrame(
        {"ecutwfc": [40.0], "ecutrho": [320.0]},
        index=["Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"],
    )
    assert qe._find_elem_cutoff("Fe") == (40.0, 320.0)


def test_find_elem_cutoff_resolves_from_upf(tmp_path):
    # Regression guard: a UPF-header hit resolves without error.
    pytest.importorskip("bs4")
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.pp_list = pd.DataFrame({"pseudopotential": ["pbe-spn-rrkjus_psl.1.0.0"]},
                              index=["Fe"])
    qe.cutoff_list = None
    qe.pseudo_dir = str(tmp_path)
    (tmp_path / "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF").write_text(
        '<UPF><PP_HEADER wfc_cutoff="40.0" rho_cutoff="320.0"/></UPF>')
    assert qe._find_elem_cutoff("Fe") == (40.0, 320.0)
```

- [ ] **Step 2: Run the tests to verify the new one fails**

Run: `python3 -m pytest tests/test_qe_pp_filename.py -v`
Expected: `test_find_elem_cutoff_missing_raises` FAILS (it returns `(0.0, 0.0)` instead of raising, so `pytest.raises` reports "DID NOT RAISE"). The two `resolves_*` guards PASS (behavior unchanged).

- [ ] **Step 3: Add the import in struct2qe.py**

In `src/cif2x/struct2qe.py`, after the existing `from cif2x.utils import dryrun_emit` line (line 12), add:

```python
from cif2x.input_validator import InputValidationError
```

- [ ] **Step 4: Replace the 0.0 fallback with a raise**

In `src/cif2x/struct2qe.py`, change `_find_elem_cutoff` so its tail reads:

```python
    def _find_elem_cutoff(self, ename):
        ecutwfc, ecutrho = None, None

        if ecutwfc is None or ecutrho is None:
            ecutwfc, ecutrho = self._find_elem_cutoff_from_table(ename)

        if ecutwfc is None or ecutrho is None:
            ecutwfc, ecutrho = self._find_elem_cutoff_from_file(ename)

        if ecutwfc is None or ecutrho is None:
            raise InputValidationError(
                "cutoff information not found for element '{}'. Provide it via "
                "optional.cutoff_file or set content.system.ecutwfc/ecutrho "
                "explicitly.".format(ename))

        return ecutwfc, ecutrho
```

- [ ] **Step 5: Remove the dead tail in `_find_elem_cutoff_from_file`**

In `src/cif2x/struct2qe.py`, change the end of `_find_elem_cutoff_from_file` from:

```python
        if ecutwfc is None or ecutrho is None:
            # raise ValueError("cutoff information not found")
            logger.error("cutoff information not found")
        return ecutwfc, ecutrho
```

to:

```python
        return ecutwfc, ecutrho
```

Leave the earlier `logger.error(f"{pseudo_file}: {e}")` file-open diagnostic untouched — it points at a specific unreadable file and is still useful.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `python3 -m pytest tests/test_qe_pp_filename.py -v`
Expected: all tests PASS, including `test_find_elem_cutoff_missing_raises`.

- [ ] **Step 7: Run the full suite as a regression guard**

Run: `python3 -m pytest -q`
Expected: all pass (no sample or test relied on the `0.0` fallback).

- [ ] **Step 8: Commit**

```bash
git add src/cif2x/struct2qe.py tests/test_qe_pp_filename.py
git commit -m "fix(qe): reject unresolved cutoffs instead of defaulting to 0.0"
```

---

## Self-Review

**Spec coverage:**
- "raise when cutoff cannot be resolved, naming the element" → Task 1 Step 4 (raises `InputValidationError` with `ename` and both remedies). ✅
- "drop the dead commented-out raise/logger.error" → Task 1 Step 5. ✅
- "import InputValidationError" → Task 1 Step 3. ✅
- "explicit content / CSV / UPF cases unchanged" → guarded by `test_find_elem_cutoff_resolves_from_table` and `test_find_elem_cutoff_resolves_from_upf`; explicit-content path is not touched because `_find_cutoff_info` is only called when the content key is empty. ✅
- "full existing suite still green" → Task 1 Step 7. ✅

**Placeholder scan:** No TBD/TODO/"handle edge cases"; every code step shows complete code. ✅

**Type consistency:** `_find_elem_cutoff`/`_find_elem_cutoff_from_table`/`_find_elem_cutoff_from_file`/`_pp_filename` names and the `(ecutwfc, ecutrho)` tuple return match the existing code and the tests. `InputValidationError` import path (`cif2x.input_validator`) matches the class location. ✅
