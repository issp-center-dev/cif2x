# cif2x RESPACK target (`-t respack`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `respack` target that emits the RESPACK workflow in one run — QE `scf`/`nscf`/`bands` (reusing `Struct2QE`) plus the RESPACK control `input.in` (new `Struct2RESPACK`).

**Architecture:** `-t respack` routes each `tasks:` entry by `mode` (QE modes → `Struct2QE`; `mode: respack` → `Struct2RESPACK`). `Struct2RESPACK` merges a namelist-only template with `content`, validates the physics it must (`N_wannier`, windows, `dense`, SCDM), auto-generates only the `&param_interpolation` k-path from `HighSymmKpath`, and renders `input.in` (namelists in RESPACK order with the k-coords block spliced after `&param_interpolation`). No production change to `Struct2QE`.

**Tech Stack:** Python, f90nml, pymatgen (`HighSymmKpath`), pytest. Run from repo root (`tests/conftest.py` adds `src/`).

---

## Background facts (verified)

- `input_validator.TARGETS` is a per-target table (aliases/required/allowed); `normalize_target` canonicalizes `-t`. `main._run` builds `struct = Cif2Struct(...)`, then `generator_cls = _generator_class(target)` and loops `tasks`. `_generator_class` resolves classes from module globals at call time.
- `Cif2Struct` stores `use_ibrav` (`cif2struct.py:33`) but NOT `use_primitive` (used inline at `:93`); the `structure:` block is `info_dict["structure"]`, available in `_run` as `struct_params`.
- `cards._generate_band_path` reads `HighSymmKpath(structure).kpath` → `{"kpoints":{label:coords}, "path":[[label,...],...]}` and surfaces pymatgen warnings through the logger (skipping `DeprecationWarning`).
- `f90nml.Namelist({name: dict}).write(stream)` renders one namelist; keys are lowercased and alphabetically ordered (RESPACK reads case-insensitively and order-independently). A list value renders as `key = a, b`.
- RESPACK `input.in` `&param_interpolation`: `N_sym_points` is a list (one positive int ≥2 per disjoint path); the k-coords block (format 0) follows the `/`, one `kx ky kz  ! label` line per high-symmetry point, segment by segment (verified in `respack-wannier-py` `input_parser.py:108-230`).
- `utils.dryrun_emit(path, content)` prints a `# === path ===` header + content (used by other generators' `--dry-run`).

## File structure

- **Create** `src/cif2x/struct2respack.py` — `Struct2RESPACK` + helpers.
- **Modify** `src/cif2x/input_validator.py` — add `respack` to `TARGETS`; respack mode-sensitive validation hook.
- **Modify** `src/cif2x/main.py` — import `Struct2RESPACK`; respack per-task routing + `use_primitive`/`use_ibrav` guard.
- **Create** `tests/test_respack.py` — generator + validation tests.
- **Modify** `tests/test_main_cli.py` — respack dispatch tests.
- **Modify** docs (`command/index.rst` en/ja) + add `sample/cif2x/respack/<case>/`.

---

## Task 1: `Struct2RESPACK` — generate `input.in`

**Files:**
- Create: `src/cif2x/struct2respack.py`
- Test: `tests/test_respack.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_respack.py`:

```python
import io

import pytest

pytest.importorskip("pymatgen")
pytest.importorskip("f90nml")

import f90nml
from pymatgen.core import Lattice, Structure

from cif2x.input_validator import InputValidationError
from cif2x.struct2respack import Struct2RESPACK


class _Struct:
    def __init__(self, structure):
        self.structure = structure


def _cubic():
    return Structure(Lattice.cubic(3.8), ["Cu"], [[0, 0, 0]])


_TEMPLATE = """&param_chiqw
Ecut_for_eps = 5.0
flg_cRPA = 1
/
&param_wannier
N_wannier = 3
Lower_energy_window = 11.0
Upper_energy_window = 14.2
/
&param_interpolation
dense = 8, 8, 8
/
&param_calc_int
/
"""


def _write(tmp_path, text=_TEMPLATE):
    p = tmp_path / "respack.in_tmpl"
    p.write_text(text)
    return str(p)


def _render(tmp_path, structure=None, content=None, template=None):
    params = {"template": template or _write(tmp_path), "content": content or {}}
    gen = Struct2RESPACK(params, _Struct(structure or _cubic()))
    return gen.render()


def test_render_passes_through_n_wannier_and_windows(tmp_path):
    out = _render(tmp_path)
    nml = f90nml.reads(out)
    assert nml["param_wannier"]["n_wannier"] == 3
    assert nml["param_wannier"]["lower_energy_window"] == 11.0
    assert nml["param_wannier"]["upper_energy_window"] == 14.2


def test_render_forces_scdm(tmp_path):
    out = _render(tmp_path)
    nml = f90nml.reads(out)
    assert nml["param_wannier"]["n_initial_guess"] == 0
    # no Gaussian initial-guess line (a bare 'dxy ...' token) is emitted
    assert "dxy" not in out


def test_render_kpath_matches_highsymmkpath_segments(tmp_path):
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    s = _cubic()
    out = _render(tmp_path, structure=s)
    nml = f90nml.reads(out)
    kpath = HighSymmKpath(s).kpath
    seg_lengths = [len(seg) for seg in kpath["path"]]
    nsp = nml["param_interpolation"]["n_sym_points"]
    nsp = [nsp] if isinstance(nsp, int) else list(nsp)
    assert nsp == seg_lengths
    # the k-coords block after &param_interpolation has sum(seg_lengths) lines
    after = out.split("&param_interpolation", 1)[1]
    coord_lines = [ln for ln in after.splitlines()
                   if ln.strip() and not ln.strip().startswith(("/", "&"))
                   and "=" not in ln]
    assert len(coord_lines) == sum(seg_lengths)


def test_missing_n_wannier_rejected(tmp_path):
    t = _TEMPLATE.replace("N_wannier = 3\n", "")
    with pytest.raises(InputValidationError, match="N_wannier"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_nonzero_n_initial_guess_rejected(tmp_path):
    t = _TEMPLATE.replace("N_wannier = 3", "N_wannier = 3\nN_initial_guess = 3")
    with pytest.raises(InputValidationError, match="N_initial_guess"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_window_order_rejected(tmp_path):
    t = _TEMPLATE.replace("Upper_energy_window = 14.2", "Upper_energy_window = 9.0")
    with pytest.raises(InputValidationError, match="window"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_non_param_namelist_rejected(tmp_path):
    t = _TEMPLATE + "&system\nnoncolin = .true.\n/\n"
    with pytest.raises(InputValidationError, match="namelist"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_trailing_content_rejected(tmp_path):
    t = _TEMPLATE + "garbage line not in a namelist\n"
    with pytest.raises(InputValidationError, match="namelist-only|outside"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_bad_dense_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = 8, 8")
    with pytest.raises(InputValidationError, match="dense"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_reading_sk_format_nonzero_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = 8, 8, 8\nreading_sk_format = 1")
    with pytest.raises(InputValidationError, match="reading_sk_format"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_write_input_dry_run(tmp_path, capsys):
    params = {"template": _write(tmp_path), "content": {}}
    gen = Struct2RESPACK(params, _Struct(_cubic()))
    gen.write_input("input.in", str(tmp_path), dry_run=True)
    out = capsys.readouterr().out
    assert "input.in" in out and "&param_wannier" in out
    assert not (tmp_path / "input.in").exists()
    gen.write_input("input.in", str(tmp_path), dry_run=False)
    assert (tmp_path / "input.in").exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_respack.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'cif2x.struct2respack'`.

- [ ] **Step 3: Implement `src/cif2x/struct2respack.py`**

```python
"""Generate a RESPACK ``input.in`` control file from a namelist-only template.

``N_wannier`` and the physics (energy windows, dense, cRPA params) come from the
template/``content``; cif2x only auto-generates the ``&param_interpolation``
high-symmetry k-path (from pymatgen) and forces SCDM (``N_initial_guess = 0``).
"""

import io
import logging
import warnings
from pathlib import Path

import f90nml

from cif2x.input_validator import InputValidationError
from cif2x.utils import dryrun_emit

logger = logging.getLogger(__name__)

_ALLOWED_NAMELISTS = (
    "param_chiqw", "param_wannier", "param_interpolation",
    "param_visualization", "param_calc_int",
)
_NAMELIST_ORDER = _ALLOWED_NAMELISTS


def _deepupdate(dst, src):
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deepupdate(dst[k], v)
        else:
            dst[k] = v


def _has_content_outside_namelists(text):
    """True if a non-blank, non-comment line sits outside a &name.../ block."""
    inside = False
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if inside:
            if line == "/" or line.endswith("/"):
                inside = False
        else:
            if line.startswith("&"):
                inside = not (line.endswith("/") or "/" in line[1:])
            else:
                return True
    return False


def _kpath(structure):
    """Return (n_sym_points list, [(label, (kx,ky,kz)), ...]) from HighSymmKpath."""
    from pymatgen.symmetry.bandstructure import HighSymmKpath

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        kpath = HighSymmKpath(structure).kpath
    for w in caught:
        if not issubclass(w.category, DeprecationWarning):
            logger.warning("band path: %s", w.message)
    coords = kpath["kpoints"]
    segments = kpath["path"]
    n_sym_points = [len(seg) for seg in segments]
    rows = []
    for seg in segments:
        for label in seg:
            c = coords[label]
            rows.append((label, (float(c[0]), float(c[1]), float(c[2]))))
    return n_sym_points, rows


def _as_three_ints(value, *, key):
    if isinstance(value, int):
        return [value, value, value]
    if isinstance(value, (list, tuple)) and len(value) == 3 and all(
        isinstance(v, int) for v in value
    ):
        return list(value)
    raise InputValidationError(
        f"respack: '{key}' must be three integers (got {value!r}).")


class Struct2RESPACK:
    def __init__(self, params, struct):
        self.params = params
        self.struct = struct
        self._data, self._n_sym_points, self._kpath_rows = self._build()

    def _build(self):
        template = self.params.get("template")
        if not template:
            raise InputValidationError("respack task: 'template' is required.")
        path = Path(template)
        if not path.exists():
            raise InputValidationError(f"respack template not found: {template}")
        text = path.read_text()
        try:
            nml = f90nml.reads(text)
        except Exception as e:
            raise InputValidationError(
                f"failed to parse RESPACK template '{template}': {e}")
        for name in nml:
            if name not in _ALLOWED_NAMELISTS:
                raise InputValidationError(
                    f"respack template: unexpected namelist '&{name}'. Only RESPACK "
                    f"&param_* namelists are allowed (namelist-only templates).")
        if _has_content_outside_namelists(text):
            raise InputValidationError(
                "respack template: content outside &param_* namelists is not "
                "allowed (templates are namelist-only; cif2x generates the "
                "k-coords block).")

        data = {k: dict(v) for k, v in nml.items()}
        content = self.params.get("content") or {}
        norm = {str(k).lower(): {str(kk).lower(): vv for kk, vv in (v or {}).items()}
                for k, v in content.items()}
        _deepupdate(data, norm)

        wann = data.setdefault("param_wannier", {})
        n_wannier = wann.get("n_wannier")
        if not isinstance(n_wannier, int) or n_wannier <= 0:
            raise InputValidationError(
                "respack: &param_wannier 'N_wannier' must be a positive integer "
                f"(got {n_wannier!r}); it is user physics and is not auto-derived.")
        if wann.get("n_initial_guess", 0) not in (0, None):
            raise InputValidationError(
                "respack: v1 is SCDM-only — 'N_initial_guess' must be 0 (got "
                f"{wann.get('n_initial_guess')!r}).")
        wann["n_initial_guess"] = 0
        lo = wann.get("lower_energy_window")
        hi = wann.get("upper_energy_window")
        if not isinstance(lo, (int, float)) or not isinstance(hi, (int, float)):
            raise InputValidationError(
                "respack: 'Lower_energy_window' and 'Upper_energy_window' are "
                "required numeric values in &param_wannier.")
        if lo >= hi:
            raise InputValidationError(
                f"respack: energy window must have lower < upper (got {lo} >= {hi}).")

        interp = data.setdefault("param_interpolation", {})
        if int(interp.get("reading_sk_format", 0)) != 0:
            raise InputValidationError(
                "respack: only reading_sk_format = 0 is supported.")
        interp["dense"] = _as_three_ints(interp.get("dense"), key="dense")

        n_sym_points, rows = _kpath(self.struct.structure)
        interp["n_sym_points"] = n_sym_points
        return data, n_sym_points, rows

    def render(self):
        out = []
        for name in _NAMELIST_ORDER:
            if name not in self._data:
                continue
            buf = io.StringIO()
            f90nml.Namelist({name: self._data[name]}).write(buf)
            out.append(buf.getvalue())
            if name == "param_interpolation":
                for label, (x, y, z) in self._kpath_rows:
                    out.append("{:.6f} {:.6f} {:.6f}  ! {}\n".format(x, y, z, label))
        return "".join(out)

    def write_input(self, output_file, output_dir, dry_run=False):
        content = self.render()
        if dry_run:
            dryrun_emit(str(Path(output_dir, output_file)), content)
            return
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        with open(Path(output_dir, output_file), "w") as fp:
            fp.write(content)
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest tests/test_respack.py -q`
Expected: PASS (all generator + validation tests).

- [ ] **Step 5: Commit**

```bash
git add src/cif2x/struct2respack.py tests/test_respack.py
git commit -m "feat(respack): add Struct2RESPACK input.in generator"
```
End the body with a blank line then:
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

---

## Task 2: Dispatch + validation (`-t respack` per-task routing)

**Files:**
- Modify: `src/cif2x/input_validator.py`
- Modify: `src/cif2x/main.py`
- Test: `tests/test_main_cli.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_main_cli.py`:

```python
def test_respack_routes_qe_and_respack_tasks(monkeypatch, tmp_path):
    pytest.importorskip("pymatgen")
    captured = {"qe": [], "respack": []}

    class _QE:
        def __init__(self, params, struct):
            pass

        def write_input(self, *a, **k):
            captured["qe"].append(True)

    class _RP:
        def __init__(self, params, struct):
            pass

        def write_input(self, *a, **k):
            captured["respack"].append(True)

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2QE", _QE)
    monkeypatch.setattr(main_mod, "Struct2RESPACK", _RP)
    yf = tmp_path / "input.yaml"
    yf.write_text(
        "structure:\n  use_primitive: true\n"
        "tasks:\n"
        "  - mode: scf\n    output_file: scf.in\n"
        "  - mode: respack\n    output_file: input.in\n    template: r.in_tmpl\n"
    )
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "respack", str(yf), str(cif)])
    main_mod.main()
    assert captured["qe"] == [True]
    assert captured["respack"] == [True]


def test_respack_rejects_non_workflow_mode(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    yf = tmp_path / "input.yaml"
    yf.write_text(
        "structure:\n  use_primitive: true\n"
        "tasks:\n  - mode: relax\n    output_file: x.in\n"
    )
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "respack", str(yf), str(cif)])
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            main_mod.main()
    assert "relax" in caplog.text


def test_respack_requires_primitive_cell(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    yf = tmp_path / "input.yaml"
    yf.write_text(
        "tasks:\n  - mode: respack\n    output_file: input.in\n    template: r\n"
    )
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "respack", str(yf), str(cif)])
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            main_mod.main()
    assert "use_primitive" in caplog.text
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_main_cli.py -q`
Expected: FAIL — `respack` is not a known target (normalize_target raises), `Struct2RESPACK` not imported, no respack routing.

- [ ] **Step 3: Add `respack` to the validator**

In `src/cif2x/input_validator.py`, add to `TARGETS` (after the `akaikkr` entry):

```python
    "respack": {
        "aliases": ("respack",),
        "required": ("mode", "output_file"),
        "allowed": _COMMON_ALLOWED | {"mode"},
    },
```

(`_COMMON_ALLOWED` already includes `template`, `content`, `output_file`, `output_dir`, `optional`, which covers both QE-mode tasks and the respack task. Per-mode legality — only `scf`/`nscf`/`bands`/`respack` — is enforced at dispatch; the respack task's required `template` is enforced by `Struct2RESPACK.__init__` (a clean `InputValidationError`), a deliberate v1 simplification rather than mode-specific validator schemas.)

- [ ] **Step 4: Wire respack routing into `main.py`**

In `src/cif2x/main.py`, add the import after `from cif2x.struct2akaikkr import Struct2AkaiKKR`:

```python
from cif2x.struct2respack import Struct2RESPACK
```

Then replace the dispatch loop. Find, in `_run`, the block:

```python
    info_optional = info_dict.get("optional") or {}
    generator_cls = _generator_class(target)

    for idx, info in enumerate(info_dict["tasks"], start=1):
        logger.info(f"start task {idx}")

        params = {}
        params.update(info)
        if params.get("optional") is None:
            params["optional"] = {}
        deepupdate(params, {"optional": info_optional})

        output_file = info.get("output_file")
        output_dir = info.get("output_dir", ".")

        generator = generator_cls(params, struct)
        generator.write_input(output_file, output_dir, dry_run=args.dry_run)
```

and replace it with:

```python
    info_optional = info_dict.get("optional") or {}
    generator_cls = None if target == "respack" else _generator_class(target)

    if target == "respack":
        if not struct_params.get("use_primitive") or struct_params.get("use_ibrav"):
            raise InputValidationError(
                "target 'respack' requires structure.use_primitive: true and "
                "use_ibrav: false (the high-symmetry k-path needs the primitive cell).")

    for idx, info in enumerate(info_dict["tasks"], start=1):
        logger.info(f"start task {idx}")

        params = {}
        params.update(info)
        if params.get("optional") is None:
            params["optional"] = {}
        deepupdate(params, {"optional": info_optional})

        output_file = info.get("output_file")
        output_dir = info.get("output_dir", ".")

        if target == "respack":
            gen_cls = _respack_generator(info.get("mode"), idx)
        else:
            gen_cls = generator_cls

        generator = gen_cls(params, struct)
        generator.write_input(output_file, output_dir, dry_run=args.dry_run)
```

Then add this helper next to `_generator_class`:

```python
def _respack_generator(mode, taskid):
    """Pick the generator for a task under -t respack, by mode."""
    if mode in ("scf", "nscf", "bands"):
        return Struct2QE
    if mode == "respack":
        return Struct2RESPACK
    raise InputValidationError(
        f"task {taskid}: unsupported mode '{mode}' for target 'respack' "
        "(use scf, nscf, bands, or respack).")
```

(The `struct_params` name already exists in `_run` — it is `info_dict.get("structure") or {}`, computed just before `struct` is built.)

- [ ] **Step 5: Run to verify it passes**

Run: `pytest tests/test_main_cli.py -q`
Expected: PASS (routing, mode-reject, and primitive-cell tests, plus the existing CLI tests).

- [ ] **Step 6: Commit**

```bash
git add src/cif2x/input_validator.py src/cif2x/main.py tests/test_main_cli.py
git commit -m "feat(respack): add -t respack target with per-task mode routing"
```
End the body with a blank line then:
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

---

## Task 3: Documentation + sample

**Files:**
- Modify: `docs/en/source/cif2x/command/index.rst`, `docs/ja/source/cif2x/command/index.rst`
- Create: `sample/cif2x/respack/SrVO3/{input.yaml, respack.in_tmpl, README.md}`

- [ ] **Step 1: Document the target (EN)**

In `docs/en/source/cif2x/command/index.rst`, in the `-t` *target* list, add a bullet after the AkaiKKR one:

```rst
    - ``respack``: generates the RESPACK workflow — the Quantum ESPRESSO inputs (``scf``/``nscf``/``bands`` tasks) and the RESPACK control file ``input.in`` (a ``mode: respack`` task). The structure must be the standardized primitive cell (``structure.use_primitive: true``, ``use_ibrav: false``). ``N_wannier`` and the energy windows are user physics supplied in the ``content``/template; initial guesses use SCDM (``N_initial_guess = 0``), targeting the ``respack-wannier-py`` port. The ``nscf`` task must set ``nosym = .true.``/``noinv = .true.`` for ``qe2respack``. Run order: ``scf`` → ``nscf`` → ``qe2respack`` → ``respack``.
```

- [ ] **Step 2: Document the target (JA)**

In `docs/ja/source/cif2x/command/index.rst`, add the equivalent bullet after the AkaiKKR one:

```rst
    - ``respack``: RESPACK ワークフロー一式を生成します。Quantum ESPRESSO 入力(``scf``/``nscf``/``bands`` タスク)と RESPACK 制御ファイル ``input.in``(``mode: respack`` タスク)を出力します。構造は標準プリミティブセル(``structure.use_primitive: true``、``use_ibrav: false``)である必要があります。``N_wannier`` とエネルギー窓は ``content``/テンプレートで与える物理量で、初期ゲスは SCDM(``N_initial_guess = 0``、``respack-wannier-py`` を対象)です。``nscf`` タスクには ``qe2respack`` 用に ``nosym = .true.``/``noinv = .true.`` を設定してください。実行順序: ``scf`` → ``nscf`` → ``qe2respack`` → ``respack``。
```

- [ ] **Step 3: Create the sample**

Create `sample/cif2x/respack/SrVO3/respack.in_tmpl`:

```
&param_chiqw
Ecut_for_eps = 5.0
flg_cRPA = 1
/
&param_wannier
N_wannier = 3
Lower_energy_window = 11.0
Upper_energy_window = 14.2
/
&param_interpolation
dense = 8, 8, 8
/
&param_calc_int
/
```

Create `sample/cif2x/respack/SrVO3/input.yaml`:

```yaml
structure:
  use_primitive: true
  use_ibrav: false

optional:
  pseudo_dir: ./pseudo
  pp_file: pp.csv

tasks:
  - mode: scf
    output_file: scf.in
    content:
      control: { calculation: scf }
      K_POINTS: { option: automatic, grid: [6, 6, 6] }
  - mode: nscf
    output_file: nscf.in
    content:
      control: { calculation: nscf }
      system: { nosym: true, noinv: true }
      K_POINTS: { option: crystal, grid: [6, 6, 6] }
  - mode: respack
    output_file: input.in
    template: respack.in_tmpl
```

Create `sample/cif2x/respack/SrVO3/README.md`:

```markdown
# SrVO3 RESPACK sample

`cif2x -t respack input.yaml SrVO3.cif` emits `scf.in`, `nscf.in`, and `input.in`.
Provide `SrVO3.cif`, `pp.csv`, and the pseudopotentials under `./pseudo`.
`N_wannier = 3` (t2g) and the energy windows are physics choices in
`respack.in_tmpl`. Run order: scf → nscf → `qe2respack` → `respack` (the
`respack-wannier-py` port; SCDM, `N_initial_guess = 0`).
```

- [ ] **Step 4: Lint (best-effort) and commit**

Run: `python3 -c "import docutils" 2>/dev/null && for f in docs/en/source/cif2x/command/index.rst docs/ja/source/cif2x/command/index.rst; do python3 -m docutils "$f" /dev/null 2>&1 | grep -iE "severe|error" || echo "$f ok"; done || echo "docutils unavailable"`

```bash
git add docs/en/source/cif2x/command/index.rst docs/ja/source/cif2x/command/index.rst sample/cif2x/respack
git commit -m "docs(respack): document the respack target and add a SrVO3 sample"
```
End the body with a blank line then:
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

---

## Task 4: Full-suite + manual verification

**Files:** none (verification only)

- [ ] **Step 1: Run the entire test suite**

Run: `pytest -q`
Expected: all tests pass (existing + `test_respack.py` + new CLI tests); no regression to QE/other targets.

- [ ] **Step 2: Manual end-to-end render (dry-run, real structure)**

Run from the sample dir (no QE/pseudos needed — `--dry-run` for the respack task path):
```bash
python3 -c "
import sys; sys.path.insert(0, 'src')
from pymatgen.core import Structure, Lattice
from cif2x.struct2respack import Struct2RESPACK
class S:
    def __init__(s, st): s.structure = st
st = Structure(Lattice.cubic(3.84), ['Sr','V','O','O','O'],
               [[0.5,0.5,0.5],[0,0,0],[0.5,0,0],[0,0.5,0],[0,0,0.5]])
tmpl='sample/cif2x/respack/SrVO3/respack.in_tmpl'
print(Struct2RESPACK({'template': tmpl, 'content': {}}, S(st)).render())
" | head -30
```
Expected: a valid `input.in` — `&param_wannier` with `n_wannier = 3`, `n_initial_guess = 0`; `&param_interpolation` with an `n_sym_points` list and a k-coords block of `kx ky kz  ! label` lines; no traceback.

- [ ] **Step 3: Final commit (if any fixups were needed)**

```bash
git add -A && git commit -m "test(respack): verify end-to-end" || echo "nothing to commit"
```

---

## Self-review notes

- **Spec coverage:** `Struct2RESPACK` namelist-only template + content merge, N_wannier/window/dense/SCDM/reading_sk_format validation, multi-segment k-path from `HighSymmKpath`, render with k-block spliced after `&param_interpolation`, dry-run (Task 1); `-t respack` per-task mode routing + primitive-cell guard + validator entry (Task 2); docs (en/ja) + sample with nosym/noinv + run order + SCDM/port notes (Task 3); full-suite + manual render (Task 4). All spec sections map to a task.
- **No QE change:** `Struct2QE._find_nbnd_info` is untouched (N_wannier is user physics); QE-mode tasks under respack reuse `Struct2QE` verbatim.
- **Name consistency:** `Struct2RESPACK`, `_respack_generator`, `_kpath`, `_has_content_outside_namelists`, `_as_three_ints`, `TARGETS["respack"]` are consistent across tasks/tests.
- **Monkeypatch safety:** `main.py` resolves `Struct2QE`/`Struct2RESPACK` from module globals at call time (via `_respack_generator`), so the CLI tests' monkeypatching works.
