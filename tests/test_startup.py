"""`cif2x -h` / `--version` (and getcif's) must not import pymatgen (issue #1).

The checks run in a fresh interpreter, because the test process itself has
pymatgen loaded by other tests.
"""

import os
import subprocess
import sys

import cif2x

HEAVY = ("pymatgen", "numpy", "pandas", "scipy", "spglib")
# works both for the src/ tree (conftest) and for an installed package
PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(cif2x.__file__)))


def _fresh(code):
    env = dict(os.environ, PYTHONPATH=PKG_ROOT)
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          env=env, check=True)
    return proc.stdout.rstrip().splitlines()[-1]


def test_importing_the_cli_modules_does_not_load_heavy_packages():
    out = _fresh("import sys, cif2x.main, getcif.main\n"
                 f"print(sorted(m for m in {HEAVY!r} if m in sys.modules))")
    assert out == "[]"


def test_help_and_version_do_not_load_pymatgen():
    for module in ("cif2x.main", "getcif.main"):
        for flag in ("-h", "--version"):
            out = _fresh(
                "import sys, runpy\n"
                f"sys.argv = ['x', '{flag}']\n"
                "try:\n"
                f"    runpy.run_module('{module}', run_name='__main__')\n"
                "except SystemExit:\n"
                "    pass\n"
                "print('LOADED' if 'pymatgen' in sys.modules else 'CLEAN')")
            assert out == "CLEAN", f"{module} {flag} imported pymatgen"


def test_help_still_lists_the_targets():
    out = subprocess.run([sys.executable, "-m", "cif2x.main", "-h"], capture_output=True,
                         text=True, env=dict(os.environ, PYTHONPATH=PKG_ROOT), check=True).stdout
    for target in ("quantum_espresso", "vasp", "openmx", "akaikkr", "respack"):
        assert target in out


def test_lazy_attributes_resolve_and_honour_patches(monkeypatch):
    import cif2x.main as main_mod

    assert main_mod.Struct2QE.__name__ == "Struct2QE"       # imported on access
    assert main_mod._generator_class("vasp").__name__ == "Struct2Vasp"
    monkeypatch.setattr(main_mod, "Struct2QE", "patched")
    assert main_mod._lazy("Struct2QE") == "patched"
    assert main_mod._generator_class("quantum_espresso") == "patched"
