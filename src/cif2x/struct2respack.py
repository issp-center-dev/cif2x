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
    # NOTE: assumes namelist value lines do not end in a bare '/' (true for
    # RESPACK's numeric/quoted namelists); an unquoted token ending in '/'
    # would falsely close the block. Fails safe (rejects) rather than mis-render.
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
    if isinstance(value, bool):
        raise InputValidationError(
            "respack: '{}' must be three integers (got {!r}).".format(key, value))
    if isinstance(value, int):
        return [value, value, value]
    if (isinstance(value, (list, tuple)) and len(value) == 3
            and all(isinstance(v, int) and not isinstance(v, bool) for v in value)):
        return list(value)
    raise InputValidationError(
        "respack: '{}' must be three integers (got {!r}).".format(key, value))


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
        norm = {}
        for k, v in content.items():
            name = str(k).lower()
            if name not in _ALLOWED_NAMELISTS:
                raise InputValidationError(
                    "respack content: unexpected namelist '{}'. Only RESPACK "
                    "&param_* namelists are allowed.".format(name))
            if v is None:
                norm[name] = {}
            elif isinstance(v, dict):
                norm[name] = {str(kk).lower(): vv for kk, vv in v.items()}
            else:
                raise InputValidationError(
                    "respack content: '{}' must be a mapping of namelist keys.".format(k))
        _deepupdate(data, norm)

        wann = data.setdefault("param_wannier", {})
        n_wannier = wann.get("n_wannier")
        if not isinstance(n_wannier, int) or isinstance(n_wannier, bool) or n_wannier <= 0:
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
        if (not isinstance(lo, (int, float)) or isinstance(lo, bool)
                or not isinstance(hi, (int, float)) or isinstance(hi, bool)):
            raise InputValidationError(
                "respack: 'Lower_energy_window' and 'Upper_energy_window' are "
                "required numeric values in &param_wannier.")
        if lo >= hi:
            raise InputValidationError(
                f"respack: energy window must have lower < upper (got {lo} >= {hi}).")

        interp = data.setdefault("param_interpolation", {})
        rsf = interp.get("reading_sk_format", 0)
        if isinstance(rsf, bool) or (isinstance(rsf, float) and not rsf.is_integer()):
            raise InputValidationError(
                "respack: 'reading_sk_format' must be an integer (got {!r}).".format(rsf))
        try:
            rsf = int(rsf)
        except (TypeError, ValueError):
            raise InputValidationError(
                "respack: 'reading_sk_format' must be an integer (got {!r}).".format(
                    interp.get("reading_sk_format")))
        if rsf != 0:
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
