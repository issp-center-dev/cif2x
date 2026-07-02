#!/usr/bin/env python3

import sys
import tempfile
from pathlib import Path
from ruamel.yaml import YAML

import logging
logger = logging.getLogger("cif2x")

from cif2x import __version__
from cif2x.cif2struct import Cif2Struct
from cif2x.struct2qe import Struct2QE
from cif2x.struct2vasp import Struct2Vasp
from cif2x.struct2openmx import Struct2OpenMX
from cif2x.struct2akaikkr import Struct2AkaiKKR
from cif2x.struct2respack import Struct2RESPACK
from cif2x.input_validator import (
    InputValidationError,
    normalize_target,
    validate_input,
    _target_choices,
)
from cif2x.mp_source import fetch_to_cif


def main():
    import argparse

    parser = argparse.ArgumentParser(prog="cif2x")
    parser.add_argument("input_file", action="store", help="input parameter file (input.yaml)")
    parser.add_argument("cif_file", action="store", nargs="?", default=None,
                        help="CIF file (data.cif); optional when --mp-id is given")
    parser.add_argument("--mp-id", action="store", default=None,
                        help="Materials Project material id (e.g. mp-149); fetch the structure instead of reading a CIF file")
    parser.add_argument("--symprec", action="store", type=float, default=0.1,
                        help="symmetry tolerance for the fetched structure's CIF (used only with --mp-id; 0 disables)")
    parser.add_argument("--api-key-file", action="store", default="materials_project.key",
                        help="file holding the Materials Project API key (used only with --mp-id)")
    parser.add_argument("--version", action="version", version="%(prog)s version {}".format(__version__))
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-t", "--target", action="store", required=True,
                        help="target application. Supported targets: {}. (case-insensitive)".format(_target_choices()))
    parser.add_argument("--dry-run", action="store_true", default=False, help="print the generated input files to stdout instead of writing them")

    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING-(args.verbose-args.quiet)*10)

    try:
        _run(args)
    except InputValidationError as e:
        logger.error(str(e))
        sys.exit(1)


def _require_one_source(args):
    """Require exactly one structure source: a positional CIF or --mp-id."""
    if args.mp_id and args.cif_file:
        raise InputValidationError("provide either a CIF file or --mp-id, not both.")
    if not args.mp_id and not args.cif_file:
        raise InputValidationError("provide a CIF file or --mp-id.")


def _generator_class(target):
    """Return the Struct2X class for a canonical target name.

    Resolved from module globals at call time so tests can monkeypatch the
    individual Struct2* classes.
    """
    return {
        "quantum_espresso": Struct2QE,
        "vasp": Struct2Vasp,
        "openmx": Struct2OpenMX,
        "akaikkr": Struct2AkaiKKR,
    }[target]


def _respack_generator(mode, taskid):
    """Pick the generator for a task under -t respack, by mode."""
    if mode in ("scf", "nscf", "bands"):
        return Struct2QE
    if mode == "respack":
        return Struct2RESPACK
    raise InputValidationError(
        f"task {taskid}: unsupported mode '{mode}' for target 'respack' "
        "(use scf, nscf, bands, or respack).")


def _run(args):
    target = normalize_target(args.target)
    _require_one_source(args)

    try:
        yaml = YAML(typ="safe")
        with open(Path(args.input_file), mode="r") as fp:
            info_dict = yaml.load(fp)
    except FileNotFoundError:
        raise InputValidationError(f"input file not found: {args.input_file}")
    except InputValidationError:
        raise
    except Exception as e:
        raise InputValidationError(
            f"failed to parse input file '{args.input_file}': {e}"
        )

    validate_input(info_dict, target)

    # `structure:` with no entries parses to None; pass {} so Cif2Struct's
    # params.get(...) does not crash on None.
    struct_params = info_dict.get("structure") or {}
    if args.mp_id:
        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = Path(tmpdir) / "structure.cif"
            fetch_to_cif(args.mp_id, cif_path,
                         symprec=args.symprec, api_key_file=args.api_key_file)
            struct = Cif2Struct(str(cif_path), struct_params)
    else:
        struct = Cif2Struct(args.cif_file, struct_params)

    # `optional:` written with no entries parses to None; treat it as empty so
    # generators that read optional.get(...) do not crash on None.
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

        if target == "respack" and info.get("mode") == "nscf":
            _nml = (info.get("content") or {}).get("namelist") or {}
            _sys = _nml.get("system") or {}
            if not _sys.get("nosym") or not _sys.get("noinv"):
                logger.warning(
                    "respack: the nscf task should set "
                    "content.namelist.system.nosym: true and noinv: true "
                    "(required by qe2respack).")
        gen_cls = _respack_generator(info.get("mode"), idx) if target == "respack" else generator_cls
        generator = gen_cls(params, struct)
        generator.write_input(output_file, output_dir, dry_run=args.dry_run)


def deepupdate(dict1, dict2):
    """
    merge dict2 into dict1; update nested dictionary recursively
    """
    for k, v in dict2.items():
        if isinstance(v, dict) and k in dict1:
            deepupdate(dict1[k], v)
        else:
            dict1[k] = v


if __name__ == '__main__':
    main()
