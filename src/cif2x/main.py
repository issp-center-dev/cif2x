#!/usr/bin/env python3

import sys
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
from cif2x.input_validator import (
    InputValidationError,
    normalize_target,
    validate_input,
)


def main():
    import argparse

    parser = argparse.ArgumentParser(prog="cif2x")
    parser.add_argument("input_file", action="store", help="input parameter file (input.yaml)")
    parser.add_argument("cif_file", action="store", help="CIF file (data.cif)")
    parser.add_argument("--version", action="version", version="%(prog)s version {}".format(__version__))
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-t", "--target", action="store", required=True, help="target application. Supported targets: quantum_espresso (qe, espresso), vasp, openmx, akaikkr. (case-insensitive)")
    parser.add_argument("--dry-run", action="store_true", default=False, help="print the generated input files to stdout instead of writing them")

    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING-(args.verbose-args.quiet)*10)

    try:
        _run(args)
    except InputValidationError as e:
        logger.error(str(e))
        sys.exit(1)


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


def _run(args):
    target = normalize_target(args.target)

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
    struct = Cif2Struct(args.cif_file, info_dict.get("structure") or {})

    # `optional:` written with no entries parses to None; treat it as empty so
    # generators that read optional.get(...) do not crash on None.
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
