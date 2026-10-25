#!/usr/bin/env python3

"""
NAME
    cif2x-pp-cutoff - extract cutoff information from pseudo-potential files

SYNOPSIS
    cif2x-pp-cutoff [-f list_file] [-o output_csv] [pp_file ...]

DESCRIPTION
    read pseudo-potential files in UPF format (e.g. for Quantum ESPRESSO)
    and extract cutoff information denoted by wfc_cutoff and rho_cutoff.
    the results are stored in CSV format.

    -f, --filelist list_file
        specify filename that contains list of pseudo-potential files.
        if "-" is specified, the list is read from standard input.

    -o, --output output_csv
        specify filename to store the results in csv format. if omitted, 
        the results are written to standard output.

    --use-suggested
        use "suggested" values in the text header given in the form, e.g.
           Suggested minimum cutoff for wavefunctions:  64. Ry
           Suggested minimum cutoff for charge density: 782. Ry
        rather than the values in pp_header attributes.

    pp_file ...
        filename(s) of pseudo-potential files.

    -h, --help
        display help information and exit.

"""

from bs4 import BeautifulSoup
import csv
import logging
logger = logging.getLogger(__name__)

from cif2x.utils import normalize_cutoff as _normalize_cutoff  # noqa: E402


def read_pseudo_cutoff(pp_file, use_header=True):
    ecutwfc = None
    ecutrho = None

    try:
        with open(pp_file, "r") as fp:
            bs = BeautifulSoup(fp, "html.parser")
    except Exception as e:
        logger.error("{}: {}".format(pp_file, e))
        return None

    header = bs.find("pp_header") if bs is not None else None
    info = bs.find("pp_info") if bs is not None else None
    if header is None:
        logger.error("{}: unknown pseudopotential file format".format(pp_file))
        return None

    if use_header:
        # read cutoff from pp_header attribute
        if header.has_attr("wfc_cutoff"):
            ecutwfc = _normalize_cutoff(header["wfc_cutoff"])
        if header.has_attr("rho_cutoff"):
            ecutrho = _normalize_cutoff(header["rho_cutoff"])
        if ecutwfc is None or ecutrho is None:
            logger.warning("{}: cutoff info not found in header".format(pp_file))

    if ecutwfc is None or ecutrho is None:
        # read cutoff from pp_info text description; fill only the missing
        # fields so header values survive an unset "Suggested ... 0.0" line
        if info is not None:
            lines = str(info).splitlines()
            sg = [s for s in lines if "Suggested" in s]
            for s in sg:
                if "wavefunctions" in s and ecutwfc is None:
                    ecutwfc = _normalize_cutoff(s.split()[5])
                elif "charge density" in s and ecutrho is None:
                    ecutrho = _normalize_cutoff(s.split()[6])

    if ecutwfc is None:
        logger.warning("{}: wfc_cutoff not found".format(pp_file))
        return None
    if ecutrho is None:
        # leave the CSV cell blank: cif2x reads it as "not specified" and
        # omits ecutrho from the input (pw.x then defaults to 4*ecutwfc)
        logger.warning("{}: rho_cutoff not found; leaving the CSV cell blank".format(pp_file))
        return [ pp_file, ecutwfc, "" ]
    return [ pp_file, ecutwfc, ecutrho ]

def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(prog="cif2x-pp-cutoff")
    parser.add_argument("-f", "--filelist", default=None, metavar="list_file", help="list file that contains list of pseudo-potential files")
    parser.add_argument("-o", "--output", default=None, metavar="output_csv", help="output file to store the result in CSV foramt")
    parser.add_argument("--use-suggested", action="store_true", help="use suggested cutoff value in text info")
    parser.add_argument("pp_file", nargs="*", default=[], help="pseudo-potential file(s)")

    args = parser.parse_args()

    pp_files = []
    if args.filelist is not None:
        if args.filelist == '-':
            pp_files += [s.strip() for s in sys.stdin.readlines()]
        else:
            with open(args.filelist, "r") as f:
                pp_files += [s.strip() for s in f.readlines()]
    if len(args.pp_file) > 0:
        pp_files += args.pp_file

    tbl = [[ "pseudofile", "ecutwfc", "ecutrho" ]]
    for f in pp_files:
        x = read_pseudo_cutoff(f, (not args.use_suggested))
        if x is not None:
            tbl.append(x)

    if args.output is None or args.output == '-':
        writer = csv.writer(sys.stdout)
        writer.writerows(tbl)
    else:
        with open(args.output, "w") as f:
            writer = csv.writer(f)
            writer.writerows(tbl)

if __name__ == "__main__":
    main()
