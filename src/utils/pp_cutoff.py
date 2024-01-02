#!/usr/bin/env python3

"""
NAME
    pp_cutoff.py - extract cutoff information from pseudo-potential files

SYNOPSIS
    python3 pp_cutoff.py [-f list_file] [-o output_csv] [pp_file ...]

DESCRIPTION
    read pseudo-potential files in UPF format (e.g. for Quantum ESPRESSO)
    and extract cutoff information denoted by wfc_cutoff and rho_cutoff.
    the results are stored in CSV format.

    -f, --filelist list_file
        specify filename that contains list of pseudo-potential files.

    -o, --output output_csv
        specify filename to store the results in csv format. if omitted, 
        the results are written to standard output.

    pp_file ...
        filename(s) of pseudo-potential files.

    -h, --help
        display help information and exit.

"""

from bs4 import BeautifulSoup
import csv
import logging
logger = logging.getLogger(__name__)

def read_pseudo_cutoff(pp_file):
    try:
        with open(pp_file, "r") as fp:
            bs = BeautifulSoup(fp, "html.parser")
    except Exception as e:
        logger.error("{}: {}".format(pp_file, e))
        return None
        
    ecutwfc = None
    ecutrho = None

    if bs and bs.upf and bs.upf.pp_header:
        if bs.upf.pp_header.has_attr("wfc_cutoff"):
            ecutwfc = float(bs.upf.pp_header["wfc_cutoff"])
        if bs.upf.pp_header.has_attr("rho_cutoff"):
            ecutrho = float(bs.upf.pp_header["rho_cutoff"])
    if ecutwfc is None or ecutrho is None:
        if bs and bs.upf and bs.upf.pp_info:
            lines = str(bs.upf.pp_info).splitlines()
            sg = [s for s in lines if "Suggested" in s]
            for s in sg:
                if "wavefunctions" in s:
                    ecutwfc = float(s.split()[5])
                elif "charge density" in s:
                    ecutrho = float(s.split()[6])

    if ecutwfc is None or ecutrho is None:
        logger.warning("{}: wfc_cutoff or rho_cutoff not found".format(pp_file))
        return None
    else:
        return [ pp_file, ecutwfc, ecutrho ]

def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(prog="pp_cutoff")
    parser.add_argument("-f", "--filelist", default=None, metavar="list_file", help="list file that contains list of pseudo-potential files")
    parser.add_argument("-o", "--output", default=None, metavar="output_csv", help="output file to store the result in CSV foramt")
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

    tbl = [["pseudofile", "ecutwfc", "ecutrho" ]]
    for f in pp_files:
        x = read_pseudo_cutoff(f)
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
