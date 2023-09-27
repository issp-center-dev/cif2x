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

import xml.etree.ElementTree as ET
import csv
import logging
logger = logging.getLogger(__name__)

def read_pseudo_cutoff(pp_file):
    try:
        tree = ET.parse(pp_file)
    except Exception as e:
        logger.error("{}: {}".format(pp_file, e))
        return None
        
    root = tree.getroot()

    wfc_cutoff = None
    rho_cutoff = None

    header = root.find("PP_HEADER")
    if header is not None:
        attr = header.attrib
        if "wfc_cutoff" in attr:
            wfc_cutoff = float(attr["wfc_cutoff"])
        if "rho_cutoff" in attr:
            rho_cutoff = float(attr["rho_cutoff"])

    if wfc_cutoff is None or rho_cutoff is None:
        logger.warning("{}: wfc_cutoff or rho_cutoff not found".format(pp_file))
        return None
    else:
        return [ pp_file, wfc_cutoff, rho_cutoff ]

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

    tbl = []
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
