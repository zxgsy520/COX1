#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fh.close()


def read_gene_desc(file):

    r = {}
    for line in read_tsv(file, "\t"):
        if len(line) <=3:
            r[line[0]] = line[1::]
        else:
            r[line[0]] = line[1:3]

    return r


def protid2taxid(desc, prot2taxid):

    data = read_gene_desc(desc)

    for line in read_tsv(prot2taxid, "\t"):
        if line[1] in data:
            data[line[1]][0] = line[2]
        elif line[0] in data:
            data[line[0]][0] = line[2]
        else:
            continue

    print("#protein_id\ttax_id\torganism")
    for i in data:
        print("%s\t%s\t%s" % (i, data[i][0], data[i][1]))

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input protein id.")
    parser.add_argument("-pt", "--prot2taxid", metavar="FILE", type=str, required=True,
        help="Input the taxid corresponding to the protein id(prot.accession2taxid.gz).")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
name:
     protid2taxid.py: Taxid of extracted protein sequence.

attention:
     protid2taxid.py gene.describe.tsv -pt prot.accession2taxid.gz >gene.new_describe.tsv

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    protid2taxid(args.input, args.prot2taxid)


if __name__ == "__main__":

    main()
