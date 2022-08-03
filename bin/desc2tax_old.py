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
    taxids = set()
    genus = set()
    species = set()

    for line in read_tsv(file, "\t"):
        if len(line) <= 3:
            r[line[0]] = line[1::]
            if line[-1] != "-":
                species.add("s__%s" % line[-1])
        else:
            r[line[0]] = line[1:3]
            genus.add("g__%s" % line[-1])
            r[line[0]].append(line[-1])

        taxids.add(line[1])

    return r, taxids, genus, species


def desc2tax(desc, taxonomy):

    data, taxids, genus, species = read_gene_desc(desc)

    r = {}
    for line in read_tsv(taxonomy, "\t"):
        if line[0] not in taxids:
            continue
        r[line[0]] = line[1].split(".")[0]

    dgs = {}
    for i in r:
        for j in species:
            if j in r[i]:
                dgs[j] = r[i]
                break
        for j in genus:
            if j in r[i]:
                dgs[j] = r[i].split("|s__")[0]
                break

    print("#protein_id\ttax\ttax_id\torganism")
    for i in data:
        taxid = data[i][0]
        skey = data[i][-1]

        tax = "k__Eukaryota"
        if taxid in r:
            tax = r[taxid]
        elif skey != "-":
            if skey in dgs:
                tax = dgs[skey]
        else:
            pass
        if ("g__" not in tax) and len(data[i])==3:
            if skey != "-":
                tax = "%s|g__%s" % (tax, skey)
        if "s__" not in tax:
            if data[i][1] != "-":
                tax = "%s|s__%s" % (tax, data[i][1])

        print("%s\t%s\t%s\t%s" % (i, tax, taxid, data[i][1]))

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input sequence and taxid corresponding list.")
    parser.add_argument("-tax", "--taxonomy", metavar="FILE", type=str, required=True,
        help="input taxonomy file(kraken.taxonomy.gz).")

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
     desc2tax.py: Get the taxonomy of a sequence

attention:
     desc2tax.py cox1.describe.tsv --taxonomy kraken.taxonomy.gz >cox1.taxonomy.tsv

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    desc2tax(args.input, args.taxonomy)


if __name__ == "__main__":

    main()
