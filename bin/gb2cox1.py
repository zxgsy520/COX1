#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

from Bio import SeqIO

LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def get_cox1(file):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    fo = open("cox1.pep.fasta", "w")

    for record in SeqIO.parse(fh, "genbank"):
        organism = record.annotations["organism"]

        for genes in record.features:
            if genes.type != "CDS":
                continue

            gene_desc = genes.qualifiers
            if "gene" in gene_desc:
                gene_id = gene_desc["gene"][0]
            elif "gene_synonym" in gene_desc:
                gene_id = gene_desc["gene_synonym"][0]
            elif "locus_tag" in gene_desc:
                gene_id = gene_desc["locus_tag"][0]
            else:
                gene_id = ""

            if gene_id not in ["COX1", "cox1", "COI"]:
                continue

            if "coded_by" in gene_desc:
                ncid = gene_desc["coded_by"][0].split(":")[0]
            elif "protein_id" in gene_desc:
                ncid = gene_desc["protein_id"][0]
            else:
                ncid = ""

            seq = genes.extract(record.seq)
            print(">%s|%s [%s]\n%s" % (ncid, gene_id.upper(), organism, seq))

            if "translation" in gene_desc:
                protein = str(gene_desc["translation"][0])
                fo.write(">%s|%s [%s]\n%s\n" % (ncid, gene_id.upper(), organism, protein))

    return 0


def gb2cox1(files):

    for file in files:
        get_cox1(file)

    return 0


def add_help_args(parser):

    parser.add_argument("genbank", nargs="+", metavar="FILE", type=str,
        help="Input the mitochondrial genebank file.")

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
     gb2cox1.py Extract COX1 gene according to genbank file.

attention:
     gb2cox1.py mitochondrion.1.*.gbff.gz >COX1.fa

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    gb2cox1(args.genbank)


if __name__ == "__main__":

    main()
