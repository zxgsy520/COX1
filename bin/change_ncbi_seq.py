#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse
import collections

from Bio import SeqIO


LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_protein_result(file):

    r = {}
    organism = ""
    protein_id = ""
    for line in open(file):
        line = line.strip()

        if not line:
            continue
        if "[" in line:
            organism = line.split("[", 1)[1].strip("]")
        if "GI:" in line:
            protein_id = line.split()[0]
            r[protein_id] = organism

    return r


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line)
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


def split_attr(attributes):

    r = collections.OrderedDict()
    contents = attributes.strip("[").strip("]").split("] [")

    for content in contents:
        content = content
        if not content:
            continue
        if "=" not in content:
            print("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def change_ncbi_seq(file, presult):

    data = read_protein_result(presult)

    fo = open("gene.describe.tsv", "w")
    fo.write("#protein_id\ttax_id\torganism\n")
    for seqid, seq in read_fasta(file):
        seqid, attribute = seqid.split(" ", 1)
        attribute = split_attr(attribute)
        protein_id = attribute["protein_id"]

        if "gene" in  attribute:
            gene = attribute["gene"]
        elif "protein" in  attribute:
            gene = attribute["protein"]
        else:
            gene = ""

        if protein_id in data:
            print(">%s|%s [organism=%s]\n%s" % (protein_id, gene.upper(), data[protein_id], seq))
            fo.write("%s\t-\t%s\n" % (protein_id, data[protein_id]))
        else:
            print(">%s|%s/n%s" % (protein_id, gene.upper(), seq))
            fo.write("%s\t-\t-\n" % (protein_id))

    fo.close()
    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input the downloaded CDS sequence.")
    parser.add_argument("-pr", "--presult", metavar="FILE", type=str, required=True,
        help="Input protein sequence description file(protein_result.txt).")

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
     change_ncbi_seq.py: Change sequence file.

attention:
     change_ncbi_seq.py all.cox1.fa -pr protein_result.txt >all.new_cox1.fa

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    change_ncbi_seq(args.input, args.presult)


if __name__ == "__main__":

    main()
