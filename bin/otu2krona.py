#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

import numpy as np
from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep="\t"):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)

    fh.close()


def read_otu(file):

    data = {}
    samples = []
    n = 0

    for line in read_tsv(file, sep="\t"):
        n += 1
        if line[0].startswith("#") or n == 1:
            samples = line[1::]
            continue
        if "k__" not in line[0]:
            continue
        if line[0] not in data:
            data[line[0]] = np.array(line[1::]).astype(float)
        else:
            data[line[0]] = data[line[0]] + np.array(line[1::]).astype(float)

    return samples, data


def split_tax(tax):

    r = OrderedDict()

    for i in tax.split("|"):
        level, value = i.split("__", 1)
        if (level == "k") and (level in r):
            continue
        if "s" == level:
            value = value.split(".")[0]
        r[level] = value

    return r


def read_tax(tax):

    data = split_tax(tax)
    index = list(data.keys())
    r = []

    for i in  ["k", "p", "c", "o", "f", "g", "s"]:
        if i in data:
            r.append("%s__%s" % (i, data[i]))
            if index[-1] == i:
                break
        else:
            r.append("%s__unclassified" % i)

    return r


def otu2krona(file):

    samples, data = read_otu(file)

    for i in range(len(samples)):
        sample = samples[i]
        fo = open("%s.krona_report" % sample, "w")

        for j in data:
            line = data[j]
            abund = int(line[i])
            if abund <= 0:
                continue
            tax = read_tax(j)
            fo.write("%s\t%s\n" % (abund, "\t".join(tax)))
        fo.close()

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input annotated otu abundance table.")

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
     otu2krona.py: Generate the krona format file of each sample according to the otu abundance table.

attention:
     otu2krona.py meta.otu_tax.tsv

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    otu2krona(args.input)


if __name__ == "__main__":

    main()
