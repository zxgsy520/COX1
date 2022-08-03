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


def read_cox1_tax(file):

    data = {}

    for line in read_tsv(file, sep="\t"):
        if "#" in line[0]:
            continue

        if line[1] not in data:
            data[line[1]] = 0
        data[line[1]] += 1

    return data


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


def cox1_tax2krona(file):

    data = read_cox1_tax(file)

    for i in data:
        abund = int(data[i])
        if abund <= 0:
            continue
        tax = read_tax(i)

        print("%s\t%s" % (abund, "\t".join(tax)))

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
     cox1_tax2krona.py: Generate the krona format file of each sample according to the otu abundance table.

attention:
     cox1_tax2krona.py cox1.taxonomy.tsv.gz

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    cox1_tax2krona(args.input)


if __name__ == "__main__":

    main()
