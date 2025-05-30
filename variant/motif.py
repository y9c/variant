#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-04-25 00:11

import sys

import pysam
from xopen import xopen

from . import utils
from .seqpy import revcomp

LOGGER = utils.get_logger(__name__)


def get_motif(fasta_file, chrom, pos, strand, lpad, rpad):
    chrom_len = fasta_file.get_reference_length(chrom)
    if not pos.isdecimal():
        LOGGER.error(f"Position {pos} is not a number!")
        sys.exit(1)
    if strand not in ["+", "-"]:
        LOGGER.error(f"Strand {strand} is not + or -!")
        sys.exit(1)

    # pos is 1-based, convert to 0-based
    pos = int(pos) - 1

    # Get the sequence of the chromosome at the given position
    if pos - lpad >= 0 and pos + rpad < chrom_len:
        start = pos - lpad
        end = pos + rpad + 1
        lfill = 0
        rfill = 0
    elif pos - lpad < 0 and pos + rpad < chrom_len:
        start = 0
        end = pos + rpad + 1
        lfill = lpad - pos
        rfill = 0
    elif pos - lpad >= 0 and pos + rpad >= chrom_len:
        start = pos - lpad
        end = chrom_len
        lfill = 0
        rfill = rpad - (chrom_len - pos) + 1
    else:
        start = 0
        end = chrom_len
        lfill = lpad
        rfill = rpad - (chrom_len - pos) + 1

    if strand == "+":
        sequence = "N" * lfill + fasta_file.fetch(chrom, start, end) + "N" * rfill
    else:
        sequence = (
            "N" * rfill + revcomp(fasta_file.fetch(chrom, start, end)) + "N" * lfill
        )

    return sequence


def run_motif(
    input,
    output,
    fasta,
    lpad,
    rpad,
    with_header,
    columns,
    to_upper=True,
    wrap_site=True,
):
    col_sep = "\t"
    columns_index = list(map(lambda x: int(x) - 1, columns.split(",")))
    columns_index_mapper = dict(zip(["chrom", "pos", "strand"], columns_index))
    strandness = "strand" in columns_index_mapper

    with xopen(input) as input_file, xopen(output, "w") as output_file, pysam.FastaFile(
        fasta
    ) as fasta_file:

        def parse_line(input_cols):
            chrom = input_cols[columns_index_mapper["chrom"]]
            strand = input_cols[columns_index_mapper["strand"]] if strandness else "+"
            pos = input_cols[columns_index_mapper["pos"]]

            if strand == "+":
                m = get_motif(fasta_file, chrom, pos, strand, lpad, rpad)
            else:
                m = get_motif(fasta_file, chrom, pos, strand, rpad, lpad)
            if to_upper:
                m = m.upper()
            if wrap_site:
                m = m[:lpad] + "[" + m[lpad] + "]" + m[lpad + 1 :]
                # if strand == "+":
                #     m = m[:lpad] + "[" + m[lpad] + "]" + m[lpad + 1 :]
                # else:
                #     m = m[:rpad] + "[" + m[rpad] + "]" + m[rpad + 1 :]

            output_cols = input_cols + [m]
            output_line = col_sep.join(output_cols) + "\n"
            output_file.write(output_line)

        # read first line
        input_cols = input_file.readline().strip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(input_cols) - 1:
            LOGGER.error(f"Input file only have {len(input_cols)} columns!")
            sys.exit(1)
        if with_header:
            input_header = input_cols
        else:
            input_header = ["."] * len(input_cols)
            for n, i in columns_index_mapper.items():
                input_header[i] = n
        # header_line = col_sep.join(input_header + ["motif"]) + "\n"
        header_line = col_sep.join(input_header + ["motif"]) + "\n"
        # output header column only if input file is with header
        if with_header:
            output_file.write(header_line)

        if not with_header:
            parse_line(input_cols)
        for line in input_file:
            input_cols = line.strip("\n").split(col_sep)
            parse_line(input_cols)
