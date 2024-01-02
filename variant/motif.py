#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-04-25 00:11

import sys

import pyfaidx

from . import utils

LOGGER = utils.get_logger(__name__)


def get_motif(chrom_obj, chrom_len, pos, strand, lpad, rpad):
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
        rfill = rpad - (chrom_len - pos)
    else:
        start = 0
        end = chrom_len
        lfill = lpad
        rfill = rpad - (chrom_len - pos)

    if strand == "+":
        sequence = "N" * lfill + chrom_obj[start:end].seq + "N" * rfill
    else:
        sequence = (
            "N" * rfill
            + chrom_obj[start:end].reverse.complement.seq
            + "N" * lfill
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

    with utils.open_file(input) as input_file, utils.open_file(
        output, "w"
    ) as output_file, pyfaidx.Fasta(fasta) as fasta_file:
        chrom_len_mapper = {k: len(v) for k, v in fasta_file.items()}

        def parse_line(input_cols):
            chrom_name = input_cols[columns_index_mapper["chrom"]]
            chrom_obj = fasta_file[chrom_name]
            chrom_len = chrom_len_mapper[chrom_name]

            strand = (
                input_cols[columns_index_mapper["strand"]]
                if strandness
                else "+"
            )
            m = get_motif(
                chrom_obj,
                chrom_len,
                input_cols[columns_index_mapper["pos"]],
                strand,
                lpad,
                rpad,
            )
            if to_upper:
                m = m.upper()
            if wrap_site:
                if strand == "+":
                    m = m[:lpad] + "[" + m[lpad] + "]" + m[lpad + 1 :]
                else:
                    m = m[:rpad] + "[" + m[rpad] + "]" + m[rpad + 1 :]

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
        header_line = col_sep.join(input_header + ["motif"]) + "\n"
        # output header column only if input file is with header
        if with_header:
            output_file.write(header_line)

        if not with_header:
            parse_line(input_cols)
        for line in input_file:
            input_cols = line.strip("\n").split(col_sep)
            parse_line(input_cols)
