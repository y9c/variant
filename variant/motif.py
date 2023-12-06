#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-04-25 00:11

import gzip
import logging
import sys

import pyfaidx
import rich_click as click


def get_motif(chrom, pos, strand, fasta, lpad, rpad, to_upper=True):
    if not pos.isdecimal():
        logging.error(f"Position {pos} is not a number!")
        sys.exit(1)
    if strand not in ["+", "-"]:
        logging.error(f"Strand {strand} is not + or -!")
        sys.exit(1)

    # pos is 1-based, convert to 0-based
    pos = int(pos) - 1
    # Get the length of the chromosome
    chrom_len = len(fasta[chrom])

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
        sequence = "N" * lfill + fasta[chrom][start:end].seq + "N" * rfill
    else:
        sequence = (
            "N" * rfill
            + fasta[chrom][start:end].reverse.complement.seq
            + "N" * lfill
        )
    if to_upper:
        sequence = sequence.upper()

    return sequence


def _open_file(filename, mode="r"):
    if filename.endswith(".gz"):
        if mode == "w":
            return gzip.open(filename, "wt")
        return gzip.open(filename, "rt")
    else:
        if mode == "w":
            return click.open_file(filename, "w")
        return click.open_file(filename, "r")


def run_motif(input, output, fasta, lpad, rpad, with_header, columns):
    col_sep = "\t"
    columns_index = list(map(lambda x: int(x) - 1, columns.split(",")))
    columns_index_mapper = dict(zip(["chrom", "pos", "strand"], columns_index))
    strandness = "strand" in columns_index_mapper

    def parse_line(input_cols):
        m = get_motif(
            input_cols[columns_index_mapper["chrom"]],
            input_cols[columns_index_mapper["pos"]],
            input_cols[columns_index_mapper["strand"]] if strandness else "+",
            fasta_file,
            lpad,
            rpad,
        )

        output_cols = input_cols + [m]
        output_line = "\t".join(output_cols) + "\n"
        output_file.write(output_line)

    with _open_file(input) as input_file, _open_file(
        output, "w"
    ) as output_file, pyfaidx.Fasta(fasta) as fasta_file:
        # read first line
        input_cols = input_file.readline().strip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(input_cols) - 1:
            logging.error(f"Input file only have {len(input_cols)} columns!")
            sys.exit(1)
        if with_header:
            input_header = input_cols
        else:
            input_header = ["."] * len(input_cols)
            for n, i in columns_index_mapper.items():
                input_header[i] = n
        header_line = "\t".join(input_header + ["motif"]) + "\n"
        # output header column only if input file is with header
        if with_header:
            output_file.write(header_line)

        if not with_header:
            parse_line(input_cols)
        for line in input_file:
            input_cols = line.strip("\n").split(col_sep)
            parse_line(input_cols)
