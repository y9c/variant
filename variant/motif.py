#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-04-25 00:11

import gzip

import pyfaidx
import rich_click as click


def get_motif(chrom, pos, strand, fasta, n, to_upper=True):
    pos = int(pos)
    # Get the length of the chromosome
    chrom_len = len(fasta[chrom])

    # Get the sequence of the chromosome at the given position
    start = max(pos - n - 1, 0)
    end = pos + n

    if strand == "+":
        sequence = fasta[chrom][start:end].seq
    else:
        sequence = fasta[chrom][start:end].reverse.complement.seq

    if to_upper:
        sequence = sequence.upper()

    # Fill with "N" if the border reaches the end or start of the sequence
    if start == 0:
        sequence = (
            "N" * (n - (pos - 1))
            if strand == "+"
            else sequence[::-1] + "N" * (n - (pos - 1))
        )
    if end >= chrom_len:
        sequence = (
            sequence + "N" * (n - (chrom_len - pos + 1))
            if strand == "+"
            else "N" * (n - (chrom_len - pos)) + sequence[::-1]
        )

    # If the position is at the end of the reference sequence, fill with "N"
    if len(sequence) < 2 * n + 1:
        sequence = (
            sequence.ljust(2 * n + 1, "N")
            if strand == "+"
            else sequence.rjust(2 * n + 1, "N")
        )
    return sequence


def _open_file(filename, mode="r"):
    if filename.endswith(".gz"):
        if mode == "w":
            return gzip.open(filename, "wb")
        return gzip.open(filename, "rb")
    else:
        if mode == "w":
            return click.open_file(filename, "w")
        return click.open_file(filename, "r")


def run_motif(input, output, fasta, npad, with_header, columns):
    col_sep = "\t"
    columns_index = list(map(lambda x: int(x) - 1, columns.split(",")))
    columns_index_mapper = dict(zip(["chrom", "pos", "strand"], columns_index))
    if "strand" in columns_index_mapper:
        strandness = True

    def parse_line(input_cols):
        m = get_motif(
            input_cols[columns_index_mapper["chrom"]],
            input_cols[columns_index_mapper["pos"]],
            input_cols[columns_index_mapper["strand"]] if strandness else "+",
            fasta_file,
            npad,
        )

        output_cols = input_cols + [m]
        output_line = "\t".join(output_cols) + "\n"
        if output.endswith(".gz"):
            output_line = output_line.encode()
        output_file.write(output_line)

    with _open_file(input) as input_file, _open_file(
        output, "w"
    ) as output_file, pyfaidx.Fasta(fasta) as fasta_file:
        if with_header:
            if input.endswith(".gz"):
                input_header = (
                    input_file.readline().decode().strip().split(col_sep)
                )
            else:
                input_header = input_file.readline().strip().split(col_sep)
        else:
            if input.endswith(".gz"):
                input_cols = (
                    input_file.readline().decode().strip().split(col_sep)
                )
            else:
                input_cols = input_file.readline().strip().split(col_sep)

            input_header = ["."] * len(input_cols)
            # rename header
            for n, i in columns_index_mapper.items():
                input_header[i] = n
        header_line = "\t".join(input_header + ["motif"]) + "\n"
        if output.endswith(".gz"):
            header_line = header_line.encode()
        if with_header:
            output_file.write(header_line)

        if not with_header:
            parse_line(input_cols)
        for line in input_file:
            if input.endswith(".gz"):
                line = line.decode()
            input_cols = line.strip().split(col_sep)
            parse_line(input_cols)
