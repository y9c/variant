#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-12-24 20:54


import logging
import os
import sys

import urllib3

from . import utils


def download_file(url, path):
    http = urllib3.PoolManager()
    r = http.request("GET", url, preload_content=False)
    if r.status != 200:
        raise ValueError("problem accessing " + url)
    with open(path, "wb") as f:
        for chunk in r.stream(1600):
            f.write(chunk)
    r.release_conn()
    http.clear()


def get_mapper(target, query, cache=None):
    # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt
    # 4 columns, first line starts with # is a comment
    # UCSC, ENSEMBL, GENEBANK, REFSEQ
    if cache is None:
        cache = os.path.expanduser("~/.cache/variant/coordinate")

    if not os.path.exists(cache):
        os.mkdir(cache)

    query = query[0].upper() + query[1:]
    target = target[0].lower() + target[1:]
    basename = "{}To{}.over.chain.gz".format(target, query)
    chain_path = os.path.join(cache, basename)

    if not os.path.exists(chain_path):
        url = (
            "https://hgdownload.cse.ucsc.edu/goldenPath/{}/liftOver/{}".format(
                target, basename
            )
        )
        download_file(url, chain_path)

    return (chain_path, target, query)


def run_coordinate(
    input, output, reference_mapping, buildin_mapping, with_header, columns
):
    col_sep = "\t"
    columns_index = list(map(lambda x: int(x) - 1, columns.split(",")))
    if len(columns_index) == 3:
        columns_index_mapper = dict(
            zip(["chrom", "pos", "strand"], columns_index)
        )
    elif len(columns_index) == 1:
        columns_index_mapper = dict(zip(["chrom"], columns_index))
    else:
        logging.error("Invalid columns!")
        sys.exit(1)

    chrom_col = columns_index_mapper.get("chrom")
    # pos_col = columns_index_mapper.get("pos")
    # strand_col = columns_index_mapper.get("strand")

    if reference_mapping:
        if buildin_mapping:
            logging.warning(
                "Both reference_mapping and buildin_mapping are provided, "
                "reference_mapping will be used!"
            )
        with utils.open_file(reference_mapping) as mapper_file:
            chrom_mapper = dict(
                (line.strip("\n").split("\t")[:2] for line in mapper_file)
            )
    elif buildin_mapping:
        if buildin_mapping == "U2E":
            chrom_mapper = dict(
                [("chr" + str(i), str(i)) for i in range(1, 100)]
                + [("chrX", "X"), ("chrY", "Y"), ("chrM", "MT")]
            )
        elif buildin_mapping == "E2U":
            chrom_mapper = dict(
                [(str(i), "chr" + str(i)) for i in range(1, 100)]
                + [("X", "chrX"), ("Y", "chrY"), ("MT", "chrM")]
            )
        else:
            logging.error("Invalid buildin_mapping!")
            sys.exit(1)
    else:
        logging.waring("No mapping provided!")
        chrom_mapper = {}

    with utils.open_file(input) as input_file, utils.open_file(
        output, "w"
    ) as output_file:

        def parse_line(input_cols):
            chrom = input_cols[chrom_col]
            # pos = input_cols[pos_col] if pos_col else None
            # strand = input_cols[strand_col] if strand_col else None
            chrom_rename = chrom_mapper.get(chrom, chrom)
            output_cols = (
                input_cols[:chrom_col]
                + [chrom_rename]
                + input_cols[chrom_col + 1 :]
            )
            output_line = col_sep.join(output_cols) + "\n"
            output_file.write(output_line)

        # read first line and check column number
        input_cols = input_file.readline().strip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(input_cols) - 1:
            logging.error(f"Input file only have {len(input_cols)} columns!")
            sys.exit(1)
        input_file.seek(0)

        if with_header:
            header_line = input_file.readline()
            output_file.write(header_line)
        for line in input_file:
            input_cols = line.strip("\n").split(col_sep)
            parse_line(input_cols)


if __name__ == "__main__":
    pass
