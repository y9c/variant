#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-12-24 20:54


import os
import sys

import urllib3

from . import utils

LOGGER = utils.get_logger(__name__)


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


def get_mapper(reference, mapper_type, cache=None):
    if reference == "hg38":
        names = ["ucsc", "assembly", "ensembl", "genbank", "refseq"]
    elif reference == "mm39":
        names = ["ucsc", "assembly", "ensembl", "genbank", "refseq"]
    # elif reference == "hg19":
    #     names = ["ucsc", "assembly", "genbank", "refseq"]
    # elif reference == "mm10":
    #     names = ["uscs", "refseq", "genbank", "assembly"]
    else:
        LOGGER.error("Invalid reference!")
        sys.exit(1)
    chrom_mapper = {}

    if cache is None:
        cache = os.path.expanduser("~/.cache/variant/coordinate")

    if not os.path.exists(cache):
        os.makedirs(cache, exist_ok=True)

    basename = "{}.chromAlias.txt".format(reference)
    reference_path = os.path.join(cache, basename)
    if not os.path.exists(reference_path):
        url = (
            "https://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/{}".format(
                reference,
                basename if reference != "hg38" else "latest/" + basename,
            )
        )
        download_file(url, reference_path)

    with utils.open_file(reference_path) as mapper_file:
        # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt
        # # sequenceName	alias names	UCSC database: hg38
        # chr1	1	CM000663.2	NC_000001.11
        # new version: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chromAlias.txt
        # # ucsc	assembly	ensembl	genbank	refseq
        # chr1	1	1	CM000663.2	NC_000001.11
        #
        # https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromAlias.txt
        # # ucsc	assembly	ensembl	genbank	refseq
        # chr1	1	1	CM000994.3	NC_000067.7
        #
        # https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chromAlias.txt
        # ucsc	assembly	genbank	refseq
        # chr1	1	CM000663.1	NC_000001.10
        #
        # https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chromAlias.txt
        # uscs	refseq	genbank	assembly
        # sequenceName	alias names	UCSC database: mm10
        # chr1	NC_000067.6	CM000994.2	1

        for line in mapper_file:
            if line.startswith("#"):
                continue
            cols = line.strip("\n").split("\t")
            mapper = dict(zip(names, cols))
            if mapper_type == "U2E":
                chrom_mapper[mapper.get("ucsc", "")] = mapper.get(
                    "ensembl", ""
                )
            elif mapper_type == "E2U":
                chrom_mapper[mapper.get("ensembl", "")] = mapper.get(
                    "ucsc", ""
                )
    return chrom_mapper


def run_coordinate(
    input,
    output,
    reference_mapping,
    buildin_mapping,
    columns,
    with_header,
    keep_original,
):
    col_sep = "\t"
    columns_index = list(map(lambda x: int(x) - 1, columns.split(",")))
    if len(columns_index) <= 3:
        columns_index_mapper = dict(
            zip(["chrom", "pos", "strand"], columns_index)
        )
    else:
        LOGGER.error("Invalid number of columns!")
        sys.exit(1)

    chrom_col = columns_index_mapper.get("chrom")
    # pos_col = columns_index_mapper.get("pos")
    # strand_col = columns_index_mapper.get("strand")

    if reference_mapping:
        if buildin_mapping:
            LOGGER.warning(
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
        elif buildin_mapping in [
            "U2E-hg38",
            "E2U-hg38",
            "U2E-mm39",
            "E2U-mm39",
        ]:
            chrom_mapper = get_mapper(
                buildin_mapping.split("-")[1], buildin_mapping.split("-")[0]
            )
        else:
            LOGGER.error("Invalid buildin_mapping!")
            sys.exit(1)
    else:
        LOGGER.warning("No mapping provided!")
        chrom_mapper = {}

    with utils.open_file(input) as input_file, utils.open_file(
        output, "w"
    ) as output_file:

        def parse_line(input_cols):
            chrom = input_cols[chrom_col]
            # pos = input_cols[pos_col] if pos_col else None
            # strand = input_cols[strand_col] if strand_col else None
            chrom_rename = chrom_mapper.get(chrom, chrom)
            if not keep_original:
                output_cols = (
                    input_cols[:chrom_col]
                    + [chrom_rename]
                    + input_cols[chrom_col + 1 :]
                )
            else:
                output_cols = input_cols + [chrom_rename]
            output_line = col_sep.join(output_cols) + "\n"
            output_file.write(output_line)

        # read first line and check column number
        input_cols = input_file.readline().strip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(input_cols) - 1:
            LOGGER.error(f"Input file only have {len(input_cols)} columns!")
            sys.exit(1)
        input_file.seek(0)

        if with_header:
            header_line = input_file.readline()
            if keep_original:
                header_line = (
                    header_line.strip("\n") + col_sep + "RenamedChrom" + "\n"
                )
            output_file.write(header_line)
        for line in input_file:
            input_cols = line.strip("\n").split(col_sep)
            parse_line(input_cols)


if __name__ == "__main__":
    pass
