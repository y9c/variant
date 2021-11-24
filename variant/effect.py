#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-11-23 22:51

"""Annotate the mutation effect of a list of sites.

- 3 column input file: chromosome, position and alternative allele.
- output file: gene name, transcript position, ...
"""


import argparse

import pyensembl
import varcode


def site2mut(chrom, pos, base, genome):
    pad = 10

    chrom = str(chrom)
    base = base.upper()
    if chrom not in genome.contigs():
        return [None] * 9
    effs = varcode.EffectCollection(
        [
            e
            for b in "ATGC"
            if b != base
            for e in varcode.Variant(
                contig=chrom, start=pos, ref=base, alt=b, ensembl=genome
            ).effects()
        ]
    )

    # show top effect only
    eff = effs.top_priority_effect()
    mut_type = type(eff).__name__
    transcript_id = None if mut_type == "Intergenic" else eff.transcript_id

    non_exon_recored = [mut_type, eff.gene_name, transcript_id] + [None] * 6
    if mut_type in ["Intergenic", "Intronic", "SpliceDonor"]:
        return non_exon_recored
    try:
        transcript_pos = eff.transcript.spliced_offset(pos) + 1
    except:
        return non_exon_recored
    # transcript motif
    s = eff.transcript.sequence
    s5 = (s[max(transcript_pos - 1 - pad, 0) : transcript_pos - 1]).rjust(
        pad, "N"
    )
    s0 = s[transcript_pos - 1]
    s3 = (s[transcript_pos : transcript_pos + pad]).ljust(pad, "N")
    transcript_motif = s5 + s0 + s3
    # codon
    if (
        mut_type not in ["NoncodingTranscript", "IncompleteTranscript"]
        and eff.transcript.first_start_codon_spliced_offset
        <= transcript_pos - 1
        <= eff.transcript.last_stop_codon_spliced_offset
    ):
        coding_pos = (
            transcript_pos - eff.transcript.first_start_codon_spliced_offset
        )
        codon_start = (coding_pos - 1) // 3 * 3
        codon_ref = eff.transcript.coding_sequence[
            codon_start : codon_start + 3
        ]
    else:
        coding_pos = None
        codon_ref = None

    if mut_type == "Silent":
        aa_pos_offset = eff.aa_pos
        aa_ref = eff.aa_ref if coding_pos else None
    elif mut_type in ["IntronicSpliceSite", "ExonicSpliceSite"]:
        aa_pos_offset = eff.alternate_effect.aa_pos
        aa_ref = eff.alternate_effect.aa_ref if coding_pos else None
    elif mut_type == "StopLoss":
        aa_pos_offset = eff.aa_mutation_start_offset
        aa_ref = "*"
    else:
        aa_pos_offset = eff.aa_mutation_start_offset
        aa_ref = eff.aa_ref if coding_pos else None
    if aa_ref == "":
        aa_ref = None
    aa_pos = int(aa_pos_offset) + 1 if coding_pos else None

    return [
        mut_type,
        eff.gene_name,
        transcript_id,
        transcript_pos,
        transcript_motif,
        coding_pos,
        codon_ref,
        aa_pos,
        aa_ref,
    ]


def usage():
    print("variant-effect -i <input> [-r <ref> -o <output>]")


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="input file")
    parser.add_argument(
        "-o", "--output", default="-", help="output file (default: stdout)"
    )
    parser.add_argument(
        "-r",
        "--reference",
        default="homo_sapiens",
        help="reference species (default: human)",
    )
    args = parser.parse_args()
    input_file = args.input
    if args.output != "-":
        output_f = open(args.output, "w")
    else:
        output_f = None

    ensembl_genome = pyensembl.EnsemblRelease(
        release="104", species=args.reference
    )
    try:
        ensembl_genome.index()
    except:
        ensembl_genome.download()
        ensembl_genome.index()
    with open(input_file, "r") as f:
        header = [
            "ref",
            "pos",
            "alt",
            "mut_type",
            "gene_name",
            "transcript_id",
            "transcript_pos",
            "transcript_motif",
            "coding_pos",
            "codon_ref",
            "aa_pos",
            "aa_ref",
        ]
        print("#" + "\t".join(header), file=output_f)
        for l in f:
            c, p, b, *_ = l.strip("\n").split("\t")
            annot = list(
                map(
                    str,
                    site2mut(c.replace("chr", ""), int(p), b, ensembl_genome),
                )
            )
            print("\t".join([c, p, b] + annot), file=output_f)
