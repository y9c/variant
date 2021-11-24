#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-11-23 22:51

"""Annotate the mutation effect of a list of sites.

- 3 column input file: chromosome, position, reference allele and alternative allele.
- output file: gene name, transcript position, ...
"""


import argparse

import pyensembl
import varcode

IUPAC = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "U": ["T"],
    "T": ["T"],
    "M": ["A", "C"],
    "R": ["A", "G"],
    "W": ["A", "T"],
    "S": ["C", "G"],
    "Y": ["C", "T"],
    "K": ["G", "T"],
    "V": ["A", "C", "G"],
    "H": ["A", "C", "T"],
    "D": ["A", "G", "T"],
    "B": ["C", "G", "T"],
    "N": ["G", "A", "T", "C"],
}

COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "N": "N",
}


def expand_base(base):
    return IUPAC[base]


def reverse_base(base):
    return "".join([COMPLEMENT[b] for b in base][::-1])


def site2mut(chrom, pos, strand, ref, alt, genome, mode="DNA"):
    pad = 10

    chrom = str(chrom).replace("chr", "")
    pos = int(pos)
    ref = ref.upper()
    alt = alt.upper()
    if chrom not in genome.contigs():
        return [None] * 9
    effs = varcode.EffectCollection(
        [
            e
            for b in expand_base(alt)
            if b != ref
            for e in varcode.Variant(
                contig=chrom, start=pos, ref=ref, alt=b, ensembl=genome
            ).effects()
        ]
    )

    # filter the correct strand
    if mode == "RNA" and strand in "+-":
        eff_list = []
        for e in effs:
            if type(e).__name__ == "Intergenic":
                eff_list.append(e)
            elif e.transcript.strand == strand:
                eff_list.append(e)
        if len(eff_list) == 0:
            return [None] * 9
        effs = varcode.EffectCollection(eff_list)

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
    parser.add_argument("--rna", action="store_true")
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
        input_header = [
            "chrom",
            "pos",
            "strand",
            "ref",
            "alt",
        ]
        annot_header = [
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
        print("#" + "\t".join(input_header + annot_header), file=output_f)
        for l in f:
            c, p, s, ref, alt, *_ = l.strip("\n").split("\t")
            mode = "RNA" if args.rna else "DNA"
            if mode == "RNA" and s == "-":
                ref = reverse_base(ref)
                alt = reverse_base(alt)
            annot = site2mut(c, p, s, ref, alt, ensembl_genome, mode=mode)
            annot_str = list(map(str, annot))
            print("\t".join([c, p, s, ref, alt] + annot_str), file=output_f)
