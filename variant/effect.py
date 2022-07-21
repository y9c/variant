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


import sys

import click
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


def mut2eff(chrom, pos, strand, ref, alt, genome, biotype):

    chrom = str(chrom).replace("chr", "")
    pos = int(pos)
    ref = ref.upper()
    alt = alt.upper()
    if chrom not in genome.contigs():
        return None
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
    if biotype == "RNA" and strand in "+-":
        eff_list = []
        for e in effs:
            if type(e).__name__ == "Intergenic":
                eff_list.append(e)
            elif e.transcript is None:
                eff_list.append(e)
            elif e.transcript.strand == strand:
                eff_list.append(e)
        if len(eff_list) == 0:
            return None
        effs = varcode.EffectCollection(eff_list)

    return effs


def parse_eff(eff, pos, pad):
    pos = int(pos)
    mut_type = type(eff).__name__
    d2s = []
    if eff.transcript is None:
        distance2splice = None
    else:
        for i, exon in enumerate(eff.transcript.exons):
            if i == 0:
                if exon.strand == "+":
                    d2s.append(pos - exon.end)
                else:
                    d2s.append(exon.start - pos)
            elif i == len(eff.transcript.exons) - 1:
                if exon.strand == "+":
                    d2s.append(pos - exon.start)
                else:
                    d2s.append(exon.end - pos)
            else:
                if exon.strand == "+":
                    d2s.append(pos - exon.start)
                    d2s.append(pos - exon.end)
                else:
                    d2s.append(exon.start - pos)
                    d2s.append(exon.end - pos)
        # minus distance to splice site
        distance2splice = sorted(d2s, key=lambda x: abs(x))[0]

    transcript_id = None if mut_type == "Intergenic" else eff.transcript_id

    non_exon_recored = [mut_type, eff.gene_name, transcript_id] + [None] * 7
    if mut_type in ["Intergenic", "Intronic", "SpliceDonor"]:
        return non_exon_recored
    try:
        transcript_pos = eff.transcript.spliced_offset(pos) + 1
    except:
        return non_exon_recored
    # transcript motif
    s = eff.transcript.sequence
    if s is None:
        transcript_motif = None
    else:
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
        aa_pos = eff.aa_pos
        aa_ref = eff.aa_ref if coding_pos else None
    elif mut_type in ["IntronicSpliceSite", "ExonicSpliceSite"]:
        aa_pos = eff.alternate_effect.aa_pos
        aa_ref = eff.alternate_effect.aa_ref if coding_pos else None
    elif mut_type == "StopLoss":
        aa_pos = eff.aa_mutation_start_offset + 1
        aa_ref = "*"
    else:
        aa_pos = eff.aa_mutation_start_offset + 1
        aa_ref = eff.aa_ref if coding_pos else None
    if aa_ref == "":
        aa_pos = None
        aa_ref = None

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
        distance2splice,
    ]


def site2mut(chrom, pos, strand, ref, alt, genome, biotype, pad=10, all=False):
    effs = mut2eff(chrom, pos, strand, ref, alt, genome, biotype)
    if effs is None:
        return [[None] * 9]
    if not all:
        eff = effs.top_priority_effect()
        return [parse_eff(eff, pos, pad)]
    return [parse_eff(eff, pos, pad) for eff in effs]


@click.command(
    help="Variant (genomic variant analysis in python)", no_args_is_help=True
)
@click.option(
    "--input", "-i", "input", help="Input position file.", required=True
)
@click.option(
    "--output",
    "-o",
    "output",
    default="-",
    help="Output annotation file",
    required=False,
)
@click.option(
    "--reference",
    "-r",
    "reference",
    default="homo_sapiens",
    help="reference species",
    required=False,
)
@click.option(
    "--type",
    "-t",
    "biotype",
    type=click.Choice(["DNA", "RNA"], case_sensitive=False),
    default="DNA",
)
@click.option(
    "--npad",
    "-n",
    "npad",
    default=10,
    type=int,
    help="Number of padding base to call motif.",
)
@click.option("--all", help="Output all effects.", is_flag=True)
@click.option("--header", help="With header line", is_flag=True)
@click.option(
    "--columns",
    "-c",
    "columns",
    default=[1, 2, 3, 4, 5],
    show_default=True,
    type=int,
    help="Sets columns for site info. (Chrom,Pos,Strand,Ref,Alt)",
    multiple=True,
)
def run(input, output, reference, biotype, npad, all, header, columns):
    ensembl_genome = pyensembl.EnsemblRelease(release="104", species=reference)
    try:
        ensembl_genome.index()
    except:
        ensembl_genome.download()
        ensembl_genome.index()

    with open(input, "r") as f:
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
            "distance2splice",
        ]
        if header:
            input_header = f.readline().strip().split()
        else:
            input_header = ["chrom", "pos", "strand", "ref", "alt"]

        if output == "-":
            output_file = sys.stdout
        else:
            output_file = open(output, "w")
        print("#" + "\t".join(input_header + annot_header), file=output_file)
        for l in f:
            input_cols = l.strip("\n").split("\t")
            c = input_cols[columns[0] - 1]
            p = input_cols[columns[1] - 1]
            s = input_cols[columns[2] - 1]
            ref = input_cols[columns[3] - 1]
            alt = input_cols[columns[4] - 1]
            if biotype == "RNA" and s == "-":
                ref = reverse_base(ref.upper())
                alt = reverse_base(alt.upper())
            annot_list = site2mut(
                c, p, s, ref, alt, ensembl_genome, biotype, pad=npad, all=all
            )
            for annot in annot_list:
                annot_str = list(map(str, annot))
                print(
                    "\t".join([c, p, s, ref, alt] + annot_str),
                    file=output_file,
                )
