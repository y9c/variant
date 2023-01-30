#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-11-23 22:51

"""Annotate the mutation effect of a list of sites.

Update in 2023-01-28: Now only chrom and pos is required, other columns is optional


Effect type  | Description
-----------: | :-----------
*AlternateStartCodon* | Replace annotated start codon with alternative  start codon (*e.g.* "ATG>CAG").
*ComplexSubstitution* | Insertion and deletion of multiple amino acids.
*Deletion* | Coding mutation which causes deletion of amino acid(s).
*ExonLoss* | Deletion of entire exon, significantly disrupts protein.
*ExonicSpliceSite* | Mutation at the beginning or end of an exon, may affect splicing.
*FivePrimeUTR* | Variant affects 5' untranslated region before start codon.
*FrameShiftTruncation* | A frameshift which leads immediately to a stop codon (no novel amino acids created).
*FrameShift* | Out-of-frame insertion or deletion of nucleotides, causes novel protein sequence and often premature stop codon.
*IncompleteTranscript* | Can't determine effect since transcript annotation is incomplete (often missing either the start or stop codon).
*Insertion* | Coding mutation which causes insertion of amino acid(s).
*Intergenic* | Occurs outside of any annotated gene.
*Intragenic* |Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA.
*IntronicSpliceSite* | Mutation near the beginning or end of an intron but less likely to affect splicing than donor/acceptor mutations.
*Intronic* | Variant occurs between exons and is unlikely to affect splicing.
*NoncodingTranscript* | Transcript doesn't code for a protein.
*PrematureStop* | Insertion of stop codon, truncates protein.
*Silent* | Mutation in coding sequence which does not change the amino acid sequence of the translated protein.
*SpliceAcceptor* | Mutation in the last two nucleotides of an intron, likely to affect splicing.
*SpliceDonor* | Mutation in the first two nucleotides of an intron, likely to affect splicing.
*StartLoss* | Mutation causes loss of start codon, likely result is that an alternate start codon will be used down-stream (possibly in a different frame).
*StopLoss* | Loss of stop codon, causes extension of protein by translation of nucleotides from 3' UTR.
*Substitution* | Coding mutation which causes simple substitution of one amino acid for another.
*ThreePrimeUTR* | Variant affects 3' untranslated region after stop codon of mRNA.
"""


import sys
from dataclasses import dataclass

import pyensembl
import rich_click as click
import varcode
from varcode.effects import Intergenic

from . import effect_ordering

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
    ".": ["G", "A", "T", "C"],
    "-": ["G", "A", "T", "C"],
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
    ".": ".",
    "-": "-",
}

FEATURE_MAPPER = {
    "AlternateStartCodon": "StartCodon",
    "StartLoss": "StartCodon",
    "StopLoss": "StopCodon",
    "ComplexSubstitution": "CDS",
    "Deletion": "CDS",
    "ExonLoss": "CDS",
    "FrameShiftTruncation": "CDS",
    "FrameShift": "CDS",
    "Insertion": "CDS",
    "PrematureStop": "CDS",
    "Substitution": "CDS",
    "Silent": "CDS",
    "ExonicSpliceSite": "SpliceSite",
    "IntronicSpliceSite": "SpliceSite",
    "SpliceAcceptor": "SpliceSite",
    "SpliceDonor": "SpliceSite",
    # "FivePrimeUTR": "5'UTR",
    # "ThreePrimeUTR": "3'UTR",
    # *IncompleteTranscript* | Can't determine effect since transcript annotation is incomplete (often missing either the start or stop codon).
    # *NoncodingTranscript* | Transcript doesn't code for a protein.
    # *Intergenic* | Occurs outside of any annotated gene.
    # *Intragenic* |Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA.
    # *Intronic* | Variant occurs between exons and is unlikely to affect splicing.
}


@dataclass
class Site:
    chrom: str = "."
    pos: int = -1
    strand: str = "."
    ref: str = "-"
    alt: str = "N"

    def to_list(self):
        return [self.chrom, self.pos, self.strand, self.ref, self.alt]


@dataclass
class Annot:
    mut_type: str | None = None
    gene_type: str | None = None
    gene_name: str | None = None
    gene_pos: int | None = None
    transcript_name: str | None = None
    transcript_pos: int | None = None
    transcript_motif: str | None = None
    coding_pos: int | None = None
    codon_ref: str | None = None
    aa_pos: int | None = None
    aa_ref: str | None = None
    distance2splice: int | None = None

    # join into string with tab
    def __str__(self):
        attrs = vars(self)
        return "\t".join([str(x) for x in attrs.values()])

    def get_values(self, as_string=False):
        attrs = vars(self)
        if as_string:
            return list(map(str, attrs.values()))
        return list(attrs.values())

    def get_names(self):
        attrs = vars(self)
        return list(attrs.keys())

    def rename_effect(self, rename_or_not=True):
        if rename_or_not:
            self.mut_type = FEATURE_MAPPER.get(
                str(self.mut_type), self.mut_type
            )
        return self


def expand_base(base):
    return IUPAC[base]


def reverse_base(base):
    return "".join([COMPLEMENT[b] for b in base][::-1])


def mut2eff(chrom, pos, strand, ref, alt, genome, strandness):
    chrom = str(chrom).replace("chr", "")
    pos = int(pos)
    ref = ref.upper()
    alt = alt.upper()

    # pyensembl will turn chrom into all upper case, this nonsenses
    if chrom.upper() not in genome.contigs():
        return varcode.EffectCollection([])

    eff_list = [
        e
        for b in expand_base(alt)
        if b != ref
        for e in varcode.Variant(
            contig=chrom, start=pos, ref=ref, alt=b, ensembl=genome
        ).effects()
    ]

    # filter the correct strand
    if strandness and strand in "+-":
        eff_list_filtered = []
        for e in eff_list:
            if (
                type(e).__name__ == "Intergenic"
                or e.transcript is None
                or e.transcript.strand == strand
            ):
                eff_list_filtered.append(e)
    else:
        eff_list_filtered = eff_list

    # set as intergenic if no transcript is found
    if len(eff_list_filtered) == 0:
        e = Intergenic(
            varcode.Variant(
                contig=chrom,
                start=pos,
                ref=ref,
                alt=[b for b in expand_base(alt)][0],
                ensembl=genome,
            )
        )
        d2g = float("inf")
        expand_intergenic_window = 10000
        for g in genome.genes_at_locus(
            chrom,
            pos - expand_intergenic_window,
            end=pos + expand_intergenic_window,
            strand=strand,
        ):
            d = min(g.start - pos, pos - g.end)
            if d < d2g:
                d2g = d
                e.gene = g
                if hasattr(g, "transcript"):
                    e.transcript = g.transcripts[0]
        eff_list_filtered = [e]

    effs = varcode.EffectCollection(eff_list_filtered)
    return effs


def parse_eff(eff, pos, pad):
    if eff is None:
        return Annot(mut_type="NotInReferenceGenome")
    pos = int(pos)
    mut_type = type(eff).__name__
    # calculate the distance to the splice site
    # TODO: check exon number
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

    try:
        gene_type = eff.gene.biotype
    except:
        gene_type = None
    # put gene name and gene id together
    gene_name = getattr(eff, "gene_id", None)
    if gene_name is not None and hasattr(eff, "gene_name") and eff.gene_name:
        gene_name = gene_name + "(" + eff.gene_name + ")"

    # put transcript name and transcript id together
    transcript_name = getattr(eff, "transcript_id", None)
    if (
        transcript_name is not None
        and hasattr(eff, "transcript_name")
        and eff.transcript_name
    ):
        transcript_name = transcript_name + "(" + eff.transcript_name + ")"

    # NOTE: report distaance to gene and transcript
    # But how to store the nubmer? Negetive for upstream? But what for downstream?
    if mut_type == "Intergenic":
        return Annot(
            mut_type=mut_type,
            gene_type=gene_type,
            gene_name=gene_name,
            transcript_name=transcript_name,
        )
    try:
        gene_pos = eff.gene.offset(pos) + 1
    except:
        gene_pos = None

    non_exon_recored = Annot(
        mut_type=mut_type,
        gene_type=gene_type,
        gene_name=gene_name,
        gene_pos=gene_pos,
        transcript_name=transcript_name,
    )
    if mut_type in ["Intronic", "SpliceDonor"]:
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
        if coding_pos == None:
            aa_pos = None
            aa_ref = None
        else:
            if mut_type == "Silent":
                aa_pos = eff.aa_pos
            elif mut_type in ["IntronicSpliceSite", "ExonicSpliceSite"]:
                aa_pos = eff.alternate_effect.aa_pos
            else:
                aa_pos = eff.aa_mutation_start_offset + 1
            if aa_pos < len(eff.original_protein_sequence):
                aa_ref = eff.original_protein_sequence[(coding_pos - 1) // 3]
            #  "StopLoss"
            elif aa_pos == len(eff.original_protein_sequence):
                aa_ref = "*"
            else:
                aa_ref = None
    else:
        coding_pos = None
        codon_ref = None
        aa_pos = None
        aa_ref = None

    return Annot(
        mut_type,
        gene_type,
        gene_name,
        gene_pos,
        transcript_name,
        transcript_pos,
        transcript_motif,
        coding_pos,
        codon_ref,
        aa_pos,
        aa_ref,
        distance2splice,
    )


def site2mut(
    site,
    genome,
    pad=10,
    strandness=True,
    all_effects=False,
    rename_effect=False,
    pU_mode=False,
):
    chrom, pos, strand, ref, alt = site.to_list()
    effs = mut2eff(chrom, pos, strand, ref, alt, genome, strandness)

    if len(effs) == 0:
        return [parse_eff(None, pos, pad).rename_effect(rename_effect)]
    if not all_effects:
        eff = effect_ordering.get_top_effect(effs, pU_mode=pU_mode)
        return [parse_eff(eff, pos, pad).rename_effect(rename_effect)]
    return [
        parse_eff(eff, pos, pad).rename_effect(rename_effect) for eff in effs
    ]


@click.command(
    help="Variant (genomic variant analysis in python)",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--input",
    "-i",
    "input",
    default="-",
    help="Input position file.",
    required=False,
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
    "--reference-gtf",
    "reference_gtf",
    help="Customized reference gtf file.",
    required=False,
)
@click.option(
    "--reference-transcript",
    "reference_transcript",
    help="Customized reference transcript fasta file.",
    multiple=True,
    required=False,
)
@click.option(
    "--reference-protein",
    "reference_protein",
    help="Customized reference protein fasta file.",
    multiple=True,
    required=False,
)
@click.option(
    "--release",
    "-e",
    "release",
    default=106,
    type=int,
    help="ensembl release",
    required=False,
)
# for backward compartable
@click.option(
    "--type",
    "-t",
    "dna_or_rna",
    type=click.Choice(["DNA", "RNA"], case_sensitive=False),
    default="DNA",
    help="(deprecated)",
)
@click.option(
    "--strandness", "-s", help="Use strand infomation or not?", is_flag=True
)
@click.option(
    "--pU-mode",
    "-u",
    "pU_mode",
    help="Make rRNA, tRNA, snoRNA into top priority.",
    is_flag=True,
)
@click.option(
    "--npad",
    "-n",
    "npad",
    default=10,
    type=int,
    help="Number of padding base to call motif.",
)
@click.option("--all-effects", "-a", help="Output all effects.", is_flag=True)
@click.option(
    "--with-header", "-H", help="With header line in input file.", is_flag=True
)
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3,4,5",
    show_default=True,
    type=str,
    help="Sets columns for site info. (Chrom,Pos,Strand,Ref,Alt)",
)
def run(
    input,
    output,
    reference,
    release,
    reference_gtf,
    reference_transcript,
    reference_protein,
    npad,
    strandness,
    dna_or_rna,
    all_effects,
    pU_mode,
    with_header,
    columns,
):
    columns_index = list(map(lambda x: int(x) - 1, columns.split(",")))
    columns_index_mapper = dict(
        zip(["chrom", "pos", "strand", "ref", "alt"], columns_index)
    )
    # if there is no alt allle, combine effect on CDS
    if len(columns_index) < 5:
        rename_effect = True
    else:
        rename_effect = False
    if dna_or_rna == "RNA":
        strandness = True

    # The max version of GRCm38 is 102
    if (
        reference_gtf
        and len(reference_transcript) > 0
        and len(reference_protein) > 0
    ):
        ensembl_genome = pyensembl.Genome(
            reference_name="customized",
            annotation_name="customized",
            gtf_path_or_url=reference_gtf,
            transcript_fasta_paths_or_urls=reference_transcript,
            protein_fasta_paths_or_urls=reference_protein,
        )

    elif reference in ["mouse", "GRCm39"]:
        ensembl_genome = pyensembl.EnsemblRelease(
            release="108", species="mus_musculus"
        )
    elif reference in ["GRCm38"]:
        ensembl_genome = pyensembl.EnsemblRelease(
            release="102", species="mus_musculus"
        )
    elif reference in ["GRCz11"]:
        ensembl_genome = pyensembl.EnsemblRelease(
            release="108", species="danio_rerio"
        )
    elif reference in ["human", "GRCh38"]:
        ensembl_genome = pyensembl.EnsemblRelease(
            release="108", species="homo_sapiens"
        )
    else:
        if release:
            ensembl_genome = pyensembl.EnsemblRelease(
                release=release, species=reference
            )
        else:
            ensembl_genome = pyensembl.EnsemblRelease(species=reference)
    try:
        ensembl_genome.index()
    except:
        ensembl_genome.download()
        ensembl_genome.index()

    with click.open_file(input, "r") as input_file, click.open_file(
        output, "w"
    ) as output_file:
        if with_header:
            input_header = input_file.readline().strip().split()
        else:
            input_header = ["chrom", "pos", "strand", "ref", "alt"]

        output_file.write("\t".join(input_header + Annot().get_names()) + "\n")
        for l in input_file:
            input_cols = l.strip("\n").split("\t")
            site = Site()
            for n, i in columns_index_mapper.items():
                setattr(site, n, input_cols[i])
            site.ref = site.ref.upper()
            site.alt = site.alt.upper()
            if strandness:
                if site.strand == "-":
                    site.ref = reverse_base(site.ref)
                    site.alt = reverse_base(site.alt)
                elif site.strand == "+":
                    pass
                else:
                    raise ValueError("Strand must be + or -")

            annot_list = site2mut(
                site,
                ensembl_genome,
                npad,
                strandness,
                all_effects,
                rename_effect,
                pU_mode,
            )
            for annot in annot_list:
                output_file.write(
                    "\t".join(input_cols + annot.get_values(as_string=True))
                    + "\n"
                )


if __name__ == "__main__":
    run()
