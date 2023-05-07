import rich_click as click


@click.group(
    invoke_without_command=False,
    help="Variant (genomic variant analysis in python)",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option(None, "-v", "--version")
@click.pass_context
def cli(ctx):
    pass


@cli.command(
    help="Annotation genomic variant effect.",
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
    default=109,
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
def effect(
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
    from .effect import run_effect

    run_effect(
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
    )


@cli.command(
    help="Fetch genomic motif.",
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
    help="Output annotation file.",
    required=False,
)
@click.option(
    "--fasta",
    "-f",
    "fasta",
    help="reference fasta file.",
    required=True,
)
@click.option(
    "--npad",
    "-n",
    "npad",
    default=10,
    type=int,
    help="Number of padding base to call motif.",
)
@click.option(
    "--with-header", "-H", help="With header line in input file.", is_flag=True
)
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3",
    show_default=True,
    type=str,
    help="Sets columns for site info. (Chrom,Pos,Strand)",
)
def motif(input, output, fasta, npad, with_header, columns):
    from .motif import run_motif

    run_motif(input, output, fasta, npad, with_header, columns)


if __name__ == "__main__":
    cli()
