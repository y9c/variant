# Python pakcage for genomic variant analysis

[![Pypi Releases](https://img.shields.io/pypi/v/variant.svg)](https://pypi.python.org/pypi/variant)
[![Downloads](https://pepy.tech/badge/variant)](https://pepy.tech/project/variant)

## `variant-effect` command can infer the effect of a mutation

- `-i/--input` to sepecify the input file. The input file has 5 columns: `chromosome`, `position`, `strand`, `reference allele`, `alternative allele`.

  - No header is required.
  - The 3rd column (strand) is not used by default, just for compatibility with RNA mode.
  - By default, the base of reference and alternative allele are based on DNA information
  - For RNA mode (through `--rna` argument), the base of reference and alternative allele is reverse complement if the strand is negative(-).

- `-o/--output` to specify the output file, leave empty for stdout.
- `-r/--reference` to specify reference name, can be human / mouse / dog / cat / chicken ...
- `-t/--type [DNA|RNA]` to run in DNA or RNA mode. If RNA is specified, the ref base will be complemented.
- `--all-effects` output all effects of the variant.

> demo:

Store the following table in sites.tsv.

```
chr3    10301112        -       G       T
chr7    94669540        +       G       N
chr2    215361150       -       A       T
chr15   72199549        -       G       T
chr17   81843580        -       C       T
chr2    84906537        +       C       T
chr14   23645352        +       G       T
chr20   37241351        +       G       T
chrX    153651037       +       G       T
chr17   81844010        -       A       T
```

Run command `variant-effect -i sites.tsv -r human -t RNA` to get the following output.

```
#chrom  pos     strand  ref     alt     mut_type        gene_name       transcript_id   transcript_pos  transcript_motif        coding_pos      codon_ref       aa_pos  aa_ref
chr3    10301112        -       C       A       Silent  SEC13   ENST00000397117 1441    TTGATCATCTGCCTTAACGTG   849     CTG     284     L
chr7    94669540        +       G       N       ThreePrimeUTR   PEG10   ENST00000612941 6240    TTTTACCCCTGTCAGTAGCCC   None    None    None    None
chr2    215361150       -       T       A       ThreePrimeUTR   FN1     ENST00000323926 8012    GGCCCGCAATACTGTAGGAAC   None    None    None    None
chr15   72199549        -       C       A       ThreePrimeUTR   PKM     ENST00000319622 2197    GCTGTAACGTGGCACTGGTAG   None    None    None    None
chr17   81843580        -       G       A       ThreePrimeUTR   P4HB    ENST00000681020 3061    AGAAGCTTGTCCCCCGTGTGG   None    None    None    None
chr2    84906537        +       C       T       ThreePrimeUTR   TMSB10  ENST00000233143 327     CCTGGGCACTCCGCGCCGATG   None    None    None    None
chr14   23645352        +       G       T       ThreePrimeUTR   DHRS2   ENST00000344777 1391    CTGCCATTCTGCCAGACTAGC   None    None    None    None
chr20   37241351        +       G       T       ThreePrimeUTR   RPN2    ENST00000237530 1959    AAAACTGAATGTCAAGAAAAG   None    None    None    None
chrX    153651037       +       G       T       ThreePrimeUTR   DUSP9   ENST00000342782 2145    CTGCTACTTTGGGGGGTGGGG   None    None    None    None
chr17   81844010        -       T       A       ThreePrimeUTR   P4HB    ENST00000681020 2631    GAACTGTAATACGCAAAGCCA   None    None    None    None
```

TODO:

- support GRCh37
