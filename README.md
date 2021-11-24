# Python pakcage for genomic variant analysis

## `variant-effect` command can infer the effect of a mutation

The input file has 5 columns: `chromosome`, `position`, `strand`, `reference allele`, `alternative allele`.

- No header is required.
- The 3rd column (strand) is not used by default, just for compatibility with RNA mode.
- By default, the base of reference and alternative allele are based on DNA information
- For RNA mode (through `--rna` argument), the base of reference and alternative allele is reverse complement if the strand is negative(-).

eg:

```
chr16   400560      .      G       T
chr17   41690930    .      G       T
chr6    61574496    .      A       T
chr2    84906522    .      G       T
chr2    216205243   .      G       T
chr4    73455665    .      G       T
chr2    101891316   .      G       T
chr2    69820761    .      G       T
chr6    30723661    .      A       T
```

- The output can be stdout stdout, or a file.

```
#chrom  pos        strand  ref     alt     mut_type        gene_name       transcript_id   transcript_pos  transcript_motif        coding_pos      codon_ref       aa_pos  aa_ref
chr16   400560     .       G       T       ThreePrimeUTR   NME4    ENST00000219479 806     GCACCAAAGTGCCGGACAACC   None    None    None    None
chr17   41690930   .       G       T       Substitution    EIF1    ENST00000591776 515     CTTGTATAATGTAACCATTTG   363     ATG     121     M
chr6    61574496   .       A       T       Intergenic      None    None    None    None    None    None    None    None                                                                                                                        chr2    84906522        G       T       ThreePrimeUTR   TMSB10  ENST00000233143 312     AAGCTGCACTGTGAACCTGGG   None    None    None    None
chr2    216205243  .       G       T       ThreePrimeUTR   XRCC5   ENST00000392133 2701    TGCCATCGCTGTGATGCTGGG   None    None    None    None
chr4    73455665   .       G       T       Substitution    AFP     ENST00000226359 1836    TTCATTCGGTGTGAACTTTTC   1820    TGT     607     C
chr2    101891316  .       G       T       ThreePrimeUTR   MAP4K4  ENST00000350878 4267    GGAATTCCTTGTAACTGGAGC   None    None    None    None
chr2    69820761   .       G       T       Substitution    ANXA4   ENST00000394295 934     AAATTGACATGTTGGATATCC   846     ATG     282     M
chr6    30723661   .       A       T       Substitution    TUBB    ENST00000327892 754     GATGAGACCTATTGCATTGAC   599     TAT     200     Y
```
