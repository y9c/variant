# Python pakcage for genomic variant analysis

[![Pypi Releases](https://img.shields.io/pypi/v/variant.svg)](https://pypi.python.org/pypi/variant)
[![Downloads](https://pepy.tech/badge/variant)](https://pepy.tech/project/variant)

## How to use?

```
pip install variant
```

- run `variant-effect` in the command line
- more functions will be supported in the future

## `variant-effect` command can infer the effect of a mutation

```
 Usage: variant-effect [OPTIONS]

 Variant (genomic variant analysis in python)

╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --input                 -i  TEXT       Input position file.                                                          │
│ --output                -o  TEXT       Output annotation file                                                        │
│ --reference             -r  TEXT       reference species                                                             │
│ --reference-gtf             TEXT       Customized reference gtf file.                                                │
│ --reference-transcript      TEXT       Customized reference transcript fasta file.                                   │
│ --reference-protein         TEXT       Customized reference protein fasta file.                                      │
│ --release               -e  INTEGER    ensembl release                                                               │
│ --type                  -t  [DNA|RNA]  (deprecated)                                                                  │
│ --strandness            -s             Use strand infomation or not?                                                 │
│ --pU-mode               -u             Make rRNA, tRNA, snoRNA into top priority.                                    │
│ --npad                  -n  INTEGER    Number of padding base to call motif.                                         │
│ --all-effects           -a             Output all effects.                                                           │
│ --with-header           -H             With header line in input file.                                               │
│ --columns               -c  TEXT       Sets columns for site info. (Chrom,Pos,Strand,Ref,Alt) [default: 1,2,3,4,5]   │
│ --help                  -h             Show this message and exit.                                                   │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

> demo:

Store the following table in file (`sites.tsv`).

|       |           |     |     |     |
| :---- | :-------- | :-- | :-- | :-- |
| chr1  | 230703034 | -   | C   | T   |
| chr2  | 215361150 | -   | A   | T   |
| chr2  | 84906537  | +   | C   | T   |
| chr3  | 10301112  | -   | G   | T   |
| chr3  | 20301112  | -   | G   | T   |
| chr7  | 45893389  | +   | G   | T   |
| chr7  | 94669540  | +   | G   | N   |
| chr12 | 69353439  | +   | A   | T   |
| chr14 | 23645352  | +   | G   | T   |

Run command `variant-effect -i sites.tsv -r human -e 106 -t RNA` to get the following output.

| #chrom | pos       | strand | ref | alt | mut_type      | gene_name               | gene_pos | transcript_name             | transcript_pos | transcript_motif      | coding_pos | codon_ref | aa_pos | aa_ref | distance2splice |
| :----- | :-------- | :----- | :-- | :-- | :------------ | :---------------------- | :------- | :-------------------------- | :------------- | :-------------------- | :--------- | :-------- | :----- | :----- | :-------------- |
| chr1   | 230703034 | -      | G   | A   | ThreePrimeUTR | ENSG00000135744(AGT)    | 42543    | ENST00000680041(AGT-208)    | 1753           | TGTGTCACCCCCAGTCTCCCA | None       | None      | None   | None   | 295             |
| chr2   | 215361150 | -      | T   | A   | ThreePrimeUTR | ENSG00000115414(FN1)    | 74924    | ENST00000323926(FN1-201)    | 8012           | GGCCCGCAATACTGTAGGAAC | None       | None      | None   | None   | 476             |
| chr2   | 84906537  | +      | C   | T   | ThreePrimeUTR | ENSG00000034510(TMSB10) | 882      | ENST00000233143(TMSB10-201) | 327            | CCTGGGCACTCCGCGCCGATG | None       | None      | None   | None   | 148             |
| chr3   | 10301112  | -      | C   | A   | Silent        | ENSG00000157020(SEC13)  | 20001    | ENST00000397117(SEC13-209)  | 1441           | TTGATCATCTGCCTTAACGTG | 849        | CTG       | 283    | L      | 35              |
| chr3   | 20301112  | -      | C   | A   | Intergenic    | None                    | None     | None                        | None           | None                  | None       | None      | None   | None   | None            |
| chr7   | 45893389  | +      | G   | T   | ThreePrimeUTR | ENSG00000146678(IGFBP1) | 5030     | ENST00000275525(IGFBP1-201) | 1243           | CAAAGCTCCTGCGTCTGTTTT | None       | None      | None   | None   | 429             |
| chr7   | 94669540  | +      | G   | N   | ThreePrimeUTR | ENSG00000242265(PEG10)  | 13216    | ENST00000612941(PEG10-206)  | 6240           | TTTTACCCCTGTCAGTAGCCC | None       | None      | None   | None   | 5030            |
| chr12  | 69353439  | +      | A   | T   | ThreePrimeUTR | ENSG00000090382(LYZ)    | 5059     | ENST00000261267(LYZ-201)    | 695            | TAGAACTAATACTGGTGAAAA | None       | None      | None   | None   | 286             |
| chr14  | 23645352  | +      | G   | T   | ThreePrimeUTR | ENSG00000100867(DHRS2)  | 15238    | ENST00000344777(DHRS2-202)  | 1391           | CTGCCATTCTGCCAGACTAGC | None       | None      | None   | None   | 210             |

## TODO:

- imporve speed. Base on [cgranges](https://github.com/lh3/cgranges), [pyranges](https://github.com/biocore-ntnu/pyranges)?, or [BioCantor](https://github.com/InscriptaLabs/BioCantor)?
