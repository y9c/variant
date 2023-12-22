# Python pakcage for genomic variant analysis

[![Pypi Releases](https://img.shields.io/pypi/v/variant.svg)](https://pypi.python.org/pypi/variant)
[![Downloads](https://pepy.tech/badge/variant)](https://pepy.tech/project/variant)

## How to use?

```
pip install variant
```

## `variant motif` subcommand can fetch motif sequence around given site.

```
 Usage: variant motif [OPTIONS]

 Fetch genomic motif.

╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│    --input        -i  TEXT  Input position file.                                                                                                                                                                                                                                                                │
│    --output       -o  TEXT  Output annotation file.                                                                                                                                                                                                                                                             │
│ *  --fasta        -f  TEXT  reference fasta file. [required]                                                                                                                                                                                                                                                    │
│    --npad         -n  TEXT  Number of padding base to call motif. If you want to set different left and right pads, use comma to separate them. (eg. 2,3)                                                                                                                                                       │
│    --with-header  -H        With header line in input file.                                                                                                                                                                                                                                                     │
│    --columns      -c  TEXT  Sets columns for site info. (Chrom,Pos,Strand) [default: 1,2,3]                                                                                                                                                                                                                     │
│    --to-upper     -u        Convert motif to upper case.                                                                                                                                                                                                                                                        │
│    --wrap-site    -w        Wrap motif site.                                                                                                                                                                                                                                                                    │
│    --help         -h        Show this message and exit.                                                                                                                                                                                                                                                         │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

> demo:

I would like to get the 2 bases before the given sites, and 3 bases after the given sites, meanwhile, wrap the give sites with bracket. Moreover, the strand information should be taken into account. 

use `-n 2,3 -w`

## `variant effect` subcommand can infer the effect of a mutation

```
 Usage: variant effect [OPTIONS]

 Annotation genomic variant effect.

╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --input                 -i  TEXT     Input position file.                                                                                                                                                                                                                                                       │
│ --output                -o  TEXT     Output annotation file                                                                                                                                                                                                                                                     │
│ --reference             -r  TEXT     reference species                                                                                                                                                                                                                                                          │
│ --reference-gtf             TEXT     Customized reference gtf file.                                                                                                                                                                                                                                             │
│ --reference-transcript      TEXT     Customized reference transcript fasta file.                                                                                                                                                                                                                                │
│ --reference-protein         TEXT     Customized reference protein fasta file.                                                                                                                                                                                                                                   │
│ --reference-mapping         TEXT     Mapping file for chrom name, first column is chrom in the input, second column is chrom in the reference db (sep by tab)                                                                                                                                                   │
│ --release               -e  INTEGER  ensembl release                                                                                                                                                                                                                                                            │
│ --strandness            -s           Use strand infomation or not?                                                                                                                                                                                                                                              │
│ --pU-mode               -u           Make rRNA, tRNA, snoRNA into top priority.                                                                                                                                                                                                                                 │
│ --npad                  -n  INTEGER  Number of padding base to call motif.                                                                                                                                                                                                                                      │
│ --all-effects           -a           Output all effects.                                                                                                                                                                                                                                                        │
│ --with-header           -H           With header line in input file.                                                                                                                                                                                                                                            │
│ --columns               -c  TEXT     Sets columns for site info. (Chrom,Pos,Strand,Ref,Alt) [default: 1,2,3,4,5]                                                                                                                                                                                                │
│ --help                  -h           Show this message and exit.                                                                                                                                                                                                                                                │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

> demo:

Store the following table in file (`sites.tsv`).

| Chrom | Position  | Strand | Ref | Alt |
| ----- | --------- | ------ | --- | --- |
| chr1  | 230703034 | -      | C   | T   |
| chr12 | 69353439  | +      | A   | T   |
| chr14 | 23645352  | +      | G   | T   |
| chr2  | 215361150 | -      | A   | T   |
| chr2  | 84906537  | +      | C   | T   |
| chr22 | 39319077  | -      | T   | A   |
| chr22 | 39319095  | -      | T   | A   |
| chr22 | 39319098  | -      | T   | A   |

Run command:

```bash
variant-effect -i sites.tsv -H -r human -e 108 -t RNA -H -c 1,2,3
```

- `-i` specify the input file
- `-H` means the file is with header line, and the first row will be skipped;
- `-r` use the specific genome, default is human
- `-e` specify the Ensembl release version
- `-c` means only use some of the columns in the input file. default will use the first 5 columns.

You will have this output

| Chrom | Position  | Strand | Ref | Alt | mut_type      | gene_type      | gene_name               | gene_pos | transcript_name             | transcript_pos | transcript_motif      | coding_pos | codon_ref | aa_pos | aa_ref | distance2splice |
| :---- | :-------- | :----- | :-- | :-- | :------------ | :------------- | :---------------------- | :------- | :-------------------------- | :------------- | :-------------------- | :--------- | :-------- | :----- | :----- | --------------- |
| chr1  | 230703034 | -      | C   | T   | ThreePrimeUTR | protein_coding | ENSG00000135744(AGT)    | 42543    | ENST00000680041(AGT-208)    | 1753           | TGTGTCACCCCCAGTCTCCCA | None       | None      | None   | None   | 295             |
| chr12 | 69353439  | +      | A   | T   | ThreePrimeUTR | protein_coding | ENSG00000090382(LYZ)    | 5059     | ENST00000261267(LYZ-201)    | 695            | TAGAACTAATACTGGTGAAAA | None       | None      | None   | None   | 286             |
| chr14 | 23645352  | +      | G   | T   | ThreePrimeUTR | protein_coding | ENSG00000100867(DHRS2)  | 15238    | ENST00000344777(DHRS2-202)  | 1391           | CTGCCATTCTGCCAGACTAGC | None       | None      | None   | None   | 210             |
| chr2  | 215361150 | -      | A   | T   | ThreePrimeUTR | protein_coding | ENSG00000115414(FN1)    | 74924    | ENST00000323926(FN1-201)    | 8012           | GGCCCGCAATACTGTAGGAAC | None       | None      | None   | None   | 476             |
| chr2  | 84906537  | +      | C   | T   | ThreePrimeUTR | protein_coding | ENSG00000034510(TMSB10) | 882      | ENST00000233143(TMSB10-201) | 327            | CCTGGGCACTCCGCGCCGATG | None       | None      | None   | None   | 148             |
| chr22 | 39319077  | -      | T   | A   | Intronic      | protein_coding | ENSG00000100316(RPL3)   | 1313     | ENST00000216146(RPL3-201)   | None           | None                  | None       | None      | None   | None   | None            |
| chr22 | 39319095  | -      | T   | A   | Intronic      | protein_coding | ENSG00000100316(RPL3)   | 1295     | ENST00000216146(RPL3-201)   | None           | None                  | None       | None      | None   | None   | None            |
| chr22 | 39319098  | -      | T   | A   | Intronic      | protein_coding | ENSG00000100316(RPL3)   | 1292     | ENST00000216146(RPL3-201)   | None           | None                  | None       | None      | None   | None   | None            |

## ⏳⏳⏳ more functions will be supported in the future

## TODO:

- imporve speed. Base on [cgranges](https://github.com/lh3/cgranges), [pyranges](https://github.com/biocore-ntnu/pyranges)?, or [BioCantor](https://github.com/InscriptaLabs/BioCantor)?
