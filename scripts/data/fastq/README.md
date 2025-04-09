# Simulated Samples

## DNA Samples 

Ground Truth

| Sample ID     | Variant ID | Read Type           | Chromosome 1 | Position 1 | Orientation 1 | Chromosome 2 | Position 2  | Orientation 2 | Variant Type   | Variant Sequence               | Description                  |
|---------------|------------|---------------------|--------------|------------|---------------|--------------|-------------|---------------|----------------|--------------------------------|------------------------------|
| dna-001-tumor | 1          | long-read           | chr17        | 7674224    | D             | chr17        | 7674226     | U             | SNV            | A                              | TP53                         |
| dna-002-tumor | 2          | long-read           | chr17        | 7674225    | D             | chr17        | 7674226     | U             | INS            | ACGTACGTGGTATGCATGCTGAGACTGAGG | TP53                         |
| dna-003-tumor | 3          | long-read           | chr17        | 7674200    | D             | chr17        | 7674231     | U             | DEL            |                                | TP53                         |
| dna-004-tumor | 4          | long-read           | chr17        | 7670399    | D             | chr17        | 7680500     | D             | BND            |                                | Inversion (TP53)             |
| dna-004-tumor | 5          | long-read           | chr17        | 7670400    | U             | chr17        | 7680501     | U             | BND            |                                | Inversion (TP53)             |
| dna-005-tumor | 6          | long-read           | chr17        | 4637154    | D             | chr17        | 7674880     | U             | BND            | TATATACGAGCGTACGTGACTGGTACGTTA | Translocation (ALOX15-TP53)  |
| dna-006-tumor | 7          | long-read           | chr17        | 7676155    | D             | chr18        | 5170100     | U             | BND            |                                | Translocation (TP53-AKAIN1)  |
| dna-007-tumor | 8          | long-read           | chr17        | 7676155    | D             | chr18        | 5170100     | D             | BND            |                                | Translocation (TP53-AKAIN1)  |

## RNA Samples

Ground Truth

| Sample ID     | Variant ID | Read Type           | Chromosome 1 | Position 1 | Orientation 1 | Chromosome 2 | Position 2 | Orientation 2 | Variant Type | Variant Sequence | Description                            |
|---------------|------------|---------------------|--------------|------------|---------------|--------------|------------|---------------|--------------|------------------|----------------------------------------|
| rna-100-tumor | 100        | long-read           | chr17        | 7674224    | D             | chr17        | 7674226    | U             | SNV          | A                | TP53                                   |
| rna-101-tumor | 101        | long-read           | chr17        | 7676400    | D             | chr17        | 7676401    | U             | INS          | GGGGGTTTTT       | TP53                                   |
| rna-102-tumor | 102        | long-read           | chr17        | 7673750    | D             | chr17        | 7673761    | U             | DEL          |                  | TP53                                   |
| rna-103-tumor | 103        | long-read           | chr17        | 7701730    | D             | chr17        | 7727200    | U             | FUS          |                  | Fusion gene (WRAP53-DNAH2)             |
| rna-104-tumor | 104        | long-read           | chr17        | 7676099    | *             | chr17        | 7676099    | *             | A5P          |                  | Alternative 5 prime splice site (TP53) |
| rna-105-tumor | 105        | long-read           | chr17        | 7676201    | *             | chr17        | 7676201    | *             | A3P          |                  | Alternative 3 prime splice site (TP53) |
| rna-106-tumor | 106        | long-read           | chr17        | 7675994    | *             | chr17        | 7676270    | *             | SKP          |                  | Exon skipping (TP53)                   |
| rna-107-tumor | 107        | long-read           | chr17        | 7675601    | *             | chr17        | 7675641    | *             | CRX          |                  | Cryptic exon (TP53)                    |
| rna-108-tumor | 108        | long-read           | chr17        | 7675976    | *             | chr17        | 7675993    | *             | IRT          |                  | Intron retention (TP53)                |
