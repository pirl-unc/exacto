# Simulated Samples

## DNA Samples 

Ground Truth

| Sample ID     | Variant ID | Read Type           | Chromosome 1 | Position 1 | Operation 1 | Chromosome 2 | Position 2  | Operation 2 | Variant Type | Sequence                        | Description                                               |
|---------------|------------|---------------------|--------------|------------|-------------|--------------|-------------|-------------|--------------|---------------------------------|-----------------------------------------------------------|
| dna-001-tumor | 1          | long-read           | chr17        | 7674224    | D           | chr17        | 7674226     | U           | SNV          | A                               | TP53 SNV                                                  |
| dna-002-tumor | 2          | long-read           | chr17        | 7674225    | D           | chr17        | 7674226     | U           | INS          | ACGTACGTGGTATGCATGCTGAGACTGAGG  | TP53 insertion                                            |
| dna-003-tumor | 3          | long-read           | chr17        | 7674200    | D           | chr17        | 7674231     | U           | DEL          |                                 | TP53 deletion                                             |
| dna-004-tumor | 4          | long-read           | chr17        | 7670399    | D           | chr17        | 7680500     | D           | BND          |                                 | TP53 inversion                                            |
| dna-004-tumor | 5          | long-read           | chr17        | 7670400    | U           | chr17        | 7680501     | U           | BND          |                                 | TP53 inversion                                            |
| dna-005-tumor | 6          | long-read           | chr17        | 4637154    | D           | chr17        | 7674880     | U           | BND          | TATATACGAGCGTACGTGACTGGTACGTTA  | Translocation (ALOX15-TP53)                               |
| dna-006-tumor | 7          | long-read           | chr17        | 7676155    | D           | chr18        | 5170100     | U           | TRA          |                                 | Translocation (TP53-AKAIN1) - interchromosomal breakpoint |
| dna-007-tumor | 8          | long-read           | chr17        | 7676155    | D           | chr18        | 5170100     | D           | TRA          |                                 | Translocation (TP53-AKAIN1) - interchromosomal breakpoint |
| dna-008-tumor | 9          | long-read           | chr17        | 7668421    | U           | chr17        | 7687490     | D           | DUP          |                                 | TP53 duplication                                          |
| dna-009-tumor | 10         | long-read           | chr17        | 7687490    | D           | chr17        | 7687490     | D           | INVDUP       |                                 | TP53 inverted duplication                                 |
| dna-009-tumor | 11         | long-read           | chr17        | 7668421    | U           | chr17        | 7687490     | U           | INVDUP       |                                 | TP53 inverted duplication                                 |
| dna-010-tumor | 12         | long-read           | chr17        | 7687490    | D           | chr17        | 7668421     | U           | INVDUP       |                                 | TP53 inverted duplication                                 |
| dna-010-tumor | 13         | long-read           | chr17        | 7687490    | D           | chr17        | 7687490     | D           | INVDUP       |                                 | TP53 inverted duplication                                 |
| dna-010-tumor | 14         | long-read           | chr17        | 7668421    | U           | chr17        | 7687490     | U           | INVDUP       |                                 | TP53 inverted duplication                                 |
| dna-011-tumor | 15         | long-read           | chr17        | 7687490    | U           | chr17        | 7687491     | U           | INS          | TATCTCGCGAATTCAGCTACTACTACGGGA  | TP53 insertion (softclipped)                              |
| dna-011-tumor | 16         | long-read           | chr17        | 7668420    | U           | chr17        | 7668421     | U           | INS          | AGCGGCGAATATCAGCTACCTCTTAAGATC  | TP53 insertion (softclipped)                              |

## RNA Samples

Ground Truth

| Sample ID     | Variant ID | Read Type           | Chromosome 1 | Position 1 | Operation 1 | Chromosome 2 | Position 2 | Operation 2 | Variant Type | Variant Sequence | Description                                              |
|---------------|------------|---------------------|--------------|------------|---------------|--------------|------------|---------------|--------------|------------------|----------------------------------------------------------|
| rna-100-tumor | 100        | long-read           | chr17        | 7674224    | D             | chr17        | 7674226    | U             | SNV          | A                | TP53 (exonic single-nucleotide variant)                  |
| rna-101-tumor | 101        | long-read           | chr17        | 7676400    | D             | chr17        | 7676401    | U             | INS          | GGGGGTTTTT       | TP53 (exonic insertion)                                  |
| rna-102-tumor | 102        | long-read           | chr17        | 7673750    | D             | chr17        | 7673761    | U             | DEL          |                  | TP53 (exonic deletion                                    |
| rna-103-tumor | 103        | long-read           | chr17        | 7701730    | D             | chr17        | 7727200    | U             | FUS          |                  | Fusion gene (WRAP53-DNAH2) - adjacent genes              |
| rna-104-tumor | 104        | long-read           | chr17        | 7676099    | *             | chr17        | 7676099    | *             | A5P          |                  | Alternative 5 prime splice site (TP53) - exon truncation |
| rna-105-tumor | 105        | long-read           | chr17        | 7676201    | *             | chr17        | 7676201    | *             | A3P          |                  | Alternative 3 prime splice site (TP53) - exon truncation |
| rna-106-tumor | 106        | long-read           | chr17        | 7675994    | *             | chr17        | 7676270    | *             | SKP          |                  | Exon skipping (TP53) - exon 4                            |
| rna-107-tumor | 107        | long-read           | chr17        | 7675601    | *             | chr17        | 7675641    | *             | CRX          |                  | Cryptic exon (TP53) - between exons 4 and 5              |
| rna-108-tumor | 108        | long-read           | chr17        | 7675976    | *             | chr17        | 7675993    | *             | IRT          |                  | Intron retention (TP53)                                  |
| rna-109-tumor | 1090       | long-read           | chr17        | 1295600    | D             | chr17        | 3801188    | U             | FUS          |                  | Fusion gene (TRARG1-ITGAE)                               |
| rna-109-tumor | 1091       | long-read           | chr17        | 3761100    | U             | chr17        | 7727200    | U             | FUS          |                  | Fusion gene (ITGAE-DNAH2)                                |
| rna-110-tumor | 1100       | long-read           | chr17        | 2464300    | U             | chr17        | 4433940    | U             | FUS          |                  | Fusion gene (METTL16-SPNS3)                              |
| rna-110-tumor | 1101       | long-read           | chr17        | 4453100    | D             | chr17        | 7603800    | U             | FUS          |                  | Fusion gene (SPNS3-FXR2)                                 |
| rna-111-tumor | 111        | long-read           | chr17        | 1844686    | D             | chr17        | 1842803    | U             | CIR          |                  | Circular RNA (RPA1 exons 2-4)                            |
