# Simulated Samples

## Ground Truth

The scripts in `01_fasta/` simulate one case/control pair per sample and write each sample's truth to
`test/data/simulation/ground_truth/<sample>_ground_truth.tsv`; the tables below list those files.
A strand of `*` means the truth row is unstranded. RNA sequences are written on the transcript strand,
so TP53 (reverse strand) variants read as the reverse complement of the chr17 + strand change.

### DNA Samples

| Sample ID | Chromosome 1 | Position 1 | Strand 1 | Operation 1 | Chromosome 2 | Position 2 | Strand 2 | Operation 2 | Variant Type | Sequence | Description |
|:---|---|---|---|---|---|---|---|---|---|:---|:---|
| scga-mini-dna-001-tumor | chr17 | 7674224 | * | D | chr17 | 7674226 | * | U | SNV | `A` | TP53 SNV (C>A) |
| scga-mini-dna-002-tumor | chr17 | 7674224 | * | D | chr17 | 7674225 | * | U | INS | `CAGGCGGATGGG` | TP53 INS (12 bp) |
| scga-mini-dna-003-tumor | chr17 | 7674200 | * | D | chr17 | 7674231 | * | U | DEL |  | TP53 DEL (30 bp) |
| scga-mini-dna-004-tumor | chr17 | 7674224 | * | D | chr17 | 7674225 | * | U | INS | `CAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGGCAGGCGGATGGG` | TP53 INS (120 bp) |
| scga-mini-dna-005-tumor | chr17 | 7673800 | * | D | chr17 | 7674901 | * | U | DEL |  | TP53 DEL (1,100 bp) |
| scga-mini-dna-006-tumor | chr17 | 7701200 | * | D | chr17 | 7717000 | * | D | BND |  | WRAP53 INV |
| scga-mini-dna-006-tumor | chr17 | 7701201 | * | U | chr17 | 7717001 | * | U | BND |  | WRAP53 INV |
| scga-mini-dna-007-tumor | chr17 | 3491600 | * | D | chr17 | 6085001 | * | U | BND |  | ASPA-WSCD1 fusion (chr17:3491601-6085000 deleted) |
| scga-mini-dna-008-tumor | chr17 | 7672799 | * | D | chr17 | 7672801 | * | U | SNV | `C` | TP53 SNV (T>C) |
| scga-mini-dna-009-tumor | chr17 | 7673533 | * | D | chr17 | 7673535 | * | U | SNV | `A` | TP53 SNV (C>A) |
| scga-mini-dna-010-tumor | chr17 | 7674222 | * | D | chr17 | 7674225 | * | U | MNV | `CC` | TP53 MNV (TT>CC) |
| scga-mini-dna-011-tumor | chr17 | 7673533 | * | D | chr17 | 7673535 | * | U | SNV | `A` | TP53 SNV (C>A; same variant as dna-009) |
| scga-mini-dna-012-tumor | chr17 | 7673300 | * | U | chr17 | 7679900 | * | D | DUP |  | TP53 DUP (chr17:7673300-7679900) |
| scga-mini-dna-013-tumor | chr17 | 7679900 | * | D | chr17 | 7679900 | * | D | BND |  | TP53 inverted DUP (chr17:7673301-7679900) |
| scga-mini-dna-013-tumor | chr17 | 7673301 | * | U | chr17 | 7679901 | * | U | BND |  | TP53 inverted DUP (chr17:7673301-7679900) |
| scga-mini-dna-014-tumor | chr17 | 7669650 | * | D | chr17 | 7669651 | * | U | INS | `TTAACTGAGTCTCAAAAAAATAAA` | TP53 INS (24 bp) |
| scga-mini-dna-015-tumor | chr17 | 3491600 | * | D | chr17 | 6085001 | * | U | BND |  | ASPA-WSCD1-ACAP1 fusion (two junctions) |
| scga-mini-dna-015-tumor | chr17 | 6115000 | * | D | chr17 | 7340400 | * | U | BND |  | ASPA-WSCD1-ACAP1 fusion (two junctions) |
| scga-mini-dna-016-tumor | chr17 | 3491600 | * | D | chr17 | 6085001 | * | U | BND | `CCCATCCGCCTG` | ASPA-WSCD1 fusion with a 12 bp insertion at the junction |

### RNA Samples

| Sample ID | Chromosome 1 | Position 1 | Strand 1 | Operation 1 | Chromosome 2 | Position 2 | Strand 2 | Operation 2 | Variant Type | Sequence | Description |
|:---|---|---|---|---|---|---|---|---|---|:---|:---|
| scga-mini-rna-001-tumor | chr17 | 7674224 | - | D | chr17 | 7674226 | - | U | SNV | `T` | TP53 SNV (ATG>ATT; p.M246I) |
| scga-mini-rna-002-tumor | chr17 | 7674224 | - | D | chr17 | 7674225 | - | U | INS | `CCCATCCGCCTG` | TP53 INS (12 bp; in-frame PIRL) |
| scga-mini-rna-003-tumor | chr17 | 7674200 | - | D | chr17 | 7674231 | - | U | DEL |  | TP53 DEL (30 bp; GMNRRPILTI deleted) |
| scga-mini-rna-004-tumor | chr17 | 7674224 | - | D | chr17 | 7674225 | - | U | INS | `CCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTGCCCATCCGCCTG` | TP53 INS (120 bp; in-frame PIRL x 10) |
| scga-mini-rna-005-tumor | chr17 | 7673800 | - | D | chr17 | 7674901 | - | U | NCS |  | TP53 exon 7 skipping with partial deletions of exons 6 and 8 |
| scga-mini-rna-006-tumor | chr17 | 7700829 | + | D | chr17 | 7716251 | - | D | BND |  | WRAP53 INV |
| scga-mini-rna-007-tumor | chr17 | 3489342 | + | D | chr17 | 6087990 | + | U | FUS |  | ASPA-WSCD1 FUS |
| scga-mini-rna-008-tumor | chr17 | 7672205 | - | I | chr17 | 7672484 | - | I | CRX |  | TP53 cryptic exon (280 bp) |
| scga-mini-rna-009-tumor | chr17 | 7673527 | - | I | chr17 | 7673534 | - | I | IRT |  | TP53 intron retention (8 bp) |
| scga-mini-rna-010-tumor | chr17 | 7674222 | - | D | chr17 | 7674225 | - | U | MNV | `GG` | TP53 MNV (TT>CC on chr17 +) |
| scga-mini-rna-011-tumor | chr17 | 7670609 | - | U | chr17 | 7674290 | - | D | CIR |  | TP53 circular RNA (exons 7-10) |
| scga-mini-rna-012-tumor | chr17 | 7673701 | * | U | chr17 | 7675236 | * | D | DUP |  | TP53 DUP (exons 5-8) |
| scga-mini-rna-013-tumor | chr17 | 7673701 | - | U | chr17 | 7673701 | + | U | BND |  | TP53 inverted DUP (exons 5-8) |
| scga-mini-rna-013-tumor | chr17 | 7673608 | - | D | chr17 | 7675236 | + | D | BND |  | TP53 inverted DUP (exons 5-8) |
| scga-mini-rna-014-tumor | chr17 | 7669650 | - | D | chr17 | 7669651 | - | U | INS | `TTTATTTTTTTGAGACTCAGTTAA` | TP53 INS (24 bp in exon 11; premature stop) |
| scga-mini-rna-015-tumor | chr17 | 3489342 | + | D | chr17 | 6087990 | + | U | FUS |  | ASPA-WSCD1-ACAP1 FUS (two junctions) |
| scga-mini-rna-015-tumor | chr17 | 6110935 | + | D | chr17 | 7341948 | + | U | FUS |  | ASPA-WSCD1-ACAP1 FUS (two junctions) |
| scga-mini-rna-016-tumor | chr17 | 3489342 | + | D | chr17 | 6087990 | + | U | FUS | `CCCATCCGCCTG` | ASPA-WSCD1 FUS with a 12 bp insertion at the junction |
