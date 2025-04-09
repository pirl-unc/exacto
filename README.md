# Exacto

Exacto (**EX**acto **A**utomated **C**aller for **T**ransformations in gen**O**mes / transcript**O**mes) 
identifies variants in aligned long-read DNA and RNA reads.

[![build](https://github.com/pirl-unc/exacto/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/exacto/actions/workflows/main.yml)

## 01. Dependencies
- python3 (3.10 tested)
- numpy (>=1.22.3)
- pandas (>=2.0.3)
- polars (>=1.12.0)
- pyarrow (>=18.0.0)
- pysam (>=0.22.0)
- rust

## 02. Installation
```
conda create -n exacto python=3.10
conda activate exacto
pip install exacto-<version>.tar.gz --verbose
```

## 03. Usage

Calling DNA variants in a long-read BAM:
```
exacto call-dna-vars [-h] 
    --bam-file BAM_FILE 
    --bam-bai-file BAM_BAI_FILE 
    --mode {all,case-specific} 
    --output-tsv-file OUTPUT_TSV_FILE 
    [--control-bam-files CONTROL_BAM_FILES [CONTROL_BAM_FILES ...]]
    [--control-bam-bai-files CONTROL_BAM_BAI_FILES [CONTROL_BAM_BAI_FILES ...]]
    [--num-threads NUM_THREADS]
    [--gzip GZIP]
    [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]
    [--min-reads MIN_READS]
    [--min-mapping-quality MIN_MAPPING_QUALITY]
    [--min-average-base-quality MIN_AVERAGE_BASE_QUALITY]
    [--min-size-proportion MIN_SIZE_PROPORTION]
    [--max-ins-norm-edit-distance MAX_INS_NORM_EDIT_DISTANCE]
    [--max-intrachromosomal-distance MAX_INTRACHROMOSOMAL_DISTANCE]
    [--max-intrachromosomal-distance-tau MAX_INTRACHROMOSOMAL_DISTANCE_TAU]
    [--max-interchromosomal-distance MAX_INTERCHROMOSOMAL_DISTANCE]
    [--apply-infinite-sites-assumption APPLY_INFINITE_SITES_ASSUMPTION]
    [--temp-dir TEMP_DIR]
```

Calling RNA variants in a long-read BAM (transcriptome assembly):
```
exacto call-rna-vars [-h] 
    --bam-file BAM_FILE 
    --bam-bai-file BAM_BAI_FILE 
    --reference-genome-fasta-file REFERENCE_GENOME_FASTA_FILE 
    --gene-annotation-file GENE_ANNOTATION_FILE 
    --gene-annotation-source GENE_ANNOTATION_SOURCE
     --output-exons-tsv-file OUTPUT_EXONS_TSV_FILE 
     --output-sj-tsv-file OUTPUT_SJ_TSV_FILE 
     --output-variants-tsv-file OUTPUT_VARIANTS_TSV_FILE 
     [--num-threads NUM_THREADS]
     [--gzip GZIP]
     [--min-mapping-quality MIN_MAPPING_QUALITY]
     [--min-average-base-quality MIN_AVERAGE_BASE_QUALITY]
     [--temp-dir TEMP_DIR]
```

## 04. DNA / RNA Variant Types Identified by Exacto

### DNA

- Single-nucleotide variant
- Insertion
- Deletion
- Translocation
- Breakpoint (inversion, duplication)

### RNA

- Single-nucleotide variant
- Insertion
- Deletion
- Breakpoint
- Alternative 5' and 3' splice sites
- Exon skipping
- Fusion gene
- Intron retention
- Cryptic exon
