# Exacto

Exacto (**E**Xacto **A**utomated **C**aller for **T**ransformations in gen**O**mes / transcript**O**mes) 
identifies variants in aligned long-read DNA and RNA reads.

[![build](https://github.com/pirl-unc/exacto/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/exacto/actions/workflows/main.yml)

## 01. Dependencies
- python3 (3.10 tested)
- pandas
- pysam
- rust

## 02. Installation
```
pip install . --verbose
```

## 03. Usage

Calling DNA variants in long-read BAM:
```
exacto call-dna-vars [-h] 
    --bam-file BAM_FILE 
    --sample-id SAMPLE_ID 
    --output-tsv-file OUTPUT_TSV_FILE 
    [--num-threads NUM_THREADS] 
    [--chromosomes CHROMOSOMES [CHROMOSOMES ...]] 
    [--min-reads MIN_READS]
    [--min-mapping-quality MIN_MAPPING_QUALITY]
    [--min-ins-size-proportion MIN_INS_SIZE_PROPORTION]
    [--max-ins-norm-edit-distance MAX_INS_NORM_EDIT_DISTANCE]
    [--min-del-size-proportion MIN_DEL_SIZE_PROPORTION]
    [--max-bnd-distance MAX_BND_DISTANCE]
    [--clustering-grid-size CLUSTERING_GRID_SIZE]
```

Calling RNA variants in long-read BAM:
```
exacto call-rna-vars [-h] 
    --bam-file BAM_FILE 
    --sample-id SAMPLE_ID 
    --output-tsv-file OUTPUT_TSV_FILE 
    [--num-threads NUM_THREADS] 
    [--chromosomes CHROMOSOMES [CHROMOSOMES ...]] 
    [--min-reads MIN_READS]
    [--min-mapping-quality MIN_MAPPING_QUALITY]
    [--min-ins-size-proportion MIN_INS_SIZE_PROPORTION]
    [--max-ins-norm-edit-distance MAX_INS_NORM_EDIT_DISTANCE]
    [--min-del-size-proportion MIN_DEL_SIZE_PROPORTION]
    [--max-bnd-distance MAX_BND_DISTANCE]
    [--clustering-grid-size CLUSTERING_GRID_SIZE]
```

## 04. Example

```
exacto call-rna-vars \
    --bam-file BAM_FILE \
    --output-tsv-file OUTPUT_TSV_FILE \
    --num-threads 8 \
    --chromosomes chr1 chr2 chr3 \
    --min-reads 3
```

## 05. DNA / RNA Variant Types Identified by Exacto

### DNA

- Single-nucleotide variant
- Insertion
- Deletion

### RNA

- Single-nucleotide variant
- Insertion
- Deletion
- Splicing
