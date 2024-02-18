# Exacto

Exacto (**E**xacto **A**lterations and **C**himeric **T**ranscripts **O**perations) identifies variant 
transcripts in aligned long-read RNA reads.

**E**Xacto **A**ccurate **C**alling of **T**ransformations in **O**ncology
**E**Xacto **A**ccurate **C**aller for  ****ransformations in genomes and transcriptomes

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
```
exacto identify [-h]
    --bam-file BAM_FILE 
    --output-tsv-file OUTPUT_TSV_FILE 
    [--num-threads NUM_THREADS]
    [--chromosomes CHROMOSOMES]
    [--min-reads MIN_READS]
    [--min-mapping-quality MIN_MAPPING_QUALITY]
    [--min-ins-size-proportion MIN_INS_SIZE_PROPORTION]
    [--max-ins-norm-edit-distance MAX_INS_NORM_EDIT_DISTANCE]
    [--min-del-size-proportion MIN_DEL_SIZE_PROPORTION]
```

## 04. Example

```
exacto identify \
    --bam-file BAM_FILE \
    --output-tsv-file OUTPUT_TSV_FILE \
    --num-threads 8 \
    --chromosomes chr1 chr2 chr3 \
    --min-reads 3
```
