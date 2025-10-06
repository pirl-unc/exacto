# Exacto

Exacto (**EX**acto **A**ccurate **C**haracterization of **T**ranscriptomes and gen**O**mes) 
identifies genomic and transcriptomic variants from reference-aligned long-read DNA and RNA sequences.

[![build](https://github.com/pirl-unc/exacto/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/exacto/actions/workflows/main.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## 01. Docker Container
Docker images for Exacto can be found here: <br/>
https://hub.docker.com/r/ajslee/exacto

## 02. Dependencies
- python3 (3.10 tested)
- numpy (>=1.22.3)
- pandas (>=2.0.3)
- polars (>=1.12.0)
- pyarrow (>=18.0.0)
- pysam (>=0.22.0)
- pytz (>=2024.1)
- rust

## 03. Installation
```
conda create -n exacto python=3.10
conda activate exacto
pip install pysam==0.23.0
conda install -c conda-forge rust==1.86.0
conda install -c anaconda pandas==2.2.3
conda install -c conda-forge polars==1.26.0
conda install -c conda-forge pyarrow==19.0.1
pip install exacto-<version>.tar.gz --verbose
```

## 04. Usage

| Command         | Description                                     |
|-----------------|-------------------------------------------------|
| `annotate`      | Annotate variants                               | 
| `call-dna-vars` | Perform somatic or germline DNA variant calling |
| `call-rna-vars` | Perform RNA variant calling                     |
| `translate`     | Translate full-length RNA reads                 |

Example scripts for running Exacto can be found [here](https://github.com/pirl-unc/exacto/tree/main/examples).

### Somatic DNA variant calling

Identify somatic DNA variants in a tumor long-read DNA BAM file against a set of normal long-read DNA BAM files:

```
exacto call-dna-vars [-h] 
    --bam-file BAM_FILE 
    --bam-bai-file BAM_BAI_FILE 
    --mode case-specific 
    --output-tsv-file OUTPUT_TSV_FILE 
    --control-bam-files CONTROL_BAM_FILES [CONTROL_BAM_FILES ...]
    --control-bam-bai-files CONTROL_BAM_BAI_FILES [CONTROL_BAM_BAI_FILES ...]
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

### Germline DNA variant calling

Identify germline DNA variants in a long-read DNA BAM file:

```
exacto call-dna-vars [-h] 
    --bam-file BAM_FILE 
    --bam-bai-file BAM_BAI_FILE 
    --mode all
    --output-tsv-file OUTPUT_TSV_FILE 
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

### RNA variant calling

Identify RNA variants in a long-read assembled transcripts BAM file:

```
exacto call-rna-vars [-h] 
    --bam-file BAM_FILE
    --bam-bai-file BAM_BAI_FILE
    --reference-genome-fasta-file REFERENCE_GENOME_FASTA_FILE
    --reference-gene-annotation-file REFERENCE_GENE_ANNOTATION_FILE 
    --reference-gene-annotation-source REFERENCE_GENE_ANNOTATION_SOURCE 
    --reference-gene-annotation-assembly REFERENCE_GENE_ANNOTATION_ASSEMBLY 
    --reference-gene-annotation-version REFERENCE_GENE_ANNOTATION_VERSION 
    --output-dir OUTPUT_DIR 
    --output-prefix OUTPUT_PREFIX 
    [--num-threads NUM_THREADS] 
    [--min-mapping-quality MIN_MAPPING_QUALITY]
    [--reference-transcript-scoring-method REFERENCE_TRANSCRIPT_SCORING_METHOD]
    [--reference-transcript-selection-strategy REFERENCE_TRANSCRIPT_SELECTION_STRATEGY] 
    [--reference-transcript-top-k REFERENCE_TRANSCRIPT_TOP_K]
    [--reference-transcript-threshold REFERENCE_TRANSCRIPT_THRESHOLD] 
    [--min-average-base-quality MIN_AVERAGE_BASE_QUALITY] 
    [--temp-dir TEMP_DIR]
    [--gene-types GENE_TYPES [GENE_TYPES ...]] 
    [--gene-levels GENE_LEVELS [GENE_LEVELS ...]] 
    [--transcript-types TRANSCRIPT_TYPES [TRANSCRIPT_TYPES ...]]
    [--transcript-levels TRANSCRIPT_LEVELS [TRANSCRIPT_LEVELS ...]]
```

### Translation

Translate RNA sequences to peptide sequences.

```
exacto translate [-h] 
    (--fastq-file FASTQ_FILE | --fasta-file FASTA_FILE | --sequence SEQUENCE) 
    --strategy {longest_orf,all_orfs} 
    [--output-tsv-file OUTPUT_TSV_FILE]
    [--output-fasta-file OUTPUT_FASTA_FILE]
    [--num-threads NUM_THREADS]
    [--temp-dir TEMP_DIR]
    [--gzip GZIP]
```

## 05. Input Preparation

### DNA variant calling

Exacto performs DNA variant identification using the [cs tag](https://github.com/lh3/minimap2#cs) produced by [minimap2](https://github.com/lh3/minimap2) alignments.

For automated long-read DNA alignment, you can use the Nexus workflow manager:<br/>
https://github.com/pirl-unc/nexus/tree/main/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2

If you prefer to run `minimap2` outside of Nexus, use the following parameters:
```
-ax map-hifi --cs --eqx -Y -L --secondary=no
```

### RNA variant calling

Exacto identifies RNA variants from assembled transcript models rather than directly from raw reads. 
You can generate transcript models from long RNA-seq reads using [RNAbloom2](https://github.com/bcgsc/RNA-Bloom). A corresponding `Nexus` workflow for RNAbloom2 is available here:<br/>
https://github.com/pirl-unc/nexus/tree/main/src/nexuslib/pipelines/assembly/transcriptome_assembly_rnabloom2

After assembling the transcriptome, align the asembled transcripts back to the reference genome using `minimap2` with the following parameters:
```
-ax splice:hq -uf --cs --eqx -Y -L --secondary=no
```

## 06. DNA / RNA Variant Types Identified by Exacto

### DNA (Somatic and Germline)

Sequence variant types:
- Breakpoint (inversion, duplication)
- Deletion
- Insertion
- Multi-nucleotide variant
- Single-nucleotide variant
- Translocation

### RNA

Sequence variant types:
- Breakpoint
- Deletion
- Insertion
- Multi-nucleotide variant
- Single-nucleotide variant

Splice variant types:
- Circular RNA
- Cryptic exon
- Exon skipping
- Fusion gene
- Intron retention
- UTR extension
