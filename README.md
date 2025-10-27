# Exacto

Exacto (**EX**acto **A**ccurate **C**haracterization of **T**ranscriptomes and gen**O**mes) performs the following 
primary tasks using long-read sequencing data for mutant proteoform prediction:
* Identification of somatic and germline DNA variants.
* Identification of RNA variants.
* Integration of DNA and RNA variants.
* Translation of full-length transcripts with any underlying DNA and RNA variant annotation at the amino-acid level.

[![build](https://github.com/pirl-unc/exacto/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/exacto/actions/workflows/main.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## 01. Docker Container
Docker images of Exacto can be found here: <br/>
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

| Command                 | Description                                     |
|-------------------------|-------------------------------------------------|
| `annotate-vars`         | Annotate variants                               | 
| `call-dna-vars`         | Perform somatic or germline DNA variant calling |
| `call-rna-vars`         | Perform RNA variant calling                     |
| `integrate-vars`        | Integrate DNA and RNA variants                  |
| `remove-unspliced-rnas` | Removed unspliced RNAs                          |
| `translate-seqs`        | Translate transcript sequences                  |
| `translate-structs`     | Translate transcript structures                 |

Example scripts for running Exacto can be found [here](https://github.com/pirl-unc/exacto/tree/main/examples).

## 05. Mutant Peptide Prediction Pipeline

Align tumor and normal DNA reads to a reference genome:

```
minimap2 \
    -ax map-hifi --cs --eqx -Y -L --secondary=no \
    reference_genome.fasta tumor_dna.fastq.gz \
| samtools view -b - \
| samtools sort -o {tumor,normal}_dna.sorted.bam
```

Identify tumor-specific (somatic) DNA variants:
```
exacto call-dna-vars \
    --bam-file tumor_dna.sorted.bam \
    --bam-bai-file tumor_dna.sorted.bam.bai \
    --mode case-specific \
    --control-bam-files normal_dna.sorted.bam \
    --control-bam-bai-files normal_dna.sorted.bam.bai \
    --output-tsv-file tumor_specific_dna_variants.tsv
```

Annotate the tumor-specific (somatic) DNA variants:

```
exacto annotate-vars \
    --tsv-file tumor_specific_dna_variants.tsv \
    --reference-gene-annotation-file gencode.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly assembly_name \
    --reference-gene-annotation-version assembly_version \
    --output-tsv-file tumor_specific_dna_variants.annotated.tsv
```

Assemble tumor transcriptome using [RNAbloom2](https://github.com/bcgsc/RNA-Bloom):

```
java -jar RNA-Bloom.jar \
    -long tumor_rna.fastq.gz \
    --outdir <rnabloom2_outputs> \
    --qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf
```

Align the assembled tumor transcriptome to a reference genome using [Minimap2](https://github.com/lh3/minimap2):

```
minimap2 \
    -ax splice:hq -uf --cs --eqx -Y -L --secondary=no \
    reference_genome.fasta tumor_rna.fastq.gz \
| samtools view -b - \
| samtools sort -o tumor_transcriptome_assembly.sorted.bam
```

Remove unspliced models from the transcriptome assembly:

```
exacto remove-unspliced-rnas \
    --bam-file tumor_transcriptome_assembly.sorted.bam \
    --bam-bai-file tumor_transcriptome_assembly.sorted.bam.bai \
    --fasta-file reference_genome.fasta \
    --reference-gene-annotation-file gencode.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly assembly_name \
    --reference-gene-annotation-version assembly_version \
    --output-bam-file tumor_transcriptome_assembly.sorted.filtered.bam \
    --output-bam-bai-file tumor_transcriptome_assembly.sorted.filtered.bam.bai \ 
    --output-fasta-file tumor_transcriptome_assembly.sorted.filtered.fasta
```

Identify tumor RNA variants:

```
exacto call-rna-vars \
    --bam-file tumor_transcriptome_assembly.sorted.filtered.bam \
    --bam-bai-file tumor_transcriptome_assembly.sorted.filtered.bam.bai \
    --reference-genome-fasta-file reference_genome.fasta \
    --reference-gene-annotation-file gencode.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly assembly_name \
    --reference-gene-annotation-version assembly_version \
    --output-dir rna_variants_outputs/ \ 
    --output-prefix tumor
```

Integrate DNA and RNA variants:

```
exacto integrate-vars \ 
    --annotated-dna-vars-tsv-file tumor_specific_dna_variants.annotated.tsv \
    --rna-vars-tsv-file rna_variants_outputs/tumor_exacto_rna_variant_calls.tsv \
    --reference-gene-annotation-file gencode.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly assembly_name \
    --reference-gene-annotation-version assembly_version \
    --output-tsv-file tumor_dna_rna_variants_integrated.tsv
```

Identify primary structures:

```
exacto translate-structs \
    --transcript-structures-tsv-file rna_variants_outputs/tumor_exacto_transcript_structures.tsv \ 
    --rna-variant-calls-tsv-file rna_variants_outputs/tumor_exacto_rna_variant_calls.tsv \
    --integrated-variants-tsv-file tumor_dna_rna_variants_integrated.tsv \
    --strategy longest_orf \
    --output-tsv-file tumor_primary_structures.tsv
```

### Case-specific DNA Variant Calling

Identify case-specific (e.g. somatic) DNA variants in a case (e.g. tumor) long-read DNA BAM file against a set of control (e.g. matched normal) long-read DNA BAM files:

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

### DNA Variant Calling

Identify all (e.g. germline) DNA variants in a long-read DNA BAM file:

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

### RNA Variant Calling

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

### DNA and RNA Variant Integration

Integrate DNA and RNA variants based on genomic coordinates.

```
exacto integrate-vars [-h] 
    --annotated-dna-vars-tsv-file ANNOTATED_DNA_VARS_TSV_FILE 
    --rna-vars-tsv-file RNA_VARS_TSV_FILE 
    --reference-gene-annotation-file REFERENCE_GENE_ANNOTATION_FILE 
    --reference-gene-annotation-source REFERENCE_GENE_ANNOTATION_SOURCE 
    --reference-gene-annotation-assembly REFERENCE_GENE_ANNOTATION_ASSEMBLY 
    --reference-gene-annotation-version REFERENCE_GENE_ANNOTATION_VERSION 
    --output-tsv-file OUTPUT_TSV_FILE
    [--num-threads NUM_THREADS] 
    [--max-exon-offset MAX_EXON_OFFSET] 
    [--max-transcript-boundary-offset MAX_TRANSCRIPT_BOUNDARY_OFFSET]
    [--max-intergenic-distance MAX_INTERGENIC_DISTANCE]
```

### Unspliced Transcript Removal

Remove unspliced (nascent) RNAs in a transcriptome assembly. 
Note that assembled transcripts whose alignments overlap with any single-exon reference transcript will be retained.

```
exacto remove-unspliced-rnas [-h] 
    --bam-file BAM_FILE 
    --bam-bai-file BAM_BAI_FILE 
    --fasta-file FASTA_FILE 
    --reference-gene-annotation-file REFERENCE_GENE_ANNOTATION_FILE 
    --reference-gene-annotation-source REFERENCE_GENE_ANNOTATION_SOURCE
    --reference-gene-annotation-assembly REFERENCE_GENE_ANNOTATION_ASSEMBLY 
    --reference-gene-annotation-version REFERENCE_GENE_ANNOTATION_VERSION 
    --output-bam-file OUTPUT_BAM_FILE 
    --output-bam-bai-file OUTPUT_BAM_BAI_FILE 
    --output-fasta-file OUTPUT_FASTA_FILE 
    [--num-threads NUM_THREADS]
    [--min-mapping-quality MIN_MAPPING_QUALITY]
    [--gene-types GENE_TYPES [GENE_TYPES ...]]
    [--gene-levels GENE_LEVELS [GENE_LEVELS ...]]
    [--transcript-types TRANSCRIPT_TYPES [TRANSCRIPT_TYPES ...]]
    [--transcript-levels TRANSCRIPT_LEVELS [TRANSCRIPT_LEVELS ...]]
```

### Transcript Sequence Translation

Translate transcript sequences to peptide sequences.

```
exacto translate-seqs [-h] 
    (--fastq-file FASTQ_FILE | --fasta-file FASTA_FILE | --sequence SEQUENCE) 
    --strategy {longest_orf,all_orfs} 
    [--output-tsv-file OUTPUT_TSV_FILE]
    [--output-fasta-file OUTPUT_FASTA_FILE]
    [--num-threads NUM_THREADS]
    [--temp-dir TEMP_DIR]
    [--gzip GZIP]
```

### Transcript Structure Translation

Translate transcript structures to primary structures.

```
exacto translate-structs [-h]
    --transcript-structures-tsv-file TRANSCRIPT_STRUCTURES_TSV_FILE 
    --rna-variant-calls-tsv-file RNA_VARIANT_CALLS_TSV_FILE
    --integrated-variants-tsv-file INTEGRATED_VARIANTS_TSV_FILE 
    --strategy {longest_orf,all_orfs} 
    --output-tsv-file OUTPUT_TSV_FILE
    [--num-threads NUM_THREADS]
```

## 06. Input Preparation

### DNA Variant Calling

Exacto performs DNA variant identification using the [cs tag](https://github.com/lh3/minimap2#cs) produced by [minimap2](https://github.com/lh3/minimap2) alignments.

For automated long-read DNA alignment, you can use the Nexus workflow manager:<br/>
https://github.com/pirl-unc/nexus/tree/main/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2

If you prefer to run `minimap2` outside of Nexus, use the following parameters:
```
-ax map-hifi --cs --eqx -Y -L --secondary=no
```

### RNA Variant Calling

Exacto identifies RNA variants from assembled transcript models rather than directly from raw reads. 
You can generate transcript models from long RNA-seq reads using [RNAbloom2](https://github.com/bcgsc/RNA-Bloom). A corresponding `Nexus` workflow for RNAbloom2 is available here:<br/>
https://github.com/pirl-unc/nexus/tree/main/src/nexuslib/pipelines/assembly/transcriptome_assembly_rnabloom2

After assembling the transcriptome, align the assembled transcripts back to the reference genome using `minimap2` with the following parameters:
```
-ax splice:hq -uf --cs --eqx -Y -L --secondary=no
```

Then, remove any assembled transcripts that are likely unspliced RNAs by running `exacto remove-unspliced-rnas`.

## 07. DNA / RNA Variant Types Identified by Exacto

### DNA (Somatic and Germline)

Sequence variant types:
- Breakpoint (duplication, inversion)
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
- Exon truncation
- Fusion gene
- Intron retention
- UTR extension
