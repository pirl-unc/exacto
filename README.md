# Exacto

**EX**acto **A**ccurate **C**haracterization of **T**ranscriptomes and gen**O**mes

Exacto is a long-read toolkit for mutant proteoform prediction. Using long-read DNA and RNA sequencing data, Exacto:

* Identifies somatic and germline DNA variants.
* Performs variant-aware assembly of full-length transcripts.
* Identifies RNA variants.
* Links DNA and RNA variants.
* Translates transcripts into predicted proteoforms with amino acid-level annotations of DNA and RNA variants.
* Identifies mutant peptides against a reference proteome.

[![CI](https://github.com/pirl-unc/exacto/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/pirl-unc/exacto/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Documentation**: [https://pirl-unc.github.io/exacto/](https://pirl-unc.github.io/exacto/)

## 01. Installation

Download the latest stable release [here](https://github.com/pirl-unc/exacto/releases).

```bash
conda create -n exacto python=3.10
conda activate exacto
pip install pysam==0.23.0
conda install -c conda-forge rust==1.95.0
conda install -c conda-forge clangdev=21.1.8
conda install -c anaconda pandas==2.2.3
conda install -c conda-forge polars==1.26.0
conda install -c conda-forge pyarrow==19.0.1
pip install exacto-<version>.tar.gz --verbose
```

A Docker image is also available on
[Docker Hub](https://hub.docker.com/r/ajslee/exacto).

## 02. Dependencies

* Python (>=3.10)
* Rust (1.95.0 tested)
* libclang (>=9.0 and <22, build-time only — required by `bindgen` for the abPOA bindings)
* numpy (>=1.22.3)
* pandas (>=2.0.3)
* polars (>=1.12.0)
* pyarrow (>=18.0.0)
* pysam (>=0.22.0)
* pytz (>=2024.1)

Alignments must be made with [minimap2](https://github.com/lh3/minimap2) using `--cs`: exacto reads
variants from the `cs` tag.

## 03. Usage

### View all available subcommands
```bash
exacto --help
```

### View a subcommand's parameters
```bash
exacto <subcommand> --help
```

### Available subcommands

| Subcommand                      | Description                                                                 |
|---------------------------------|-----------------------------------------------------------------------------|
| `cluster-rna-reads`             | Cluster long RNA reads by transcript of origin and call variants in each cluster. |
| `correct-rna-reads`             | Error-correct RNA reads against the cluster each was assigned to.           |
| `determine-rna-consensus`       | Build one consensus sequence for each RNA read cluster.                     |
| `stitch-reference-transcripts`  | Restore degraded transcript ends from the best-matching reference transcripts. |
| `quantify-rna-abundances`       | Estimate per-cluster and per-reference-transcript abundances.               |
| `remove-unspliced-rnas`         | Remove unspliced RNAs from an RNA alignment.                                |
| `call-germline-dna-vars`        | Call germline DNA variants from long-read alignments.                       |
| `call-somatic-dna-vars`         | Call somatic DNA variants against one or more control samples.              |
| `call-rna-transcript-vars`      | Call RNA variants, transcript structures and nonsense-mediated decay predictions for assembled transcripts. |
| `annotate-vars`                 | Annotate DNA or RNA variants with gene-level context.                       |
| `integrate-vars`                | Link RNA variants to DNA variants based on linear genomic / transcriptomic distances. |
| `translate-transcripts`         | Translate assembled transcripts to proteoforms with DNA and RNA variants for each amino acid.     |
| `translate-seqs`                | Translate sequences in a FASTA / FASTQ file, or one sequence, to peptides.    |
| `call-peptide-vars`             | Call mutant peptides by comparing proteoforms with a reference proteome.    |
| `build-genome-var-graph`        | Build a personalized genome variation graph.                                |
| `build-transcriptome-var-graph` | Build a personalized transcriptome variation graph.                         |

See the [Commands documentation](https://pirl-unc.github.io/exacto/cli/) for full
parameter documentation, and the [Pipelines documentation](https://pirl-unc.github.io/exacto/pipelines/)
for end-to-end mutant-proteoform-prediction and variation-graph-construction
walkthroughs.

## 04. License

Licensed under the Apache License, Version 2.0.
