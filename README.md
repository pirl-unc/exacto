# Exacto

**EX**acto **A**ccurate **C**haracterization of **T**ranscriptomes and gen**O**mes

A long-read toolkit for mutant proteoform prediction. Exacto identifies somatic
and germline DNA variants, RNA variants, integrates them, and translates
full-length transcripts with variant annotation at the amino-acid level.

[![CI](https://github.com/pirl-unc/exacto/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/pirl-unc/exacto/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Documentation**: [https://pirl-unc.github.io/exacto/](https://pirl-unc.github.io/exacto/)

## 01. Installation

Download the latest stable release [here](https://github.com/pirl-unc/exacto/releases).

```bash
conda create -n exacto python=3.10
conda activate exacto
pip install pysam==0.23.0
conda install -c conda-forge rust==1.88.0
conda install -c anaconda pandas==2.2.3
conda install -c conda-forge polars==1.26.0
conda install -c conda-forge pyarrow==19.0.1
pip install exacto-<version>.tar.gz --verbose
```

A Docker image is also available on
[Docker Hub](https://hub.docker.com/r/ajslee/exacto).

## 02. Dependencies

* Python (>=3.10)
* Rust (1.88.0 tested)
* numpy (>=1.22.3)
* pandas (>=2.0.3)
* polars (>=1.12.0)
* pyarrow (>=18.0.0)
* pysam (>=0.22.0)
* pytz (>=2024.1)

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

| Subcommand                      | Description                                                  |
|---------------------------------|--------------------------------------------------------------|
| `annotate-vars`                 | Annotate DNA or RNA variants with gene-level context.        |
| `build-genome-var-graph`        | Build a personalized genome variation graph.                 |
| `build-transcriptome-var-graph` | Build a personalized transcriptome variation graph.          |
| `call-germline-dna-vars`        | Call germline DNA variants from long-read alignments.        |
| `call-somatic-dna-vars`         | Call somatic DNA variants against matched control samples.   |
| `call-rna-vars`                 | Call RNA variants from assembled transcript alignments.      |
| `call-peptide-vars`             | Call peptide-level variants from translated proteoforms.     |
| `integrate-vars`                | Integrate DNA and RNA variants into a unified callset.       |
| `remove-unspliced-rnas`         | Filter out unspliced (nascent) RNAs from a transcriptome assembly. |
| `translate-seqs`                | Translate transcript sequences into peptide sequences.       |
| `translate-structs`             | Translate transcript structures into mutant proteoforms.     |

See the [Commands documentation](https://pirl-unc.github.io/exacto/cli/) for full
parameter documentation, and the [Pipelines documentation](https://pirl-unc.github.io/exacto/pipelines/)
for end-to-end mutant-proteoform-prediction and variation-graph-construction
walkthroughs.

## 04. License

Licensed under the Apache License, Version 2.0.
