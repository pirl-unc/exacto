# Exacto

Exacto (EXtracting And Counting Transcripts in Oncology) identifies and quantifies wildtype and variant transcripts in tumor samples.

[![build](https://github.com/pirl-unc/exacto/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/exacto/actions/workflows/main.yml)

## 01. Dependencies
- python3 (3.10 tested)
- pandas
- pysam
- pyensembl

## 02. Installation
```
pip install . --ignore-installed --verbose
```

## 03. Packaging
```
python -m build
```

## 04. Testing
```
bash lint.sh
bash unittest.sh
```
