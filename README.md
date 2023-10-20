# Exacto

Exacto (EXtracting And Counting Transcripts in Oncology) identifies and quantifies augmented transcripts in tumor samples.

## 01. Dependencies
- python3 (3.9 tested)
- pandas
- pysam
- pyvcf
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
