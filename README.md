# Exacto

Expression Quantification and Identification of Augmented Transcripts of Somatic Variant Origin

## 01. Dependencies
- python3 (3.9 tested)
- pysam
- pyvcf
- pyensembl

## 02. Installation
```
python setup.py install
```

## 03. Usage
Refine a structural variant callset.
```
exacto_refine_sv_callset.py 
    --vcf_file <sample_id>.vcf 
    --sv_calling_method {sniffles,sniffles2,svim,cutesv,delly2,lumpy,pbsv}
    --sequencing_platform {illumina,pacbio-hifi-ccs,ont}
    --blacklisted_regions_tsv_file ucsc_hg38_gap_table.txt
    --gap_padding 1000000
    --filter_values_to_include PASS
    --min_total_coverage 7
    --min_variant_reads_count 3
    --keep_only_precise 1
    --output_tsv_file <sample_id>_refined.tsv 
    --chromosomes chr1 chr2 chr3
```

Annotate a structural variant callset.
```
exacto_annotate_sv_callset.py
    --ensembl_release 106 
    --tsv_file <sample_id>_refined.tsv 
    --output_tsv_file <sample_id>_annotated.tsv
```

Merge structural variant callsets.
```
exacto_merge_sv_callsets.py 
    --tsv_files <sample_id>_1.tsv <sample_id>_2.tsv <sample_id>_3.tsv 
    --output_merged_tsv_file <sample_id>_merged.tsv
    --output_merged_deduped_tsv_file <sample_id>_merged_deduped.tsv
    --methods_priority_list pacbio_sniffles2 pacbio_cutesv pacbio_svim
    --max_cluster_distance 10
```