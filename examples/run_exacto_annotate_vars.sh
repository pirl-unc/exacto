gtime -v \
  exacto annotate-vars \
    --tsv-file outputs/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv \
    --reference-gene-annotation-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly hg38 \
    --reference-gene-annotation-version v41 \
    --output-tsv-file outputs/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv

