gtime -v \
  exacto integrate-vars \
    --annotated-dna-vars-tsv-file ../test/data/tsv/variant_callset/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv \
    --rna-vars-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor_minimap2_mdtagged_sorted_exacto_variants.tsv \
    --reference-gene-annotation-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
    --reference-gene-annotation-source gencode \
    --output-tsv-file outputs/rna-100_dna-001_integrated_callset.tsv
