exacto integrate-vars \
  --annotated-dna-vars-tsv-file ../test/data/tsv/variant_callset/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv \
  --rna-vars-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv \
  --reference-gene-annotation-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-tsv-file outputs/rna-100_dna-001_integration.tsv
