mkdir -p outputs/integrate-vars/

exacto integrate-vars \
  --dna-variants-tsv-file ../test/data/tsv/variant_callset/integrate_variants/dna-001-tumor_minimap2_mdtagged_sorted_exacto_dna_variant_records.tsv \
  --rna-variants-tsv-file ../test/data/tsv/variant_callset/integrate_variants/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_records.tsv \
  --reference-gene-annotation-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-tsv-file outputs/integrate-vars/rna-100-tumor_dna-001-tumor_variants_integrated.tsv
