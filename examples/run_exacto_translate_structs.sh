exacto translate-structs \
  --transcript-structures-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv \
  --rna-variant-calls-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv \
  --integrated-variants-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor_dna-001-tumor_variants_integrated.tsv \
  --strategy longest_orf \
  --num-threads 2 \
  --output-tsv-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_primary_structures.tsv \
  --output-fasta-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_primary_structures.fasta
