exacto translate-structs \
  --transcript-structures-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv \
  --rna-variant-calls-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv \
  --integrated-variants-tsv-file ../test/data/tsv/variant_callset/rna-100_dna-001_integration.tsv \
  --strategy longest_orf \
  --num-threads 4 \
  --output-tsv-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_primary_structures.tsv
