mkdir -p outputs/translate-structs/

exacto translate-structs \
  --rna-assembly-support-tsv-file ../test/data/tsv/rna_assembly/rna-100_transcriptome_assembly_read_support.tsv \
  --transcript-model-structures-tsv-file ../test/data/tsv/transcript_model_structure/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv \
  --rna-variants-tsv-file ../test/data/tsv/variant_callset/integrate_variants/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_records.tsv \
  --dna-variants-tsv-file ../test/data/tsv/variant_callset/integrate_variants/dna-001-tumor_minimap2_mdtagged_sorted_exacto_dna_variant_records.tsv \
  --integrated-variants-tsv-file ../test/data/tsv/variant_callset/rna-100_dna-001_integration.tsv \
  --strategy longest_orf \
  --num-threads 2 \
  --output-dir outputs/translate-structs/ \
  --output-prefix rna-100-tumor_minimap2_mdtagged_sorted
