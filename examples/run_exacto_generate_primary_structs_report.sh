mkdir -p outputs/generate-primary-structs-report/

exacto generate-primary-structs-report \
  --dna-vars-tsv-file ../test/data/tsv/variant_callset/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv \
  --rna-vars-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv \
  --transcripts-read-support-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcripts_read_support.tsv \
  --transcriptome-assembly-read-support-tsv-file ../test/data/tsv/rna_assembly/sample_rna_assembly_read_support.tsv \
  --integrated-vars-tsv-file ../test/data/tsv/variant_callset/rna-100-tumor_dna-001-tumor_variants_integrated.tsv \
  --primary-structures-tsv-file ../test/data/tsv/primary_structure/rna-100-tumor_minimap2_mdtagged_sorted_exacto_primary_structures.tsv \
  --output-tsv-file outputs/generate-primary-structs-report/sample_primary_structs_report.tsv
