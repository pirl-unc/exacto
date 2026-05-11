mkdir -p outputs/call-peptide-vars/

exacto call-peptide-vars \
  --primary-structures-tsv-file ../test/data/tsv/primary_structure/rna-100-tumor_minimap2_mdtagged_sorted_exacto_primary_structures.tsv \
  --reference-fasta-file ../test/data/fasta/reference_peptides.fasta \
  --output-tsv-file outputs/call-peptide-vars/rna-100-tumor_minimap2_mdtagged_sorted_exacto_peptide_variants.tsv \
  --output-fasta-file outputs/call-peptide-vars/rna-100-tumor_minimap2_mdtagged_sorted_exacto_peptide_variants.fasta \
  --min-k 8 \
  --max-k 11 \
  --num-threads 1
