exacto diff-kmers \
  --query-fasta-file ../test/data/fasta/query_peptides.fasta \
  --reference-fasta-file ../test/data/fasta/reference_peptides.fasta \
  --min-k 8 \
  --max-k 11 \
  --output-tsv-file outputs/query_unique_kmer_peptides.tsv
