exacto translate \
  --fasta-file ../test/data/fasta/sample_tumor_long_read_rna.fasta \
  --strategy longest-orf \
  --num-threads 4 \
  --output-tsv-file outputs/sample_tumor_long_read_rna_exacto_translations.tsv \
  --output-fasta-file outputs/sample_tumor_long_read_rna_exacto_translations.fasta \
  --gzip yes
