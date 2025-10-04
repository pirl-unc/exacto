exacto translate \
  --fasta-file ../test/data/fasta/sample_tumor_long_read_rna.fasta \
  --strategy longest_orf \
  --num-threads 4 \
  --output-tsv-file outputs/sample_tumor_long_read_rna_exacto_translations_longest-orf.tsv \
  --output-fasta-file outputs/sample_tumor_long_read_rna_exacto_translations_longest-orf.fasta \
  --gzip yes

exacto translate \
  --fasta-file ../test/data/fasta/sample_tumor_long_read_rna.fasta \
  --strategy all_orfs \
  --num-threads 4 \
  --output-tsv-file outputs/sample_tumor_long_read_rna_exacto_translations_all-orfs.tsv \
  --output-fasta-file outputs/sample_tumor_long_read_rna_exacto_translations_all-orfs.fasta \
  --gzip yes

exacto translate \
  --sequence ATGGGGCCCATGCCTTAG \
  --strategy longest_orf

exacto translate \
  --sequence ATGGGGCCCATGCCTTAG \
  --strategy all_orfs