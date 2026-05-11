mkdir -p outputs/translate-seqs/

exacto translate-seqs \
  --fasta-file ../test/data/fasta/sample_tumor_long_read_rna.fasta \
  --strategy longest_orf \
  --num-threads 4 \
  --output-tsv-file outputs/translate-seqs/sample_tumor_long_read_rna_exacto_translations_longest-orf.tsv \
  --output-fasta-file outputs/translate-seqs/sample_tumor_long_read_rna_exacto_translations_longest-orf.fasta \
  --gzip yes

exacto translate-seqs \
  --fasta-file ../test/data/fasta/sample_tumor_long_read_rna.fasta \
  --strategy all_orfs \
  --num-threads 4 \
  --output-tsv-file outputs/translate-seqs/sample_tumor_long_read_rna_exacto_translations_all-orfs.tsv \
  --output-fasta-file outputs/translate-seqs/sample_tumor_long_read_rna_exacto_translations_all-orfs.fasta \
  --gzip yes

exacto translate-seqs \
  --sequence ATGGGGCCCATGCCTTAG \
  --strategy longest_orf

exacto translate-seqs \
  --sequence ATGGGGCCCATGCCTTAG \
  --strategy all_orfs
