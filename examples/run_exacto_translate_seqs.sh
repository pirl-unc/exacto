EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/translate-seqs/

exacto translate-seqs \
  --fastx-file ${EXACTO_TEST_DATA}/assembly/rna/stitched/scga-mini-rna-001-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts.fasta \
  --strategy longest_orf \
  --num-threads 4 \
  --output-tsv-file outputs/translate-seqs/scga-mini-rna-001-tumor_exacto_stitched_transcripts_translations_longest-orf.tsv \
  --output-fasta-file outputs/translate-seqs/scga-mini-rna-001-tumor_exacto_stitched_transcripts_translations_longest-orf.fasta

exacto translate-seqs \
  --fastx-file ${EXACTO_TEST_DATA}/assembly/rna/stitched/scga-mini-rna-001-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts.fasta \
  --strategy all_orfs \
  --num-threads 4 \
  --output-tsv-file outputs/translate-seqs/scga-mini-rna-001-tumor_exacto_stitched_transcripts_translations_all-orfs.tsv \
  --output-fasta-file outputs/translate-seqs/scga-mini-rna-001-tumor_exacto_stitched_transcripts_translations_all-orfs.fasta

exacto translate-seqs \
  --sequence ATGGGGCCCATGCCTTAG \
  --strategy longest_orf

exacto translate-seqs \
  --sequence ATGGGGCCCATGCCTTAG \
  --strategy all_orfs
