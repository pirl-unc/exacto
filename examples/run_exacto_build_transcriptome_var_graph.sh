mkdir -p outputs/build-transcriptome-var-graph/

exacto build-transcriptome-var-graph \
  --transcript-structures-tsv-file ../test/data/tsv/transcript_structure/sample_transcript_structure.tsv \
  --fasta-file ../test/data/fasta/sample3.fa \
  --output-fasta-file outputs/build-transcriptome-var-graph/sample_transcriptome.fasta