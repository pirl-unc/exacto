EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/build-transcriptome-var-graph/

exacto build-transcriptome-var-graph \
  --transcript-structures-tsv-file ${EXACTO_TEST_DATA}/variant_calling/rna/scga-mini-rna-001-tumor_exacto_assembled_transcript_model_alignments.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz \
  --output-fasta-file outputs/build-transcriptome-var-graph/scga-mini-rna-001-tumor_transcriptome_var_graph.fasta
