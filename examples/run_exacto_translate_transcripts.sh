EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/translate-transcripts/

exacto translate-transcripts \
  --assembled-transcript-model-alignments-tsv-file ${EXACTO_TEST_DATA}/variant_calling/rna/scga-mini-rna-001-tumor_exacto_assembled_transcript_model_alignments.tsv \
  --assembled-transcript-variants-tsv-file ${EXACTO_TEST_DATA}/variant_calling/rna/scga-mini-rna-001-tumor_exacto_assembled_transcript_variants.tsv \
  --dna-variants-tsv-file ${EXACTO_TEST_DATA}/variant_calling/dna/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv \
  --integrated-variants-tsv-file ${EXACTO_TEST_DATA}/integration/scga-mini-rna-001-tumor_exacto_integrated_variants.tsv \
  --rna-consensus-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_consensus.tsv \
  --stitched-transcripts-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/stitched/scga-mini-rna-001-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts.tsv \
  --strategy longest_orf \
  --num-threads 2 \
  --output-dir outputs/translate-transcripts/ \
  --output-prefix scga-mini-rna-001-tumor
