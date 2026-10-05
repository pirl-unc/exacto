EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/quantify-rna-abundances/

exacto quantify-rna-abundances \
  --clusters-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_clusters.tsv \
  --cluster-reference-matches-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_clusters_reference_transcripts.tsv \
  --cluster-splice-junctions-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_clusters_splice_junctions.tsv \
  --cluster-variants-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_clusters_variants_passed.tsv \
  --consensus-rna-reference-matches-tsv-file ${EXACTO_TEST_DATA}/variant_calling/rna/scga-mini-rna-001-tumor_exacto_assembled_transcript_reference_transcript_matches.tsv \
  --output-dir outputs/quantify-rna-abundances/ \
  --output-prefix scga-mini-rna-001-tumor
