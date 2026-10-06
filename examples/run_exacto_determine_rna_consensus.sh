EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/determine-rna-consensus/

exacto determine-rna-consensus \
  --tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_clusters.tsv \
  --fastq-file ${EXACTO_TEST_DATA}/assembly/rna/pass1/scga-mini-rna-001-tumor_pass1_exacto_rna_corrected_reads.fastq.gz \
  --output-dir outputs/determine-rna-consensus/ \
  --output-prefix scga-mini-rna-001-tumor_pass2
