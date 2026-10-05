EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/call-peptide-vars/

exacto call-peptide-vars \
  --proteoforms-tsv-file ${EXACTO_TEST_DATA}/translation/scga-mini-rna-001-tumor_exacto_proteoforms.tsv \
  --reference-fasta-file ${EXACTO_TEST_DATA}/exacto/exactolib/reference_peptides.fasta \
  --output-tsv-file outputs/call-peptide-vars/scga-mini-rna-001-tumor_exacto_peptide_variants.tsv \
  --output-fasta-file outputs/call-peptide-vars/scga-mini-rna-001-tumor_exacto_peptide_variants.fasta \
  --min-k 8 \
  --max-k 11 \
  --num-threads 1
