EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/call-germline-dna-vars/

exacto call-germline-dna-vars \
  --bam-file ${EXACTO_TEST_DATA}/alignment/scga-mini-dna-001-normal_minimap2_sorted.bam \
  --bai-file ${EXACTO_TEST_DATA}/alignment/scga-mini-dna-001-normal_minimap2_sorted.bam.bai \
  --fasta-file ${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz \
  --output-tsv-file outputs/call-germline-dna-vars/scga-mini-dna-001-normal_exacto_germline_dna_variants.tsv \
  --preset pb
