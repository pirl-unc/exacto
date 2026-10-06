EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/call-somatic-dna-vars/

exacto call-somatic-dna-vars \
  --bam-file ${EXACTO_TEST_DATA}/alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam \
  --bai-file ${EXACTO_TEST_DATA}/alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai \
  --control-bam-files ${EXACTO_TEST_DATA}/alignment/scga-mini-dna-001-normal_minimap2_sorted.bam \
  --control-bai-files ${EXACTO_TEST_DATA}/alignment/scga-mini-dna-001-normal_minimap2_sorted.bam.bai \
  --fasta-file ${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz \
  --output-tsv-file outputs/call-somatic-dna-vars/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv \
  --preset pb
