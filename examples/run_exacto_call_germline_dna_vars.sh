mkdir -p outputs/call-germline-dna-vars/

exacto call-germline-dna-vars \
  --bam-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam \
  --bam-bai-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai \
  --fasta-file ../test/data/fasta/hg38_chr17-18.fa.gz \
  --output-tsv-file outputs/call-germline-dna-vars/dna-001-tumor_minimap2_mdtagged_sorted_exacto_variants.tsv
