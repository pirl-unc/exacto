mkdir -p outputs/call-somatic-dna-vars/

exacto call-somatic-dna-vars \
  --bam-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam \
  --bam-bai-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai \
  --control-bam-files ../test/data/bam/dna-001-normal_minimap2_mdtagged_sorted.bam \
  --control-bam-bai-files ../test/data/bam/dna-001-normal_minimap2_mdtagged_sorted.bam.bai \
  --fasta-file ../test/data/fasta/hg38_chr17-18.fa.gz \
  --output-tsv-file outputs/call-somatic-dna-vars/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv
