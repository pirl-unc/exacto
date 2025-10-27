exacto call-dna-vars \
  --bam-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam \
  --bam-bai-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai \
  --mode all \
  --output-tsv-file outputs/dna-001-tumor_minimap2_mdtagged_sorted_exacto_variants.tsv

exacto call-dna-vars \
  --bam-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam \
  --bam-bai-file ../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai \
  --control-bam-files ../test/data/bam/dna-001-normal_minimap2_mdtagged_sorted.bam \
  --control-bam-bai-files ../test/data/bam/dna-001-normal_minimap2_mdtagged_sorted.bam.bai \
  --mode case-specific \
  --output-tsv-file outputs/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv
