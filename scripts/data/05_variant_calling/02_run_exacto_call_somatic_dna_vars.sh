run() {
  local sample_id="$1"
  local control_sample_id="${sample_id%-tumor}-normal"

  exacto call-somatic-dna-vars \
    --bam-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/alignment/${sample_id}_minimap2_sorted.bam \
    --bai-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/alignment/${sample_id}_minimap2_sorted.bam.bai \
    --control-bam-files /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/alignment/${control_sample_id}_minimap2_sorted.bam \
    --control-bai-files /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/alignment/${control_sample_id}_minimap2_sorted.bam.bai \
    --fasta-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz \
    --output-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/dna/${sample_id}_exacto_somatic_dna_variants.tsv \
    --preset pb \
    --num-threads 16
}

mkdir -p /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/dna/

sample_ids=(
    "scga-mini-dna-001-tumor"
    "scga-mini-dna-002-tumor"
    "scga-mini-dna-003-tumor"
    "scga-mini-dna-004-tumor"
    "scga-mini-dna-005-tumor"
    "scga-mini-dna-006-tumor"
    "scga-mini-dna-007-tumor"
    "scga-mini-dna-008-tumor"
    "scga-mini-dna-009-tumor"
    "scga-mini-dna-010-tumor"
    "scga-mini-dna-011-tumor"
    "scga-mini-dna-012-tumor"
    "scga-mini-dna-013-tumor"
    "scga-mini-dna-014-tumor"
    "scga-mini-dna-015-tumor"
    "scga-mini-dna-016-tumor"
)

for sample_id in "${sample_ids[@]}"; do
    echo $sample_id
    run "$sample_id"
done