run() {
  local sample_id="$1"

  exacto cluster-rna-reads \
    --bam-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_corrected_minimap2_sorted.bam \
    --bai-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_corrected_minimap2_sorted.bam.bai \
    --reference-genome-fasta-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz \
    --reference-gene-annotation-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/gencode.v41.annotation.chr17-18.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly hg38 \
    --reference-gene-annotation-version v41 \
    --output-dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass2/ \
    --output-prefix ${sample_id}_pass2 \
    --preset corrected \
    --allowed-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_clusters_variants_passed.tsv \
    --num-threads 16
}

mkdir -p /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass2/

sample_ids=(
    "scga-mini-rna-001-tumor"
    "scga-mini-rna-002-tumor"
    "scga-mini-rna-003-tumor"
    "scga-mini-rna-004-tumor"
    "scga-mini-rna-005-tumor"
    "scga-mini-rna-006-tumor"
    "scga-mini-rna-007-tumor"
    "scga-mini-rna-008-tumor"
    "scga-mini-rna-009-tumor"
    "scga-mini-rna-010-tumor"
    "scga-mini-rna-011-tumor"
    "scga-mini-rna-012-tumor"
    "scga-mini-rna-013-tumor"
    "scga-mini-rna-014-tumor"
    "scga-mini-rna-015-tumor"
    "scga-mini-rna-016-tumor"
)

for sample_id in "${sample_ids[@]}"; do
    echo $sample_id
    run "$sample_id"
done