run() {
  local sample_id="$1"

  exacto correct-rna-reads \
    --bam-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/alignment/${sample_id}_minimap2_sorted.bam \
    --clusters-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_clusters.tsv \
    --cluster-reference-transcripts-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_clusters_reference_transcripts.tsv \
    --cluster-splice-junctions-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_clusters_splice_junctions.tsv \
    --cluster-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/${sample_id}_pass1_exacto_rna_clusters_variants_passed.tsv \
    --output-dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass1/ \
    --reference-gene-annotation-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/gencode.v41.annotation.chr17-18.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly hg38 \
    --reference-gene-annotation-version v41 \
    --output-prefix ${sample_id}_pass1 \
    --num-threads 16
}

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