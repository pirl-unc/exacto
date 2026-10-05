run() {
  local sample_id="$1"
  local dna_sample_id="${sample_id/-rna-/-dna-}"

  exacto integrate-vars \
    --dna-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/dna/${dna_sample_id}_exacto_somatic_dna_variants.tsv \
    --rna-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/rna/${sample_id}_exacto_assembled_transcript_variants.tsv \
    --reference-gene-annotation-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/gencode.v41.annotation.chr17-18.gtf.gz \
    --reference-gene-annotation-source gencode \
    --reference-gene-annotation-assembly hg38 \
    --reference-gene-annotation-version v41 \
    --output-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/integration/${sample_id}_exacto_integrated_variants.tsv \
    --num-threads 16
}

mkdir -p /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/integration/

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
