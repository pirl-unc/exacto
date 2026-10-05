run() {
  local sample_id="$1"
  local dna_sample_id="${sample_id/-rna-/-dna-}"

  exacto translate-transcripts \
    --assembled-transcript-model-alignments-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/rna/${sample_id}_exacto_assembled_transcript_model_alignments.tsv \
    --assembled-transcript-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/rna/${sample_id}_exacto_assembled_transcript_variants.tsv \
    --dna-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/variant_calling/dna/${dna_sample_id}_exacto_somatic_dna_variants.tsv \
    --integrated-variants-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/integration/${sample_id}_exacto_integrated_variants.tsv \
    --rna-consensus-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/pass2/${sample_id}_pass2_exacto_rna_consensus.tsv \
    --stitched-transcripts-tsv-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/assembly/rna/stitched/${sample_id}_pass2_exacto_rna_consensus_exacto_stitched_transcripts.tsv \
    --strategy longest_orf \
    --output-dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/translation/ \
    --output-prefix ${sample_id} \
    --num-threads 16
}

mkdir -p /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/translation/

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
