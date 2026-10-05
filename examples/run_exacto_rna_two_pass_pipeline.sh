# Two-pass RNA pipeline: cluster raw reads, correct them against their cluster's
# variant list, realign, re-cluster the corrected reads, and build one consensus
# per transcript.

# Requires alongside exacto: minimap2, samtools.
set -e

EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

SAMPLE=scga-mini-rna-001-tumor
REFERENCE=${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz
GTF=${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz
BAM=${EXACTO_TEST_DATA}/alignment/${SAMPLE}_minimap2_sorted.bam
OUT=outputs/rna-two-pass/${SAMPLE}

mkdir -p ${OUT}/pass1/ ${OUT}/pass2/ ${OUT}/consensus/

# Pass 1: cluster the raw reads.
exacto cluster-rna-reads \
  --bam-file ${BAM} \
  --bai-file ${BAM}.bai \
  --reference-genome-fasta-file ${REFERENCE} \
  --reference-gene-annotation-file ${GTF} \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-dir ${OUT}/pass1/ \
  --output-prefix ${SAMPLE}-pass1 \
  --preset pb

# Correct the raw reads against their pass-1 cluster's variant list.
exacto correct-rna-reads \
  --bam-file ${BAM} \
  --clusters-tsv-file ${OUT}/pass1/${SAMPLE}-pass1_exacto_rna_clusters.tsv \
  --cluster-reference-transcripts-tsv-file ${OUT}/pass1/${SAMPLE}-pass1_exacto_rna_clusters_reference_transcripts.tsv \
  --cluster-splice-junctions-tsv-file ${OUT}/pass1/${SAMPLE}-pass1_exacto_rna_clusters_splice_junctions.tsv \
  --cluster-variants-tsv-file ${OUT}/pass1/${SAMPLE}-pass1_exacto_rna_clusters_variants_passed.tsv \
  --reference-gene-annotation-file ${GTF} \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-dir ${OUT}/pass1/ \
  --output-prefix ${SAMPLE}

# Realign the corrected reads.
minimap2 -ax splice:hq -uf --cs --eqx -Y -L --secondary=no \
  ${REFERENCE} \
  ${OUT}/pass1/${SAMPLE}_exacto_rna_corrected_reads.fastq.gz \
  | samtools sort -o ${OUT}/pass2/${SAMPLE}-corrected_minimap2_sorted.bam
samtools index ${OUT}/pass2/${SAMPLE}-corrected_minimap2_sorted.bam

# Pass 2: re-cluster the corrected reads.
exacto cluster-rna-reads \
  --bam-file ${OUT}/pass2/${SAMPLE}-corrected_minimap2_sorted.bam \
  --bai-file ${OUT}/pass2/${SAMPLE}-corrected_minimap2_sorted.bam.bai \
  --reference-genome-fasta-file ${REFERENCE} \
  --reference-gene-annotation-file ${GTF} \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-dir ${OUT}/pass2/ \
  --output-prefix ${SAMPLE}-pass2 \
  --preset corrected \
  --allowed-variants-tsv-file ${OUT}/pass1/${SAMPLE}-pass1_exacto_rna_clusters_variants_passed.tsv

# One consensus per pass-2 cluster, built from the corrected reads.
exacto determine-rna-consensus \
  --tsv-file ${OUT}/pass2/${SAMPLE}-pass2_exacto_rna_clusters.tsv \
  --fastq-file ${OUT}/pass1/${SAMPLE}_exacto_rna_corrected_reads.fastq.gz \
  --output-dir ${OUT}/consensus/ \
  --output-prefix ${SAMPLE}

# Align the consensus sequences. -un: consensus orientation is not guaranteed
# transcript-forward.
minimap2 -ax splice:hq -un --cs --eqx -Y -L --secondary=no \
  ${REFERENCE} \
  ${OUT}/consensus/${SAMPLE}_exacto_rna_consensus.fasta \
  | samtools sort -o ${OUT}/consensus/${SAMPLE}-consensus_minimap2_sorted.bam
samtools index ${OUT}/consensus/${SAMPLE}-consensus_minimap2_sorted.bam

# Next: run_exacto_stitch_reference_transcripts.sh on the consensus alignment,
# realign the stitched FASTA, then run_exacto_call_rna_transcript_vars.sh.
