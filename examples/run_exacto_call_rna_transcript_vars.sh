EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/call-rna-transcript-vars/

exacto call-rna-transcript-vars \
  --bam-file ${EXACTO_TEST_DATA}/assembly/rna/stitched/scga-mini-rna-001-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts_minimap2_sorted.bam \
  --reference-genome-fasta-file ${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz \
  --reference-gene-annotation-file ${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-dir outputs/call-rna-transcript-vars/ \
  --output-prefix scga-mini-rna-001-tumor
