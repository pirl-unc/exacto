EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/correct-rna-reads/

exacto correct-rna-reads \
  --bam-file ${EXACTO_TEST_DATA}/alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam \
  --clusters-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass1/scga-mini-rna-001-tumor_pass1_exacto_rna_clusters.tsv \
  --cluster-reference-transcripts-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass1/scga-mini-rna-001-tumor_pass1_exacto_rna_clusters_reference_transcripts.tsv \
  --cluster-splice-junctions-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass1/scga-mini-rna-001-tumor_pass1_exacto_rna_clusters_splice_junctions.tsv \
  --cluster-variants-tsv-file ${EXACTO_TEST_DATA}/assembly/rna/pass1/scga-mini-rna-001-tumor_pass1_exacto_rna_clusters_variants_passed.tsv \
  --reference-gene-annotation-file ${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-dir outputs/correct-rna-reads/ \
  --output-prefix scga-mini-rna-001-tumor_pass1
