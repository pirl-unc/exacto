EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/cluster-rna-reads/

exacto cluster-rna-reads \
  --bam-file ${EXACTO_TEST_DATA}/alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam \
  --bai-file ${EXACTO_TEST_DATA}/alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai \
  --reference-genome-fasta-file ${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz \
  --reference-gene-annotation-file ${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-dir outputs/cluster-rna-reads/ \
  --output-prefix scga-mini-rna-001-tumor_pass1 \
  --preset pb
