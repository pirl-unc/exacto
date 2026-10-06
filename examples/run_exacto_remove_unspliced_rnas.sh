EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/remove-unspliced-rnas/

exacto remove-unspliced-rnas \
  --bam-file ${EXACTO_TEST_DATA}/alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam \
  --bam-bai-file ${EXACTO_TEST_DATA}/alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai \
  --reference-gene-annotation-file ${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-bam-file outputs/remove-unspliced-rnas/scga-mini-rna-001-tumor_minimap2_sorted.spliced.bam \
  --output-bam-bai-file outputs/remove-unspliced-rnas/scga-mini-rna-001-tumor_minimap2_sorted.spliced.bam.bai
