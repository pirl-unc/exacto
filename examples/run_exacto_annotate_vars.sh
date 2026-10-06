EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/annotate-vars/

exacto annotate-vars \
  --tsv-file ${EXACTO_TEST_DATA}/variant_calling/dna/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv \
  --reference-gene-annotation-file ${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-tsv-file outputs/annotate-vars/scga-mini-dna-001-tumor_exacto_somatic_dna_variants_annotated.tsv
