EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/integrate-vars/

exacto integrate-vars \
  --dna-variants-tsv-file ${EXACTO_TEST_DATA}/variant_calling/dna/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv \
  --rna-variants-tsv-file ${EXACTO_TEST_DATA}/variant_calling/rna/scga-mini-rna-001-tumor_exacto_assembled_transcript_variants.tsv \
  --reference-gene-annotation-file ${EXACTO_TEST_DATA}/references/gencode.v41.annotation.chr17-18.gtf.gz \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v41 \
  --output-tsv-file outputs/integrate-vars/scga-mini-rna-001-tumor_exacto_integrated_variants.tsv
