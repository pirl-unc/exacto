EXACTO_TEST_DATA=${EXACTO_TEST_DATA:-../test/data}

mkdir -p outputs/build-genome-var-graph/

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/variant_calling/dna/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/references/hg38_chr17-18.fa.gz \
  --output-fasta-file outputs/build-genome-var-graph/scga-mini-dna-001-tumor_genome_var_graph.fasta \
  --sequence-prefix scga-mini-dna-001-tumor

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample_dna_variant_callset_1.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_1.fasta \
  --sequence-prefix sample_dna_variant_callset_1 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample_dna_variant_callset_2.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_2.fasta \
  --sequence-prefix sample_dna_variant_callset_2 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample_dna_variant_callset_3.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_3.fasta \
  --sequence-prefix sample_dna_variant_callset_3 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample_dna_variant_callset_4.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_4.fasta \
  --sequence-prefix sample_dna_variant_callset_4 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample_dna_variant_callset_5.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample2.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_5_1.fasta \
  --sequence-prefix sample_dna_variant_callset_5 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample_dna_variant_callset_5.tsv \
  --fasta-file ${EXACTO_TEST_DATA}/exacto/exacto-graph/sample2.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_5_2.fasta \
  --sequence-prefix sample_dna_variant_callset_5 \
  --remove-unknown-bases no
