mkdir -p outputs/build-genome-var-graph/

exacto build-genome-var-graph \
  --variants-tsv-file ../test/data/tsv/variant_callset/sample_dna_variant_callset_1.tsv \
  --fasta-file ../test/data/fasta/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_1.fasta \
  --sequence-prefix sample_dna_variant_callset_1 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ../test/data/tsv/variant_callset/sample_dna_variant_callset_2.tsv \
  --fasta-file ../test/data/fasta/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_2.fasta \
  --sequence-prefix sample_dna_variant_callset_2 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ../test/data/tsv/variant_callset/sample_dna_variant_callset_3.tsv \
  --fasta-file ../test/data/fasta/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_3.fasta \
  --sequence-prefix sample_dna_variant_callset_3 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ../test/data/tsv/variant_callset/sample_dna_variant_callset_4.tsv \
  --fasta-file ../test/data/fasta/sample.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_4.fasta \
  --sequence-prefix sample_dna_variant_callset_4 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ../test/data/tsv/variant_callset/sample_dna_variant_callset_4.tsv \
  --fasta-file ../test/data/fasta/sample2.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_4_1.fasta \
  --sequence-prefix sample_dna_variant_callset_4 \
  --remove-unknown-bases yes

exacto build-genome-var-graph \
  --variants-tsv-file ../test/data/tsv/variant_callset/sample_dna_variant_callset_4.tsv \
  --fasta-file ../test/data/fasta/sample2.fa \
  --output-fasta-file outputs/build-genome-var-graph/sample_dna_variant_callset_4_2.fasta \
  --sequence-prefix sample_dna_variant_callset_4 \
  --remove-unknown-bases no

