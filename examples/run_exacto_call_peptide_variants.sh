#gtime -v \
#  exacto call-peptide-vars \
#    --fasta-file ../test/data/fasta/sample_peptide_sequences.fa.gz \
#    --rna-bam-file "" \
#    --rna-bam-bai-file "" \
#    --reference-fasta-file ../test/data/fasta/reference_peptide_sequences.fa.gz \
#    --translations-tsv-file ../test/data/tsv/translations/sample_translations.tsv.gz \
#    --rna-variants-tsv-file ../test/data/tsv/variants/sample_rna_variants.tsv.gz \
#    --dna-variants-tsv-file ../test/data/tsv/variants/sample_dna_variants.tsv.gz \
#    --excl-bed-file "" \
#    --output-tsv-file outputs/variant_peptide_sequences.tsv \
#    --output-fasta-file outputs/variant_peptide_sequences.fasta \
#    --min-reads 3 \
#    --k 8 \
#    --dna-variant-padding 1000 \
#    --num-threads 2 \
#    --gzip yes
