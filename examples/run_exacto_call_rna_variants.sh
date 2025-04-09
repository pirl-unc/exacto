gtime -v \
  exacto call-rna-vars \
    --bam-file ../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam \
    --bam-bai-file ../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai \
    --reference-genome-fasta-file ../test/data/fasta/hg38_chr17-18.fa.gz \
    --gene-annotation-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
    --gene-annotation-source gencode \
    --output-ref-transcript-matches-tsv-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_reference_transcript_matches.tsv.gz \
    --output-exons-tsv-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_exons.tsv.gz \
    --output-sj-tsv-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_splice_junctions.tsv.gz \
    --output-variants-tsv-file outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_variants.tsv.gz
