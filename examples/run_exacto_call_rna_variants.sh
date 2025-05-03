gtime -v \
  exacto call-rna-vars \
    --bam-file ../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam \
    --bam-bai-file ../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai \
    --reference-genome-fasta-file ../test/data/fasta/hg38_chr17-18.fa.gz \
    --gene-annotation-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
    --gene-annotation-source gencode \
    --output-dir outputs/ \
    --output-prefix rna-100-tumor_minimap2_mdtagged_sorted \
    --reference-transcript-scoring-method cosine_similarity \
    --reference-transcript-selection-strategy top_k \
    --reference-transcript-top-k 3
