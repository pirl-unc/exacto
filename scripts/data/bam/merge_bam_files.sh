samtools merge \
  -o ../../../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted_dna-001-tumor_minimap2_mdtagged_sorted.bam \
  ../../../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam \
  ../../../test/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam

samtools index ../../../test/data/bam/rna-100-tumor_minimap2_mdtagged_sorted_dna-001-tumor_minimap2_mdtagged_sorted.bam

gzcat ../../../test/data/fastq/rna-100-tumor_long-read.fastq.gz \
     ../../../test/data/fastq/dna-001-tumor_long-read.fastq.gz \
| awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' \
> ../../../test/data/fasta/rna-100-tumor_long-read_dna-001-tumor_long-read.fasta