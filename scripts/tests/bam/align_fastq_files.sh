#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -n 16 # Number of cores
#SBATCH --mem=64G
#SBATCH --job-name=align-fastq-files
#SBATCH -o align_fastq_files_slurm.out


HG38_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa


align_dna() {
  local fastq_file=$1
  local sam_file=$2
  local bam_file=$3
  local sample_id=$4

  minimap2 \
    -ax map-hifi --cs --eqx -Y -L \
    -t 16 \
    -R "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:$sample_id\\tLB:$sample_id\\tPU:$sample_id" \
    $HG38_FASTA_FILE \
    $fastq_file > $sam_file
  samtools sort $sam_file -o $bam_file
  samtools index -b $bam_file
}

align_rna() {
  local fastq_file=$1
  local sam_file=$2
  local bam_file=$3
  local sample_id=$4

  minimap2 \
    -ax splice:hq -uf --cs --eqx -Y -L \
    -t 16 \
    -R "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:$sample_id\\tLB:$sample_id\\tPU:$sample_id" \
    $HG38_FASTA_FILE \
    $fastq_file > $sam_file
  samtools sort $sam_file -o $bam_file
  samtools index -b $bam_file
}

align_dna \
  ../../../test/data/fastq/hg38_tumor_long_read_dna_1.fastq.gz \
  ../../../test/data/bam/hg38_tumor_long_read_dna_1.sam \
  ../../../test/data/bam/hg38_tumor_long_read_dna_1.bam \
  hg38_tumor_long_read_dna_1

align_dna \
  ../../../test/data/fastq/hg38_tumor_long_read_dna_2.fastq.gz \
  ../../../test/data/bam/hg38_tumor_long_read_dna_2.sam \
  ../../../test/data/bam/hg38_tumor_long_read_dna_2.bam \
  hg38_tumor_long_read_dna_2

align_dna \
  ../../../test/data/fastq/hg38_normal_long_read_dna_2.fastq.gz \
  ../../../test/data/bam/hg38_normal_long_read_dna_2.sam \
  ../../../test/data/bam/hg38_normal_long_read_dna_2.bam \
  hg38_normal_long_read_dna_2

align_rna \
  ../../../test/data/fastq/hg38_tumor_long_read_rna_1.fastq.gz \
  ../../../test/data/bam/hg38_tumor_long_read_rna_1.sam \
  ../../../test/data/bam/hg38_tumor_long_read_rna_1.bam \
  hg38_tumor_long_read_rna_1