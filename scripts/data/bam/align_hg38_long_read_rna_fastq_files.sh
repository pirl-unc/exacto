nexus run --nf-workflow long_read_alignment_minimap2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_exacto/data/processed/work/long_read_alignment_minimap2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_exacto/exacto/scripts/data/bam/samples_long_read_rna_fastq_files.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M.fa.gz \
  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M.fa.gz.fai \
  --params_minimap2 '"-ax splice:hq -uf --cs --eqx -Y -L"' \
  --output_dir /Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/bam/

cp /Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/bam/*.bam /Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/bam/
cp /Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/bam/*.bam.bai /Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/bam/
