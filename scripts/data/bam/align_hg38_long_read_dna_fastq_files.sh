nexus run --nf-workflow alignment_minimap2.nf \
  -c /Users/ajslee/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/ajslee/Documents/Research/projects/project_exacto/data/processed/work/alignment_minimap2/ \
  --samples_tsv_file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/scripts/data/bam/samples_long_read_dna_fastq_files.tsv \
  --params_minimap2 '"-ax map-hifi --cs --eqx -Y -L"' \
  --reference_genome_fasta_file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/fasta/hg38_chr17-18.fa.gz \
  --output_dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/bam/

cp /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/bam/*.bam /Users/ajslee/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/bam/
cp /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/bam/*.bam.bai /Users/ajslee/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/bam/