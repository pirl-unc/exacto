nexus run --nf-workflow utilities_pbccs.nf \
  -c /Users/ajslee/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/ajslee/Documents/Research/projects/project_exacto/data/processed/work/utilities_pbccs/ \
  --nextflow /opt/miniconda3/envs/exacto/bin/nextflow \
  --samples_tsv_file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/scripts/data/01_simulation/03_fastq/samples/samples_pbccs-rna.tsv \
  --output_dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/simulation/fastq/pbccs/