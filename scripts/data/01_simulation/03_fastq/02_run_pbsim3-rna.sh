nexus run --nf-workflow sequencing_simulation_pbsim3-rna.nf \
  -c /Users/ajslee/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/ajslee/Documents/Research/projects/project_exacto/data/processed/work/sequencing_simulation_pbsim3-rna/ \
  --nextflow /opt/miniconda3/envs/exacto/bin/nextflow \
  --samples_tsv_file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/scripts/data/01_simulation/03_fastq/samples/samples_pbsim3-rna.tsv \
  --output_dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/simulation/fastq/pbsim3/ \
  --pbsim3_model_file /Users/ajslee/Documents/Research/projects/seqdata/tool-resources/pbsim3/pbsim3-3.0.5/data/ERRHMM-SEQUEL.model \
  --params_pbsim3_mode '"errhmm"' \
  --params_pbsim3 '"--length-mean 16379 --length-sd 2043 --pass-num 10 --strategy trans"'
