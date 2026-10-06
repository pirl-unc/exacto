nexus run --nf-workflow haplotagging_long-read-dna.nf \
  -c /Users/ajslee/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/ajslee/Documents/Research/projects/project_exacto/data/processed/work/haplotagging_long-read-dna/ \
  --nextflow /opt/miniconda3/envs/exacto/bin/nextflow \
  -params-file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/scripts/data/03_haplotagging/dna/params.yaml
