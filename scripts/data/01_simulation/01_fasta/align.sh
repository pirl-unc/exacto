nexus run --nf-workflow alignment_minimap2.nf \
  -c /Users/ajslee/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/ajslee/Documents/Research/projects/project_exacto/data/processed/work/alignment_minimap2/ \
  --samples_tsv_file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/scripts/data/01_simulation/01_fasta/samples_minimap2.tsv \
  --params_minimap2 '"-ax splice:hq -uf --cs --eqx -Y -L -k 19 --secondary=no --splice-flank=no"' \
  --reference_genome_fasta_file /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz \
  --output_dir /Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/simulation/ground_truth/
