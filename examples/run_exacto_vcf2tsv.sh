exacto vcf2tsv \
  --vcf-file ../test/data/hg002_cutesv.vcf \
  --variant-calling-method cutesv \
  --sequencing-platform pacbio \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_cutesv.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_deepvariant.vcf \
  --variant-calling-method deepvariant \
  --sequencing-platform pacbio \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_deepvariant.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_pbsv.vcf \
  --variant-calling-method pbsv \
  --sequencing-platform pacbio \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_pbsv.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_sniffles2.vcf \
  --variant-calling-method sniffles2 \
  --sequencing-platform pacbio \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_sniffles2.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_svim.vcf \
  --variant-calling-method svim \
  --sequencing-platform pacbio \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_svim.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_delly2.vcf \
  --variant-calling-method delly2-somatic \
  --sequencing-platform ilmn \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_hg001_delly2.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_gatk4_mutect2.vcf \
  --variant-calling-method gatk4-mutect2 \
  --sequencing-platform ilmn \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_hg001_gatk4_mutect2.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_lumpy.vcf \
  --variant-calling-method lumpy-somatic \
  --sequencing-platform ilmn \
  --source-id hg002 \
  --output-tsv-file ../test/data/hg002_hg001_lumpy.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_strelka2_indels.vcf \
  --variant-calling-method strelka2-somatic \
  --sequencing-platform ilmn \
  --source-id hg002 \
  --case-id hg002 \
  --control-id hg001 \
  --output-tsv-file ../test/data/hg002_hg001_strelka2_indels.tsv

exacto vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_strelka2_snvs.vcf \
  --variant-calling-method strelka2-somatic \
  --sequencing-platform ilmn \
  --source-id hg002 \
  --case-id hg002 \
  --control-id hg001 \
  --output-tsv-file ../test/data/hg002_hg001_strelka2_snvs.tsv