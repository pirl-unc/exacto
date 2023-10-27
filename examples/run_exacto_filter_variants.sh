exacto filter-variants \
  --tsv-file ../test/data/hg002_all_variants.tsv \
  --case-sample-id HG002 \
  --output-filtered-tsv-file ../test/data/hg002_all_varaints_refined_hq.tsv \
  --output-rejected-tsv-file ../test/data/hg002_all_variants_refined_lq.tsv \
  --filter 'case all filter == PASS' \
  --filter 'case all precise == true' \
  --filter 'case all alternate_allele_read_count >= 3' \
  --filter 'case all chromosome_1 in ["chr1","chr2","chr3"]' \
  --filter 'case all chromosome_2 in ["chr1","chr2","chr3"]' \
  --excluded-regions-tsv-file ../test/data/hg38_ucsc_gap_table.tsv \
  --excluded-regions-padding 100000 \
  --excluded-variants-tsv-file ../test/data/hg002_pbsv.tsv \
  --excluded-variants-padding 1000 \
  --num-threads 4

exacto filter-variants \
  --tsv-file ../test/data/hg002_sniffles2.tsv \
  --case-sample-id HG002 \
  --output-filtered-tsv-file hg002_sniffles2_refined_hq.tsv \
  --output-rejected-tsv-file hg002_sniffles2_refined_lq.tsv \
  --filter 'case all filter == PASS' \
  --filter 'case all precise == true' \
  --filter 'case all alternate_allele_read_count >= 3' \
  --filter 'case all chromosome_1 in ["chr1","chr2","chr3"]' \
  --excluded-regions-tsv-file ../test/data/hg38_ucsc_gap_table.tsv \
  --excluded-regions-padding 100000 \
  --excluded-variants-tsv-file ../test/data/hg002_pbsv.tsv \
  --excluded-variants-padding 1000 \
  --num-threads 4

exacto filter-variants \
  --tsv-file ../test/data/hg002_pbsv.tsv \
  --case-sample-id HG002 \
  --output-filtered-tsv-file hg002_pbsv_refined_hq.tsv \
  --output-rejected-tsv-file hg002_pbsv_refined_lq.tsv \
  --filter 'case all filter == PASS' \
  --filter 'case all alternate_allele_read_count >= 3' \
  --filter 'case all chromosome_1 in ["chr1","chr2","chr3"]' \
  --excluded-regions-tsv-file ../test/data/hg38_ucsc_gap_table.tsv \
  --excluded-regions-padding 100000 \
  --excluded-variants-tsv-file ../test/data/hg002_sniffles2.tsv \
  --excluded-variants-padding 1000 \
  --num-threads 4
