import collections
import os
import pytest
from test.data import get_data_path, get_simulated_read_origins
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.default import CLUSTER_RNA_READS_PRESETS
from exactolib.main import cluster_rna_reads


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]

OUTPUT_FILE_SUFFIXES = [
    'rna_clusters',
    'rna_clusters_summary',
    'rna_clusters_reference_transcripts',
    'rna_clusters_splice_junctions',
    'rna_clusters_variants_passed',
    'rna_clusters_variants_failed',
    'rna_clusters_template_switch'
]


def _cluster_rna_reads(bam_file: str, preset: str, output_dir: str, output_prefix: str = 'scga-mini', output_type: OutputType = OutputType.DATAFRAME, **options):
    # The tables are written to output_dir in dataframe mode too, so it is always a temporary directory.
    return cluster_rna_reads(
        bam_file=bam_file,
        bai_file=bam_file + '.bai',
        reference_genome_fasta_file=get_data_path(name='references/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir=output_dir,
        output_prefix=output_prefix,
        num_threads=4,
        output_type=output_type,
        **CLUSTER_RNA_READS_PRESETS[preset],
        **options
    )


def _check_clusters_match_simulated_transcripts(sample_number: str, df_clusters):
    # Every numbered row of the .transcript file (200 reads each) is the majority of exactly one
    # cluster; the partially spliced rows (nascent RNA, 10 reads each) are the majority of none.
    read_origins = get_simulated_read_origins(sample_id='scga-mini-rna-%s-tumor' % sample_number)
    with open(get_data_path(name='simulation/transcript/scga-mini-rna-%s-tumor.transcript' % sample_number)) as handle:
        transcript_names = [line.split('\t')[0] for line in handle if line.strip()]
    mature_transcript_names = [name for name in transcript_names if 'partially-spliced' not in name]

    majority_transcript_names = []
    for cluster_id in df_clusters['cluster_id'].unique().to_list():
        read_names = df_clusters.filter(df_clusters['cluster_id'] == cluster_id)['read_name'].to_list()
        counts = collections.Counter(read_origins[read_name] for read_name in read_names)
        majority_transcript_names.append(counts.most_common(1)[0][0])
    assert sorted(majority_transcript_names) == sorted(mature_transcript_names)


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_cluster_rna_reads_pass_1_finds_simulated_transcripts(sample_number, tmp_path):
    # Pass 1, as scripts/data/04_assembly/rna/01_run_exacto_cluster_rna_reads_pass1.sh runs it (--preset pb)
    frames = _cluster_rna_reads(
        bam_file=get_data_path(name='alignment/scga-mini-rna-%s-tumor_minimap2_sorted.bam' % sample_number),
        preset='pb',
        output_dir=str(tmp_path)
    )

    _check_clusters_match_simulated_transcripts(sample_number=sample_number, df_clusters=frames[0])


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_cluster_rna_reads_pass_2_finds_simulated_transcripts(sample_number, tmp_path):
    # Pass 2 on the corrected reads, as scripts/data/04_assembly/rna/04_run_exacto_cluster_rna_reads_pass2.sh
    # runs it (--preset corrected, pass-1 passed variants as the allow list)
    frames = _cluster_rna_reads(
        bam_file=get_data_path(name='assembly/rna/pass1/scga-mini-rna-%s-tumor_pass1_exacto_rna_corrected_minimap2_sorted.bam' % sample_number),
        preset='corrected',
        output_dir=str(tmp_path),
        allowed_variants_tsv_file=get_data_path(name='assembly/rna/pass1/scga-mini-rna-%s-tumor_pass1_exacto_rna_clusters_variants_passed.tsv' % sample_number)
    )

    _check_clusters_match_simulated_transcripts(sample_number=sample_number, df_clusters=frames[0])


def test_cluster_rna_reads_finds_tp53_m246i(tmp_path):
    (df_clusters,
     df_summary,
     df_reference_transcripts,
     df_splice_junctions,
     df_variants_passed,
     df_variants_failed,
     df_template_switch) = _cluster_rna_reads(
        bam_file=get_data_path(name='alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam'),
        preset='pb',
        output_dir=str(tmp_path)
    )

    # Two TP53 ENST00000269305.9 clusters; only the mutant one passes a variant, the c.738G>A SNV.
    assert df_summary['cluster_id'].to_list() == df_clusters['cluster_id'].unique(maintain_order=True).to_list()
    assert df_summary['num_reads'].sum() == len(df_clusters)
    assert set(df_reference_transcripts['reference_transcript_id'].to_list()) == {'ENST00000269305.9'}
    assert df_variants_passed.select(['chromosome_1', 'position_1', 'position_2', 'variant_type']).rows() == [('chr17', 7674224, 7674226, 'SNV')]
    # Every failed variant is reported with the reason it failed.
    assert df_variants_failed['failure_reason'].null_count() == 0


def test_cluster_rna_reads_writes_every_table(tmp_path):
    # In file mode the tables are written, not returned; each file holds the frame of the same position.
    bam_file = get_data_path(name='alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam')
    _cluster_rna_reads(bam_file=bam_file, preset='pb', output_dir=str(tmp_path), output_prefix='scga-mini-rna-001-tumor_pass1', output_type=OutputType.FILE)
    frames = _cluster_rna_reads(bam_file=bam_file, preset='pb', output_dir=str(tmp_path), output_prefix='dataframe')

    assert len(frames) == len(OUTPUT_FILE_SUFFIXES)
    for df, suffix in zip(frames, OUTPUT_FILE_SUFFIXES):
        output_tsv_file = str(tmp_path / ('scga-mini-rna-001-tumor_pass1_exacto_%s.tsv' % suffix))
        assert os.path.exists(output_tsv_file), suffix
        with open(output_tsv_file) as handle:
            lines = handle.read().splitlines()
        # Same columns; the splice junctions file orders them chromosome_1, chromosome_2, position_1, ...
        # while its frame keeps each breakpoint's chromosome, position and strand together.
        assert sorted(lines[0].split('\t')) == sorted(df.columns), suffix
        assert len(lines) - 1 == len(df), suffix
