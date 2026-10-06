import os
import pysam
import pytest
from test.data import get_data_path
from exactolib.constants import GeneAnnotationSource
from exactolib.main import remove_unspliced_rnas


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]
MIN_MAPPING_QUALITY = 30


def _remove_unspliced_rnas(sample_number: str, output_bam_file: str, output_bam_bai_file: str, write_output_bam_file: bool = True):
    bam_file = get_data_path(name='alignment/scga-mini-rna-%s-tumor_minimap2_sorted.bam' % sample_number)
    return remove_unspliced_rnas(
        bam_file=bam_file,
        bam_bai_file=bam_file + '.bai',
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        gene_types=['protein_coding'],
        gene_levels=[1,2],
        transcript_types=['protein_coding'],
        transcript_levels=[1,2],
        output_bam_file=output_bam_file,
        output_bam_bai_file=output_bam_bai_file,
        num_threads=2,
        min_mapping_quality=MIN_MAPPING_QUALITY,
        write_output_bam_file=write_output_bam_file
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_remove_unspliced_rnas_keeps_spliced_reads(sample_number, tmp_path):
    output_bam_file = str(tmp_path / 'spliced.bam')
    read_names_to_keep = _remove_unspliced_rnas(sample_number, output_bam_file, output_bam_file + '.bai')

    # A read is kept when a primary or supplementary record reaches the minimum mapping quality and
    # splices (no scga-mini read overlaps a single-exon transcript, the other way to be kept).
    # Nascent reads are partially spliced, so they stay too.
    spliced_read_names = set()
    with pysam.AlignmentFile(get_data_path(name='alignment/scga-mini-rna-%s-tumor_minimap2_sorted.bam' % sample_number)) as bam:
        records = [record for record in bam if not record.is_unmapped and not record.is_secondary]
    for record in records:
        if record.mapping_quality >= MIN_MAPPING_QUALITY and any(operation == pysam.CREF_SKIP for operation, _ in record.cigartuples):
            spliced_read_names.add(record.query_name)
    assert set(read_names_to_keep) == spliced_read_names
    assert 0 < len(spliced_read_names) < len(set(record.query_name for record in records))

    # The output BAM holds every primary and supplementary record of the kept reads, in input order.
    expected_records = [(record.query_name, record.reference_start) for record in records if record.query_name in spliced_read_names]
    with pysam.AlignmentFile(output_bam_file) as bam:
        assert [(record.query_name, record.reference_start) for record in bam] == expected_records
    assert os.path.exists(output_bam_file + '.bai')


def test_remove_unspliced_rnas_without_writing_output_bam_file(tmp_path):
    output_bam_file = str(tmp_path / 'spliced.bam')
    read_names_written = _remove_unspliced_rnas('008', str(tmp_path / 'written.bam'), str(tmp_path / 'written.bam.bai'))
    read_names_to_keep = _remove_unspliced_rnas('008', output_bam_file, output_bam_file + '.bai', write_output_bam_file=False)

    assert sorted(read_names_to_keep) == sorted(read_names_written)
    assert not os.path.exists(output_bam_file)
