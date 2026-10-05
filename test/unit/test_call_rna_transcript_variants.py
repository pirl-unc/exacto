import os
import pandas as pd
import pytest
from test.data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import identify_rna_transcript_variants


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]

# The calls put fusion and back-splice junctions up to 2 bp from the simulated breakpoint
# (007, 011, 015, 016); every other variant type is matched exactly.
JUNCTION_POSITION_TOLERANCE = 2
JUNCTION_VARIANT_TYPES = {'FUS', 'CIR', 'BND'}

# Ground truth rows (variant_id) that no call matches, pinned so that a fix shows up here:
# 112, the 012 duplication, is called as an insertion; 113 and 114, the 013 inversion breakends,
# come out elsewhere; 117, the second 015 fusion, is called from 6110799 to 6165489.
KNOWN_MISSES = {112, 113, 114, 117}

# The frames in the order identify_rna_transcript_variants returns them, by the suffix of the
# file each is written to.
OUTPUT_FILE_SUFFIXES = [
    'assembled_transcripts',
    'assembled_transcript_exons',
    'assembled_transcript_splice_junctions',
    'assembled_transcript_reference_transcript_matches',
    'assembled_transcript_filter_status',
    'assembled_transcript_model_alignments',
    'assembled_transcript_variants',
    'assembled_transcript_nmd_predictions'
]


def _identify_rna_transcript_variants(sample_number: str, output_dir: str = '', output_type: OutputType = OutputType.DATAFRAME, dna_variants_tsv_files=[]):
    # The stitched transcripts, as scripts/data/05_variant_calling/03_run_exacto_call_rna_transcript_vars.sh runs it
    return identify_rna_transcript_variants(
        bam_file=get_data_path(name='assembly/rna/stitched/scga-mini-rna-%s-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts_minimap2_sorted.bam' % sample_number),
        reference_genome_fasta_file=get_data_path(name='references/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir=output_dir,
        output_prefix='scga-mini-rna-%s-tumor' % sample_number,
        dna_variants_tsv_files=dna_variants_tsv_files,
        num_threads=4,
        output_type=output_type
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_identify_rna_transcript_variants_matches_ground_truth(sample_number):
    df_rna_variants = _identify_rna_transcript_variants(sample_number=sample_number)[6]
    df_ground_truth = pd.read_csv(get_data_path(name='simulation/ground_truth/scga-mini-rna-%s-tumor_ground_truth.tsv' % sample_number), sep='\t', keep_default_na=False)

    for _, row in df_ground_truth.iterrows():
        tolerance = JUNCTION_POSITION_TOLERANCE if row['variant_type'] in JUNCTION_VARIANT_TYPES else 0
        # Calls are written genomic-ascending; a truth row may name its breakpoints the other way round.
        df_matched = df_rna_variants[
            ((df_rna_variants['chromosome_1'] == row['chromosome_1']) &
             df_rna_variants['position_1'].between(row['position_1'] - tolerance, row['position_1'] + tolerance) &
             (df_rna_variants['operation_1'] == row['operation_1']) &
             (df_rna_variants['chromosome_2'] == row['chromosome_2']) &
             df_rna_variants['position_2'].between(row['position_2'] - tolerance, row['position_2'] + tolerance) &
             (df_rna_variants['operation_2'] == row['operation_2'])) |
            ((df_rna_variants['chromosome_1'] == row['chromosome_2']) &
             df_rna_variants['position_1'].between(row['position_2'] - tolerance, row['position_2'] + tolerance) &
             (df_rna_variants['operation_1'] == row['operation_2']) &
             (df_rna_variants['chromosome_2'] == row['chromosome_1']) &
             df_rna_variants['position_2'].between(row['position_1'] - tolerance, row['position_1'] + tolerance) &
             (df_rna_variants['operation_2'] == row['operation_1']))
        ]
        if row['variant_id'] in KNOWN_MISSES:
            assert len(df_matched) == 0, 'ground truth row %i is now called; take it out of KNOWN_MISSES' % row['variant_id']
        else:
            # One transcript per simulated mutant isoform carries the variant.
            assert len(df_matched) == 1, 'ground truth row %i' % row['variant_id']
            assert df_matched['variant_type'].iloc[0] == row['variant_type']


def test_identify_rna_transcript_variants_writes_every_table(tmp_path):
    # In file mode the tables are written, not returned; each file holds the frame of the same position.
    _identify_rna_transcript_variants(sample_number='001', output_dir=str(tmp_path), output_type=OutputType.FILE)
    frames = _identify_rna_transcript_variants(sample_number='001')

    assert len(frames) == len(OUTPUT_FILE_SUFFIXES)
    for df, suffix in zip(frames, OUTPUT_FILE_SUFFIXES):
        output_tsv_file = str(tmp_path / ('scga-mini-rna-001-tumor_exacto_%s.tsv' % suffix))
        assert os.path.exists(output_tsv_file), suffix
        df_file = pd.read_csv(output_tsv_file, sep='\t')
        assert list(df_file.columns) == list(df.columns), suffix
        assert len(df_file) == len(df), suffix


def test_identify_rna_transcript_variants_returns_the_two_tp53_transcripts():
    (df_assembled_transcripts,
     df_exons,
     df_splice_junctions,
     df_reference_transcript_matches,
     df_filter_status,
     df_model_alignments,
     df_rna_variants,
     df_nmd_predictions) = _identify_rna_transcript_variants(sample_number='001')

    # scga-mini-rna-001 simulates TP53 ENST00000269305.9 with and without the M246I SNV; each is one
    # stitched transcript with 11 exons, and only the mutant carries a variant.
    assert len(df_assembled_transcripts) == 2
    assert df_assembled_transcripts['num_exons'].tolist() == [11, 11]
    assert sorted(df_assembled_transcripts['is_variant'].tolist()) == [False, True]
    assert set(df_reference_transcript_matches['reference_transcript_id']) == {'ENST00000269305.9'}
    assert df_rna_variants[['variant_type', 'reference_gene_name']].values.tolist() == [['SNV', 'TP53']]
    assert not df_nmd_predictions['nmd_predicted'].any()


def test_identify_rna_transcript_variants_marks_variants_seen_in_dna():
    dna_variants_tsv_file = get_data_path(name='variant_calling/dna/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv')

    df_rna_variants = _identify_rna_transcript_variants(sample_number='001')[6]
    assert df_rna_variants['origin'].tolist() == ['']

    df_rna_variants = _identify_rna_transcript_variants(sample_number='001', dna_variants_tsv_files=[dna_variants_tsv_file])[6]
    assert df_rna_variants['origin'].tolist() == ['somatic']
