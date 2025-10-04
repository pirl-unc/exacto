import pandas as pd
from .data import get_data_path
from exactolib.constants import OutputType
from exactolib.main import identify_dna_variants,identify_case_specific_dna_variants


def test_call_dna_variants_1():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-001-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-001-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] == row['position_1']) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] == row['position_2']) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_dna_variants_2():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-002-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-002-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] == row['position_1']) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] == row['position_2']) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_dna_variants_3():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-003-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-003-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-003-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] == row['position_1']) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] == row['position_2']) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_dna_variants_4():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-004-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )
    assert(len(df_variants) == 2)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-004-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 1000) &
            (df_variants['position_1'] >= row['position_1'] - 1000) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 1000) &
            (df_variants['position_2'] >= row['position_2'] - 1000) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_dna_variants_5():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-005-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-005-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-005-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_dna_variants_6():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-006-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-006-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_dna_variants_7():
    df_variants = identify_dna_variants(
        bam_file=get_data_path(name='bam/dna-007-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/dna-007-tumor_minimap2_mdtagged_sorted.bam.bai'),
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-007-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_1():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-001-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-001-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-001-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-001-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_2():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-002-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-002-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-002-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-002-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_3():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-003-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-003-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-003-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-003-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-003-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_4():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-004-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-004-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-004-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 2)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-004-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 1000) &
            (df_variants['position_1'] >= row['position_1'] - 1000) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 1000) &
            (df_variants['position_2'] >= row['position_2'] - 1000) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_5():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-005-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-005-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-005-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-005-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-005-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_6():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-006-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-006-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-006-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-006-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_case_specific_dna_variants_7():
    df_variants = identify_case_specific_dna_variants(
        case_bam_file=get_data_path(name='bam/dna-007-tumor_minimap2_mdtagged_sorted.bam'),
        case_bam_bai_file=get_data_path(name='bam/dna-007-tumor_minimap2_mdtagged_sorted.bam.bai'),
        control_bam_files=[get_data_path(name='bam/dna-007-normal_minimap2_mdtagged_sorted.bam')],
        control_bam_bai_files=[get_data_path(name='bam/dna-007-normal_minimap2_mdtagged_sorted.bam.bai')],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/dna-007-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1
