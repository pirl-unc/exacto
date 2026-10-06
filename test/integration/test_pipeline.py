import pandas as pd
import pysam
import pytest
from test.data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType, TranslationStrategy
from exactolib.default import CALL_DNA_VARS_PRESETS, CLUSTER_RNA_READS_PRESETS
from exactolib.main import (
    cluster_rna_reads,
    correct_rna_reads,
    determine_rna_consensus,
    identify_peptide_variants,
    identify_rna_transcript_variants,
    identify_somatic_dna_variants,
    integrate_variants,
    quantify_rna_abundances,
    stitch_reference_transcripts,
    translate_transcripts
)


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]

FASTA_FILE = get_data_path(name='references/hg38_chr17-18.fa.gz')
GTF_FILE = get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz')
GENE_ANNOTATION = dict(
    reference_gene_annotation_file=GTF_FILE,
    reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
    reference_gene_annotation_assembly='hg38',
    reference_gene_annotation_version='v41'
)


def _check_alignment_holds(bam_file: str, sequences: dict):
    """
    The pipeline realigns with minimap2 between some stages; the test reads the committed alignment
    instead (scripts/data/04_assembly/rna/0{3,6,8}_run_minimap2_*.sh). That stands in for realigning
    this run's output only if its primary records are exactly these sequences.
    """
    aligned_sequences = {}
    with pysam.AlignmentFile(bam_file) as bam:
        for record in bam:
            if record.is_secondary or record.is_supplementary or record.is_unmapped:
                continue
            sequence = record.get_forward_sequence()
            aligned_sequences[record.query_name] = sequence
    assert {name: sequence.upper() for name, sequence in aligned_sequences.items()} == \
           {name: sequence.upper() for name, sequence in sequences.items() if name in aligned_sequences}
    assert len(aligned_sequences) > 0


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_pipeline(sample_number, tmp_path):
    dna_sample_id = 'scga-mini-dna-%s-tumor' % sample_number
    rna_sample_id = 'scga-mini-rna-%s-tumor' % sample_number
    output_dir = str(tmp_path)

    # Step 1. Call somatic DNA variants (scripts/data/05_variant_calling/02)
    somatic_dna_variants_tsv_file = str(tmp_path / ('%s_exacto_somatic_dna_variants.tsv' % dna_sample_id))
    identify_somatic_dna_variants(
        bam_file=get_data_path(name='alignment/%s_minimap2_sorted.bam' % dna_sample_id),
        bam_bai_file=get_data_path(name='alignment/%s_minimap2_sorted.bam.bai' % dna_sample_id),
        control_bam_files=[get_data_path(name='alignment/scga-mini-dna-%s-normal_minimap2_sorted.bam' % sample_number)],
        control_bam_bai_files=[get_data_path(name='alignment/scga-mini-dna-%s-normal_minimap2_sorted.bam.bai' % sample_number)],
        fasta_file=FASTA_FILE,
        output_tsv_file=somatic_dna_variants_tsv_file,
        regions=[],
        num_threads=4,
        output_type=OutputType.FILE,
        **CALL_DNA_VARS_PRESETS['pb']
    )

    # Step 2. Cluster the RNA reads, pass 1 (scripts/data/04_assembly/rna/01)
    rna_bam_file = get_data_path(name='alignment/%s_minimap2_sorted.bam' % rna_sample_id)
    pass1_prefix = '%s_pass1' % rna_sample_id
    cluster_rna_reads(
        bam_file=rna_bam_file,
        bai_file=rna_bam_file + '.bai',
        reference_genome_fasta_file=FASTA_FILE,
        output_dir=output_dir,
        output_prefix=pass1_prefix,
        num_threads=4,
        output_type=OutputType.FILE,
        **GENE_ANNOTATION,
        **CLUSTER_RNA_READS_PRESETS['pb']
    )
    pass1_clusters = str(tmp_path / ('%s_exacto_rna_clusters' % pass1_prefix))

    # Step 3. Correct the clustered reads (scripts/data/04_assembly/rna/02)
    corrected_reads_fastq_file = correct_rna_reads(
        bam_file=rna_bam_file,
        clusters_tsv_file=pass1_clusters + '.tsv',
        cluster_reference_transcripts_tsv_file=pass1_clusters + '_reference_transcripts.tsv',
        cluster_splice_junctions_tsv_file=pass1_clusters + '_splice_junctions.tsv',
        cluster_variants_tsv_file=pass1_clusters + '_variants_passed.tsv',
        output_dir=output_dir,
        output_prefix=pass1_prefix,
        num_threads=4,
        **GENE_ANNOTATION
    )

    # Step 4. Realign the corrected reads (minimap2: the committed alignment)
    corrected_bam_file = get_data_path(name='assembly/rna/pass1/%s_exacto_rna_corrected_minimap2_sorted.bam' % pass1_prefix)
    with pysam.FastxFile(corrected_reads_fastq_file) as fastq:
        _check_alignment_holds(bam_file=corrected_bam_file, sequences={record.name: record.sequence for record in fastq})

    # Step 5. Cluster the corrected reads, pass 2 (scripts/data/04_assembly/rna/04)
    pass2_prefix = '%s_pass2' % rna_sample_id
    cluster_rna_reads(
        bam_file=corrected_bam_file,
        bai_file=corrected_bam_file + '.bai',
        reference_genome_fasta_file=FASTA_FILE,
        output_dir=output_dir,
        output_prefix=pass2_prefix,
        allowed_variants_tsv_file=pass1_clusters + '_variants_passed.tsv',
        num_threads=4,
        output_type=OutputType.FILE,
        **GENE_ANNOTATION,
        **CLUSTER_RNA_READS_PRESETS['corrected']
    )
    pass2_clusters = str(tmp_path / ('%s_exacto_rna_clusters' % pass2_prefix))

    # Step 6. Build one consensus per cluster (scripts/data/04_assembly/rna/05)
    determine_rna_consensus(
        tsv_file=pass2_clusters + '.tsv',
        fastq_file=corrected_reads_fastq_file,
        output_dir=output_dir,
        output_prefix=pass2_prefix,
        num_threads=2,
        output_type=OutputType.FILE
    )
    rna_consensus_tsv_file = str(tmp_path / ('%s_exacto_rna_consensus.tsv' % pass2_prefix))

    # Step 7. Realign the consensus sequences (minimap2: the committed alignment)
    consensus_bam_file = get_data_path(name='assembly/rna/pass2/%s_exacto_rna_consensus_minimap2_sorted.bam' % pass2_prefix)
    df_consensus = pd.read_csv(rna_consensus_tsv_file, sep='\t')
    _check_alignment_holds(bam_file=consensus_bam_file, sequences=dict(zip(df_consensus['cluster_id'].astype(str), df_consensus['consensus_sequence'])))

    # Step 8. Stitch reference sequence onto degraded ends (scripts/data/04_assembly/rna/07)
    stitch_prefix = '%s_exacto_rna_consensus' % pass2_prefix
    stitch_reference_transcripts(
        bam_file=consensus_bam_file,
        reference_genome_fasta_file=FASTA_FILE,
        output_dir=output_dir,
        output_prefix=stitch_prefix,
        num_threads=2,
        output_type=OutputType.FILE,
        **GENE_ANNOTATION
    )
    stitched_transcripts_tsv_file = str(tmp_path / ('%s_exacto_stitched_transcripts.tsv' % stitch_prefix))

    # Step 9. Realign the stitched transcripts (minimap2: the committed alignment)
    stitched_bam_file = get_data_path(name='assembly/rna/stitched/%s_exacto_stitched_transcripts_minimap2_sorted.bam' % stitch_prefix)
    df_stitched = pd.read_csv(stitched_transcripts_tsv_file, sep='\t')
    _check_alignment_holds(bam_file=stitched_bam_file, sequences=dict(zip(df_stitched['read_name'].astype(str), df_stitched['stitched_sequence'])))

    # Step 10. Call RNA transcript variants (scripts/data/05_variant_calling/03)
    identify_rna_transcript_variants(
        bam_file=stitched_bam_file,
        reference_genome_fasta_file=FASTA_FILE,
        output_dir=output_dir,
        output_prefix=rna_sample_id,
        num_threads=4,
        output_type=OutputType.FILE,
        **GENE_ANNOTATION
    )
    rna_prefix = str(tmp_path / ('%s_exacto_assembled_transcript' % rna_sample_id))

    # Step 11. Quantify the clusters
    df_abundances, _ = quantify_rna_abundances(
        clusters_tsv_file=pass2_clusters + '.tsv',
        cluster_reference_matches_tsv_file=pass2_clusters + '_reference_transcripts.tsv',
        cluster_splice_junctions_tsv_file=pass2_clusters + '_splice_junctions.tsv',
        cluster_variants_tsv_file=pass2_clusters + '_variants_passed.tsv',
        consensus_rna_reference_matches_tsv_file=rna_prefix + '_reference_transcript_matches.tsv',
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    # Step 12. Integrate DNA and RNA variants (scripts/data/06_integration/01)
    integrated_variants_tsv_file = str(tmp_path / ('%s_exacto_integrated_variants.tsv' % rna_sample_id))
    integrate_variants(
        dna_variants_tsv_file=somatic_dna_variants_tsv_file,
        rna_variants_tsv_file=rna_prefix + '_variants.tsv',
        output_tsv_file=integrated_variants_tsv_file,
        num_threads=1,
        output_type=OutputType.FILE,
        **GENE_ANNOTATION
    )

    # Step 13. Translate the transcripts (scripts/data/07_translation/01)
    translate_transcripts(
        assembled_transcript_model_alignments_tsv_file=rna_prefix + '_model_alignments.tsv',
        assembled_transcript_variants_tsv_file=rna_prefix + '_variants.tsv',
        dna_variants_tsv_file=somatic_dna_variants_tsv_file,
        integrated_variants_tsv_file=integrated_variants_tsv_file,
        rna_consensus_tsv_file=rna_consensus_tsv_file,
        stitched_transcripts_tsv_file=stitched_transcripts_tsv_file,
        strategy=TranslationStrategy.LONGEST_ORF,
        output_dir=output_dir,
        output_prefix=rna_sample_id,
        num_threads=2,
        output_type=OutputType.FILE
    )
    proteoforms_tsv_file = str(tmp_path / ('%s_exacto_proteoforms.tsv' % rna_sample_id))

    # Step 14. Find mutant peptides absent from a reference proteome
    df_peptide_variants = identify_peptide_variants(
        proteoforms_tsv_file=proteoforms_tsv_file,
        reference_proteome_fasta_file=get_data_path(name='exacto/exactolib/reference_peptides.fasta'),
        min_k=8,
        max_k=11,
        num_processes=1
    )

    # Every simulated peptide is a proteoform, and there are no others.
    df_proteoforms = pd.read_csv(proteoforms_tsv_file, sep='\t', keep_default_na=False)
    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-pep-%s-tumor.fasta' % sample_number)) as fasta:
        expected_sequences = sorted(record.sequence.replace('*', '') for record in fasta)
    assert sorted(df_proteoforms['amino_acid_sequence'].str.replace('*', '', regex=False)) == expected_sequences

    # Each cluster is quantified, and only proteoforms with mutant amino acids yield mutant peptides.
    assert sorted(df_abundances['cluster_id']) == sorted(df_consensus['cluster_id'])
    mutant_proteoform_ids = set(df_proteoforms.loc[df_proteoforms['num_mutant_amino_acids'] > 0, 'proteoform_id'])
    assert set(df_peptide_variants['proteoform_id']) <= mutant_proteoform_ids
