import pandas as pd
import pysam
import pytest
from test.data import get_data_path, get_simulated_read_origins
from exactolib.constants import GeneAnnotationSource
from exactolib.main import correct_rna_reads


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]


def _reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def _correct_rna_reads(sample_number: str, output_dir: str) -> str:
    # As scripts/data/04_assembly/rna/02_run_exacto_correct_rna_reads_pass1.sh runs it, on the pass-1 clusters
    prefix = 'assembly/rna/pass1/scga-mini-rna-%s-tumor_pass1_exacto_rna_clusters' % sample_number
    return correct_rna_reads(
        bam_file=get_data_path(name='alignment/scga-mini-rna-%s-tumor_minimap2_sorted.bam' % sample_number),
        clusters_tsv_file=get_data_path(name=prefix + '.tsv'),
        cluster_reference_transcripts_tsv_file=get_data_path(name=prefix + '_reference_transcripts.tsv'),
        cluster_splice_junctions_tsv_file=get_data_path(name=prefix + '_splice_junctions.tsv'),
        cluster_variants_tsv_file=get_data_path(name=prefix + '_variants_passed.tsv'),
        output_dir=output_dir,
        output_prefix='scga-mini-rna-%s-tumor_pass1' % sample_number,
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        num_threads=4
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_correct_rna_reads_restores_simulated_transcript_sequence(sample_number, tmp_path):
    corrected_reads_fastq_file = _correct_rna_reads(sample_number=sample_number, output_dir=str(tmp_path))

    assert corrected_reads_fastq_file == str(tmp_path / ('scga-mini-rna-%s-tumor_pass1_exacto_rna_corrected_reads.fastq.gz' % sample_number))
    with pysam.FastxFile(corrected_reads_fastq_file) as fastq:
        corrected_reads = {record.name: record.sequence for record in fastq}

    # One corrected read for each clustered read.
    df_clusters = pd.read_csv(get_data_path(name='assembly/rna/pass1/scga-mini-rna-%s-tumor_pass1_exacto_rna_clusters.tsv' % sample_number), sep='\t')
    assert set(corrected_reads) == set(df_clusters['read_name'])

    # A read is error-free when it is a stretch of the transcript it was simulated from. Before
    # correction 0-14% of each sample's reads are; after, 65% (013) to 91% (009).
    read_origins = get_simulated_read_origins(sample_id='scga-mini-rna-%s-tumor' % sample_number)
    with open(get_data_path(name='simulation/transcript/scga-mini-rna-%s-tumor.transcript' % sample_number)) as handle:
        transcript_sequences = {fields[0]: fields[3].upper() for fields in (line.rstrip('\n').split('\t') for line in handle if line.strip())}
    with pysam.AlignmentFile(get_data_path(name='alignment/scga-mini-rna-%s-tumor_minimap2_sorted.bam' % sample_number)) as bam:
        sequenced_reads = {record.query_name: (_reverse_complement(record.query_sequence) if record.is_reverse else record.query_sequence)
                           for record in bam if not (record.is_unmapped or record.is_secondary or record.is_supplementary)}

    def is_error_free(read_name: str, sequence: str) -> bool:
        transcript_sequence = transcript_sequences[read_origins[read_name]]
        return sequence in transcript_sequence or _reverse_complement(sequence) in transcript_sequence

    num_error_free_corrected = sum(is_error_free(name, sequence) for name, sequence in corrected_reads.items())
    num_error_free_sequenced = sum(is_error_free(name, sequenced_reads[name]) for name in corrected_reads)
    assert num_error_free_corrected >= 0.6 * len(corrected_reads)
    assert num_error_free_corrected > 2 * num_error_free_sequenced


def test_correct_rna_reads_keeps_each_reads_allele(tmp_path):
    # scga-mini-rna-001 simulates TP53 with (row 1) and without (row 2) one SNV. Correction rewrites
    # reads toward their own cluster, so a mutant read keeps the SNV and a wild-type read gains none.
    corrected_reads_fastq_file = _correct_rna_reads(sample_number='001', output_dir=str(tmp_path))
    read_origins = get_simulated_read_origins(sample_id='scga-mini-rna-001-tumor')
    with open(get_data_path(name='simulation/transcript/scga-mini-rna-001-tumor.transcript')) as handle:
        transcript_sequences = {fields[0]: fields[3].upper() for fields in (line.rstrip('\n').split('\t') for line in handle if line.strip())}
    mutant_sequence = transcript_sequences['scga-mini-rna-001-tumor-1']
    wild_type_sequence = transcript_sequences['scga-mini-rna-001-tumor-2']
    snv_index = next(i for i, (a, b) in enumerate(zip(mutant_sequence, wild_type_sequence)) if a != b)
    alleles = {
        'scga-mini-rna-001-tumor-1': mutant_sequence[snv_index - 10:snv_index + 11],
        'scga-mini-rna-001-tumor-2': wild_type_sequence[snv_index - 10:snv_index + 11]
    }

    num_reads = {'scga-mini-rna-001-tumor-1': 0, 'scga-mini-rna-001-tumor-2': 0}
    with pysam.FastxFile(corrected_reads_fastq_file) as fastq:
        for record in fastq:
            origin = read_origins[record.name]
            for name, allele in alleles.items():
                if allele in record.sequence or _reverse_complement(allele) in record.sequence:
                    assert name == origin, record.name
                    num_reads[name] += 1
    assert num_reads['scga-mini-rna-001-tumor-1'] > 100 and num_reads['scga-mini-rna-001-tumor-2'] > 100
