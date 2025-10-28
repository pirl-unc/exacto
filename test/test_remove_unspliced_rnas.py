import tempfile
import pysam
from .data import get_data_path
from exactolib.constants import GeneAnnotationSource
from exactolib.main import remove_unspliced_rnas


def test_remove_unspliced_rnas_1():
    output_bam_file = tempfile.mktemp(suffix=".bam")
    output_bai_file = tempfile.mktemp(suffix=".bai")
    output_fasta_file = tempfile.mktemp(suffix=".fasta")

    remove_unspliced_rnas(
        bam_file=get_data_path(name='bam/rna-100-tumor_minimap2_mdtagged_sorted_dna-001-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-100-tumor_minimap2_mdtagged_sorted_dna-001-tumor_minimap2_mdtagged_sorted.bam.bai'),
        fasta_file=get_data_path(name='fasta/rna-100-tumor_long-read_dna-001-tumor_long-read.fasta'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        gene_types=['protein_coding'],
        gene_levels=[1,2],
        transcript_types=['protein_coding'],
        transcript_levels=[1,2],
        output_bam_file=output_bam_file,
        output_bam_bai_file=output_bai_file,
        output_fasta_file=output_fasta_file,
        num_threads=1,
        min_mapping_quality=30
    )

    fasta = pysam.FastaFile(output_fasta_file)

    assert fasta.nreferences == 2
