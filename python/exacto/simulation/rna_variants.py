# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement functions related to
simulating RNA variants.
"""


import pysam
from ..utilities.utils import *
from ..constants import *
from ..utilities.gencode_utils import *
from ..wrappers.gene import *
from ..wrappers.transcript import *
from ..wrappers.exon import *
from .rna.match import Match
from .rna.deletion import Deletion
from .rna.insertion import Insertion
from .rna.substitution import Substitution
from .utils import *


def build_variant_transcript_sequence(
        genome_fasta: pysam.FastaFile,
        df_transcript_variants: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame):
    """
    Builds variant transcript sequence.

    Parameters
    ----------
    genome_fasta                :   pysam.FastaFile.
    df_transcript_variants      :   DataFrame of all variants
                                    pertaining to a specific transcript.
    df_transcripts              :   DataFrame of reference transcripts.
    df_exons                    :   DataFrame of reference exons.

    Returns
    -------

    """
    # Sort the variants
    df_transcript_variants = df_transcript_variants.sort_values(['exon_number'], ascending=True)

    variant_exons = []
    for exon_id, df_exon_variants in df_transcript_variants.groupby('exon_id'):
        # Step 1. Fetch the reference exon and transcript
        ref_exon = df_exons.loc[df_exons['exon_id'] == exon_id,:]
        ref_transcript = df_transcripts.loc[df_transcripts['transcript_id'] == ref_exon['transcript_id'].values.tolist()[0],:]

        # Step 2. Fetch the exon sequence
        exon_sequence = genome_fasta.fetch(
            ref_exon['exon_chrom'].values.tolist()[0],
            ref_exon['exon_start'].values.tolist()[0] - 1,
            ref_exon['exon_end'].values.tolist()[0]
        )

        # Step 3. Initialize a variant exon to the reference exon
        variant_exon = Match(
            gene_id=ref_exon['gene_id'].values.tolist()[0],
            transcript_id=ref_exon['transcript_id'].values.tolist()[0],
            exon_id=ref_exon['exon_id'].values.tolist()[0],
            exon_number=ref_exon['exon_number'].values.tolist()[0],
            strand=ref_transcript['transcript_strand'].values.tolist()[0],
            chrom=ref_exon['exon_chrom'].values.tolist()[0],
            start=ref_exon['exon_start'].values.tolist()[0],
            end=ref_exon['exon_end'].values.tolist()[0],
            length=ref_exon['exon_end'].values.tolist()[0] - ref_exon['exon_start'].values.tolist()[0],
            sequence=exon_sequence
        )

        # Step 4. Sort by variant position in ascending order
        df_exon_variants = df_exon_variants.sort_values(['pos_1'], ascending=True)

        # Step 5. Apply variants
        for _, row in df_exon_variants.iterrows():
            if row['variant_type'] == VariantTypes.SNV:
                variant_exon = Substitution(
                    ref_exon=variant_exon,
                    snv_pos=row['pos_1'],
                    snv_alt=row['alt']
                )
            if row['variant_type'] == VariantTypes.DELETION:
                variant_exon = Deletion(
                    ref_exon=variant_exon,
                    del_start=row['pos_1'],
                    del_end=row['pos_2']
                )
            if row['variant_type'] == VariantTypes.INSERTION:
                variant_exon = Insertion(
                    ref_exon=variant_exon,
                    ins_pos=row['pos_1'],
                    ins_sequence=row['alt']
                )
        variant_exons.append(variant_exon)

    # Build variant transcript sequence
    variant_transcript_sequence = ''
    for curr_variant_exon in variant_exons:
        variant_transcript_sequence += curr_variant_exon.sequence

    return variant_transcript_sequence


def generate_single_nucleotide_rna_variants(
        genome_fasta: pysam.FastaFile,
        num_snv: int,
        df_genes: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame,
        df_rna_variants: pd.DataFrame) -> pd.DataFrame:
    """
    Simulates single-nucleotide RNA variants.

    Parameters
    ----------
    genome_fasta        :   pysam.FastaFile object.
    num_snv             :   Number of SNVs to simulate.
    df_genes            :   DataFrame of genes.
    df_transcripts      :   DataFrame of transcripts.
    df_exons            :   DataFrame of exons.
    df_rna_variants     :   DataFrame of RNA variants.

    Returns
    -------
    df_rna_variants     :   DataFrame with the following columns:
                            'variant_id'
                            'variant_type'
                            'gene_id'
                            'transcript_id'
                            'exon_id'
                            'chr_1'
                            'pos_1'
                            'chr_2'
                            'pos_2'
                            'ref'
                            'alt'
                            'variant_sequence'
    """
    variant_idx = 1
    while True:
        # Step 1. Randomly select a SNV position
        rna_pos = randomly_select_rna_position(
            df_genes=df_genes,
            df_transcripts=df_transcripts,
            df_exons=df_exons
        )

        # Step 2. Skip if the randomly selected position is already in the list of variants
        df_rna_variants_curr_exon = df_rna_variants.loc[
            df_rna_variants['gene_id'] == rna_pos.gene_id &
            df_rna_variants['transcript_id'] == rna_pos.transcript_id &
            df_rna_variants['exon_id'] == rna_pos.exon_id,:
        ]
        overlaps = overlaps_any(df=df_rna_variants_curr_exon,
                                chrom=rna_pos.chrom,
                                start=rna_pos.pos,
                                end=rna_pos.pos)
        if overlaps:
            continue

        # Step 3. Query position nucleotide
        ref_allele = genome_fasta.fetch(
            rna_pos.chrom,
            rna_pos.pos - 1,
            rna_pos.pos
        )
        if rna_pos.strand == '-':
            ref_allele = get_complement_nucleotide(nucleotide=ref_allele)
        alt_allele = generate_single_nucleotide_variant(reference_allele=ref_allele)

        # Step 4. Append SNV to df_rna_variants
        df_temp = pd.DataFrame({
            'variant_id': ['variant_' + str(variant_idx)],
            'variant_type': [VariantTypes.SNV],
            'gene_id': [rna_pos.gene_id],
            'transcript_id': [rna_pos.transcript_id],
            'transcript_strand': [rna_pos.strand],
            'exon_id': [rna_pos.exon_id],
            'chr_1': [rna_pos.chrom],
            'pos_1': [rna_pos.pos],
            'chr_2': [rna_pos.chrom],
            'pos_2': [rna_pos.pos],
            'ref': [ref_allele],
            'alt': [alt_allele],
            'variant_sequence': [alt_allele]
        })
        df_rna_variants = pd.concat([df_rna_variants, df_temp])
        variant_idx += 1

        if variant_idx == num_snv:
            break
    return df_rna_variants


def generate_insertion_rna_variants(
        num_insertion: int,
        insertion_size_mean: int,
        insertion_size_stdev: int,
        df_genes: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame,
        df_rna_variants: pd.DataFrame) -> pd.DataFrame:
    """
    Simulates single-nucleotide RNA variants.

    Parameters
    ----------
    num_insertion           :   Number of insertions to simulate.
    insertion_size_mean      :  Mean value of insertion size.
    insertion_size_stdev    :   Standard deviation of insertion size.
    df_genes                :   DataFrame of genes.
    df_transcripts          :   DataFrame of transcripts.
    df_exons                :   DataFrame of exons.
    df_rna_variants         :   DataFrame of RNA variants.

    Returns
    -------
    df_rna_variants         :   DataFrame with the following columns:
                                'variant_id'
                                'variant_type'
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
                                'ref'
                                'alt'
                                'variant_sequence'
    """
    variant_idx = 1
    while True:
        # Step 1. Randomly select an insertion site and an insertion sequence
        rna_pos = randomly_select_rna_position(
            df_genes=df_genes,
            df_transcripts=df_transcripts,
            df_exons=df_exons
        )
        alt_allele = generate_insertion(insertion_size_mean=insertion_size_mean,
                                        insertion_size_stdev=insertion_size_stdev)

        # Step 2. Skip if the randomly selected position is already in the list of variants
        df_rna_variants_curr_exon = df_rna_variants.loc[
            df_rna_variants['gene_id'] == rna_pos.gene_id &
            df_rna_variants['transcript_id'] == rna_pos.transcript_id &
            df_rna_variants['exon_id'] == rna_pos.exon_id,:
        ]
        overlaps = overlaps_any(df=df_rna_variants_curr_exon,
                                chrom=rna_pos.chrom,
                                start=rna_pos.pos,
                                end=rna_pos.pos)
        if overlaps:
            continue

        # Step 3. Append insertion to df_rna_variants
        df_temp = pd.DataFrame({
            'variant_id': ['variant_' + str(variant_idx)],
            'variant_type': [VariantTypes.INSERTION],
            'gene_id': [rna_pos.gene_id],
            'transcript_id': [rna_pos.transcript_id],
            'transcript_strand': [rna_pos.strand],
            'exon_id': [rna_pos.exon_id],
            'chr_1': [rna_pos.chrom],
            'pos_1': [rna_pos.pos],
            'chr_2': [rna_pos.chrom],
            'pos_2': [rna_pos.pos],
            'ref': [''],
            'alt': [alt_allele],
            'variant_sequence': [alt_allele]
        })
        df_rna_variants = pd.concat([df_rna_variants, df_temp])
        variant_idx += 1

        if variant_idx == num_insertion:
            break
    return df_rna_variants


def generate_deletion_rna_variants(
        genome_fasta: pysam.FastaFile,
        num_deletion: int,
        deletion_size_mean: int,
        deletion_size_stdev: int,
        df_genes: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame,
        df_rna_variants: pd.DataFrame) -> pd.DataFrame:
    """
    Simulates single-nucleotide RNA variants.

    Parameters
    ----------
    genome_fasta            :   pysam.FastaFile object.
    num_deletion            :   Number of deletions to simulate.
    deletion_size_mean      :   Mean value of deletion size.
    deletion_size_stdev     :   Standard deviation of deletion size.
    df_genes                :   DataFrame of genes.
    df_transcripts          :   DataFrame of transcripts.
    df_exons                :   DataFrame of exons.
    df_rna_variants         :   DataFrame of RNA variants.

    Returns
    -------
    df_rna_variants     :   DataFrame with the following columns:
                            'variant_id'
                            'variant_type'
                            'gene_id'
                            'transcript_id'
                            'exon_id'
                            'chr_1'
                            'pos_1'
                            'chr_2'
                            'pos_2'
                            'ref'
                            'alt'
                            'variant_sequence'
    """
    variant_idx = 1
    while True:
        # Step 1. Randomly select deletion position, type, and size
        rna_pos = randomly_select_rna_position(
            df_genes=df_genes,
            df_transcripts=df_transcripts,
            df_exons=df_exons
        )

        # Step 2. Generate a deletion
        deletion_start, deletion_end, deletion_size = generate_rna_deletion(
            rna_pos=rna_pos,
            deletion_size_mean=deletion_size_mean,
            deletion_size_stdev=deletion_size_stdev
        )

        # Step 3. Skip if the randomly selected deletion is already in the list of variants
        df_rna_variants_curr_exon = df_rna_variants.loc[
            df_rna_variants['gene_id'] == rna_pos.gene_id &
            df_rna_variants['transcript_id'] == rna_pos.transcript_id &
            df_rna_variants['exon_id'] == rna_pos.exon_id,:
        ]
        overlaps = overlaps_any(df=df_rna_variants_curr_exon,
                                chrom=rna_pos.chrom,
                                start=deletion_start,
                                end=deletion_end)
        if overlaps:
            continue

        # Step 4. Query the deletion sequence
        ref_allele = genome_fasta.fetch(
            rna_pos.chrom,
            deletion_start - 1,
            deletion_end
        )

        # Step 5. Append deletion to df_rna_variants
        df_temp = pd.DataFrame({
            'variant_id': ['variant_' + str(variant_idx)],
            'variant_type': [VariantTypes.DELETION],
            'gene_id': [rna_pos.gene_id],
            'transcript_id': [rna_pos.transcript_id],
            'transcript_strand': [rna_pos.strand],
            'exon_id': [rna_pos.exon_id],
            'chr_1': [rna_pos.chrom],
            'pos_1': [deletion_start],
            'chr_2': [rna_pos.chrom],
            'pos_2': [deletion_end],
            'ref': [ref_allele],
            'alt': [''],
            'variant_sequence': ['']
        })
        df_rna_variants = pd.concat([df_rna_variants, df_temp])
        variant_idx += 1

        if variant_idx == num_deletion:
            break
    return df_rna_variants


def simulate_rna_variants(
        genome_fasta: pysam.FastaFile,
        df_genes: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame,
        df_target_regions: pd.DataFrame,
        df_herv_regions: pd.DataFrame,
        num_snv: int,
        num_insertion: int,
        num_deletion: int,
        num_fusion: int,
        num_inversion: int,
        num_herv: int,
        insertion_size_mean: int,
        insertion_size_stdev: int,
        deletion_size_mean: int,
        deletion_size_stdev: int,
        herv_solo_ltr_proportion: float,
        herv_truncated_proportion: float,
        herv_chimeric_proportion: float,
        herv_chimeric_max_neighboring_distance: int,
        herv_full_length_proportion: float,
        infinite_sites_assumption: bool) -> Tuple[pd.DataFrame, List[str]]:
    """
    Simulates RNA variants.

    Parameters
    ----------
    genome_fasta                            :   pysam.FastaFile object of reference genome.
    df_transcripts                          :   DataFrame of transcripts.
    df_exons                                :   DataFrame of exons.
    df_target_regions                       :   DataFrame of regions to simulate RNA variants.
    df_herv_regions                         :   HERV regions.
    num_snv                                 :   Number of SNVs to simulate.
    num_insertion                           :   Number of insertions to simulate.
    num_deletion                            :   Number of deletions to simulate.
    num_fusion                              :   Number of fusions to simulate.
    num_inversion                           :   Number of inversions to simulate.
    num_herv                                :   Number of HERVs to simulate.
    herv_solo_ltr_proportion                :   Proportion of expressed HERVs that only have solo LTR sequences.
    herv_truncated_proportion               :   Proportion of HERVs that are truncated.
    herv_chimeric_proportion                :   Proportion of HERVs that are chimeric (concatenation of neighboring HERVs).
    herv_chimeric_max_neighboring_distance  :   Maximum distance for two HERVs to be considered for simulation of a
                                                chimeric HERV.
    herv_full_length_proportion             :   Proportion of HERVs that are full-lengths.
    infinite_sites_assumption               :   If true, the simulation enforces the infinite sites assumption.

    Returns
    -------
    df_rna_variants                         :   DataFrame of RNA variants
    variant_transcript_sequences            :   List of variant transcript sequences
    """
    # Step 1. Check if there are desired genomic regions
    if df_target_regions:
        df_genes, df_transcripts, df_exons = subset_gencode_dataframes(
            df_target_regions=df_target_regions,
            df_genes=df_genes,
            df_transcripts=df_transcripts,
            df_exons=df_exons
        )

    # Step 2. Generate single-nucleotide RNA variants
    df_rna_variants = generate_single_nucleotide_rna_variants(
        genome_fasta=genome_fasta,
        num_snv=num_snv,
        df_rna_variants=pd.DataFrame(),
        df_genes=df_genes,
        df_transcripts=df_transcripts,
        df_exons=df_exons
    )

    # Step 3. Generate insertion RNA variants
    df_insertions = generate_insertion_rna_variants(
        num_insertion=num_insertion,
        insertion_size_mean=insertion_size_mean,
        insertion_size_stdev=insertion_size_stdev,
        df_genes=df_genes,
        df_transcripts=df_transcripts,
        df_exons=df_exons,
        df_rna_variants=df_rna_variants
    )
    df_rna_variants = pd.concat([df_rna_variants, df_insertions])

    # Step 4. Generate deletion RNA variants
    df_deletions = generate_deletion_rna_variants(
        genome_fasta=genome_fasta,
        num_deletion=num_deletion,
        deletion_size_mean=deletion_size_mean,
        deletion_size_stdev=deletion_size_stdev,
        df_genes=df_genes,
        df_transcripts=df_transcripts,
        df_exons=df_exons,
        df_rna_variants=df_rna_variants
    )
    df_rna_variants = pd.concat([df_rna_variants, df_deletions])

    # Step X. Build variant transcript sequences
    variant_transcript_sequences = []
    for name, group in df_rna_variants.groupby('transcript_id'):
        variant_transcript_sequence = build_variant_transcript_sequence(
            genome_fasta=genome_fasta,
            df_transcript_variants=group,
            df_genes=df_genes,
            df_transcripts=df_transcripts,
            df_exons=df_exons
        )
        variant_transcript_sequences.append(variant_transcript_sequence)

    return df_rna_variants, variant_transcript_sequences


