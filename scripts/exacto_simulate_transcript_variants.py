#!/usr/bin/python3

"""
The purpose of this python3 script is to randomly generate variant transcripts
and to output the following files:

1.  One FASTA file of reference and simulated variant transcript sequences.
2.  One TSV file of variant transcript SV breakpoints and ground truth on
    the abundance of each transcript simulated.

Author: Jin Seok (Andy) Lee

Last updated date: June 9, 2022
"""


import argparse
from exactolib.simulation.common import *
from exactolib.simulation.snv_indel import *
from exactolib.utilities.file_reader import *


logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s [%(funcName)s] %(message)s")


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""Simulates transcript variants."""
    )
    arg_parser.add_argument(
        "--gencode_pc_transcripts_fasta_file",
        dest="gencode_pc_transcripts_fasta_file",
        type=str,
        required=True,
        help="GENCODE protein-coding transcripts FASTA file."
    )
    arg_parser.add_argument(
        "--gencode_genome_reference_fasta_file",
        dest="gencode_genome_reference_fasta_file",
        type=str,
        required=True,
        help="GENCODE primary assembly genome reference (primary assembly) FASTA file."
    )
    arg_parser.add_argument(
        "--gencode_gtf_file",
        dest="gencode_gtf_file",
        type=str,
        required=True,
        help="GENCODE GTF file."
    )
    arg_parser.add_argument(
        "--gencode_refseq_file",
        dest="gencode_refseq_file",
        type=str,
        required=True,
        help="GENCODE RefSeq file."
    )
    arg_parser.add_argument(
        "--num_snvs",
        dest="num_snvs",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with SNVs (default: 200)."
    )
    arg_parser.add_argument(
        "--num_indels",
        dest="num_indel",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with INDELs (default: 200)."
    )
    arg_parser.add_argument(
        "--num_deletions",
        dest="num_deletions",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with deletions (default: 200)."
    )
    arg_parser.add_argument(
        "--num_insertions",
        dest="num_insertions",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with insertions (default: 200)."
    )
    arg_parser.add_argument(
        "--num_duplications",
        dest="num_duplications",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with duplications (default: 200)."
    )
    arg_parser.add_argument(
        "--num_inversions",
        dest="num_inversions",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with inversions (default: 200)."
    )
    arg_parser.add_argument(
        "--num_translocations",
        dest="num_translocations",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with translocations (default: 200)."
    )
    arg_parser.add_argument(
        "--num_viral_integrations",
        dest="num_viral_integrations",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with viral integrations (default: 200)."
    )
    arg_parser.add_argument(
        "--num_viral_integrations",
        dest="num_viral_integrations",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with viral integrations (default: 200)."
    )
    arg_parser.add_argument(
        "--num_intron_retentions",
        dest="num_intron_retentions",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts with intron retention (default: 200)."
    )
    arg_parser.add_argument(
        "--num_circular_rnas",
        dest="num_circular_rnas",
        type=int,
        default=200,
        required=True,
        help="Number of variant transcripts that are circular RNAs (default: 200)."
    )
    arg_parser.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    args = arg_parser.parse_args()
    return args


if __name__ == "__main__":
    # Step 1. Read input arguments
    args = parse_args()

    # Step 2. Read GTF file
    df_gtf = read_gencode_gtf_file(gencode_gtf_file=args.gencode_gtf_file)

    # Step 3. Read GENCODE RefSeq metadata file
    df_refseq = read_gencode_refseq_file(gencode_refseq_metadata_file=args.gencode_refseq_file)

    # Step 4. Read transcripts FASTA file
    df_transcript_sequences = read_gencode_transcripts_fasta_file(gencode_fasta_file=args.gencode_pc_transcripts_fasta_file)

    # Step 5. Identify common transcript IDs
    annotation_transcript_ids = df_gtf.loc[df_gtf['transcript_type'] == 'protein_coding', 'transcript_id'].unique()
    refseq_transcript_ids = df_refseq['ensembl_transcript_id'].unique()
    pc_transcript_ids = df_transcript_sequences['ensembl_transcript_id'].unique()
    eligible_transcript_ids = set.intersection(set(annotation_transcript_ids), set(refseq_transcript_ids))
    eligible_transcript_ids = set.intersection(set(eligible_transcript_ids), set(pc_transcript_ids))
    logging.info("%i protein-coding transcripts have been identified." % len(eligible_transcript_ids))

    # Step 6. Simulate transcript variants
    data = {
        'variant_id': [],
        'variant_type': [],
        'transcript_id': [],
        'chrom_1': [],
        'pos_1': [],
        'chrom_2': [],
        'pos_2': [],
        'transcript_strand': [],
        'variant_transcript_sequence': [],
        'variant_exon_number': [],
        'ref_allele': [],
        'var_allele': []
    }

    transcript_ids_list = list(eligible_transcript_ids)

    # SNVs
    for i in range(0, 200):
        # Randomly pick a transcript ID
        transcript_id = random.choice(transcript_ids_list)
        transcript_ids_list.remove(transcript_id)

        # Generate a random SNV
        variant_type, \
        strand, \
        variant_transcript_sequence, \
        variant_chrom, \
        variant_position, \
        exon_number, \
        ref_allele, \
        var_allele = randomly_generate_snv(
            transcript_id=transcript_id,
            df_gtf=df_gtf,
            reference_genome_fasta_file=args.gencode_genome_reference_fasta_file
        )

        # Append data
        data['variant_id'].append('snv_' + str(i + 1))
        data['variant_type'].append(variant_type)
        data['transcript_id'].append(transcript_id)
        data['chrom_1'].append(variant_chrom)
        data['pos_1'].append(variant_position)
        data['chrom_2'].append(variant_chrom)
        data['pos_2'].append(variant_position)
        data['transcript_strand'].append(strand)
        data['variant_transcript_sequence'].append(variant_transcript_sequence)
        data['variant_exon_number'].append(exon_number)
        data['ref_allele'].append(ref_allele)
        data['var_allele'].append(var_allele)

    df = pd.DataFrame(data)
    print(df.head(n=20))

    # INDELs
