#!/usr/bin/python3

"""
The purpose of this python3 script is to randomly generate variant transcripts
and to output the following files:

1.  One FASTA file of reference and simulated variant transcript sequences.
2.  One TSV file of variant transcript SV breakpoints and ground truth on
    the abundance of each transcript simulated.
3.  One long-read RNA-seq FASTQ file of reference and variant transcripts.

Author: Jin Seok (Andy) Lee

Last updated date: June 9, 2022
"""


import argparse
import datetime
import pandas as pd
import numpy as np
import random
import math
import multiprocessing as mp
from multiprocessing import Process, Manager


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""
        Simulates long-read RNA (PacBio Iso-seq) sequencing FASTQ file based on
        a list of variants.
        """
    )
    arg_parser.add_argument(
        "--gencode_transcripts_fasta_file",
        dest="gencode_transcripts_fasta_file",
        type=str,
        required=True,
        help="GENCODE transcripts FASTA file."
    )
    arg_parser.add_argument(
        "--bailey_et_al_cell_2018_s1_excel_file",
        dest="bailey_et_al_cell_2018_s1_excel_file",
        type=str,
        required=True,
        help="Bailey et al., Cell 2018, supplementary table S1 file."
    )
    arg_parser.add_argument(
        "--read_length_average",
        dest="read_length_average",
        type=int,
        default=10000,
        required=True,
        help="Simulated read length average (default: 10,000)."
    )
    arg_parser.add_argument(
        "--read_length_stdev",
        dest="read_length_stdev",
        type=int,
        default=1000,
        required=True,
        help="Simulated read length standard deviation (default: 1,000)."
    )
    arg_parser.add_argument(
        "--total_num_reads",
        dest="total_num_reads",
        type=int,
        default=500000,
        required=True,
        help="Total number of reads to simulate (default: 500,000)."
    )
    arg_parser.add_argument(
        "--base_phred_score_average",
        dest="base_phred_score_average",
        type=int,
        default=100,
        required=True,
        help="Simulated read base phred score average (default: 100)."
    )
    arg_parser.add_argument(
        "--base_phred_score_stdev",
        dest="base_phred_score_stdev",
        type=int,
        default=10,
        required=True,
        help="Simulated read base phred score standard deviation (default: 10)."
    )
    arg_parser.add_argument(
        "--num_cores",
        dest="num_cores",
        type=int,
        default=4,
        required=True,
        help="Number of cores to use (default: 10)."
    )
    arg_parser.add_argument(
        "--output_dir",
        dest="output_dir",
        type=str,
        required=True,
        help="Output directory."
    )
    arg_parser.add_argument(
        "--sample_id",
        dest="sample_id",
        type=str,
        required=True,
        help="Sample ID."
    )
    args = arg_parser.parse_args()
    return args


def read_gencode_fasta_file(fasta_file):
    data = {
        'Ensembl_Transcript_ID': [],
        'Ensembl_Gene_ID': [],
        'Havana_Transcript_ID': [],
        'Havana_Gene_ID': [],
        'Gene_Symbol_Versioned': [],
        'Gene_Symbol': [],
        'Transcript_Length': [],
        'Transcript_Type': [],
        'Transcript_Sequence': []
    }
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        curr_transcript_sequence = ''
        first = True
        for line in lines:
            line = line.strip()
            if '>' in line:
                # Dump previously stored sequence
                if first == True:
                    first = False
                else:
                    data['Transcript_Sequence'].append(curr_transcript_sequence)
                    curr_transcript_sequence = ''

                # Append new transcript info
                line = line[1:]
                line_elements = line.split('|')
                data['Ensembl_Transcript_ID'].append(line_elements[0])
                data['Ensembl_Gene_ID'].append(line_elements[1])
                data['Havana_Gene_ID'].append(line_elements[2])
                data['Havana_Transcript_ID'].append(line_elements[3])
                data['Gene_Symbol_Versioned'].append(line_elements[4])
                data['Gene_Symbol'].append(line_elements[5])
                data['Transcript_Length'].append(int(line_elements[6]))
                data['Transcript_Type'].append(line_elements[7])
            else:
                curr_transcript_sequence = curr_transcript_sequence + line

        # Dump sequence of the last transcript
        data['Transcript_Sequence'].append(curr_transcript_sequence)
    df_gencode = pd.DataFrame(data)
    return df_gencode


def read_bailey_oncogenes_and_tsgs_list_file(excel_file):
    df_bailey_s1 = pd.read_excel(excel_file)
    df_bailey_s1 = df_bailey_s1.loc[
        (df_bailey_s1['Tumor suppressor or oncogene prediction (by 20/20+)'].isin(
            ['oncogene', 'tsg']
        )) &
        (df_bailey_s1['Decision'] == 'official'),:
    ]
    print(len(df_bailey_s1['Gene'].unique()), 'oncogenes and TSGs from Bailey et al., Cell 2018')
    return df_bailey_s1


def generate_variant_transcript_sequences(candidate_genes,
                                          df_gencode,
                                          deletion_size_average,
                                          deletion_size_stdev,
                                          insertion_size_average,
                                          insertion_size_stdev,
                                          duplication_size_average,
                                          duplication_size_stdev):
    df_gencode_variants = pd.DataFrame()
    for curr_gene in candidate_genes:
        df_matched = df_gencode.loc[(
            (df_gencode['Gene_Symbol'] == curr_gene) &
            (df_gencode['Transcript_Type'] == 'protein_coding')),:]

        if len(df_matched) == 0:
            continue

        # Take the longest transcript
        df_matched = df_matched.sort_values(by=['Transcript_Length'], ascending=False).iloc[0]

        # Randomly choose a SV type
        curr_sv_type_rand = random.randint(0, 2)
        if curr_sv_type_rand == 0:      # deletion
            variant_transcript_sequence, variant_start, variant_end = generate_variant_transcript_with_deletion(
                transcript_sequence=df_matched['Transcript_Sequence'],
                length_average=deletion_size_average,
                length_stdev=deletion_size_stdev
            )
            curr_sv_type = 'del'
        elif curr_sv_type_rand == 1:    # insertion
            variant_transcript_sequence, variant_start, variant_end = generate_variant_transcript_with_insertion(
                transcript_sequence=df_matched['Transcript_Sequence'],
                length_average=insertion_size_average,
                length_stdev=insertion_size_stdev
            )
            curr_sv_type = 'ins'
        else:                           # duplication
            variant_transcript_sequence, variant_start, variant_end = generate_variant_transcript_with_duplication(
                transcript_sequence=df_matched['Transcript_Sequence'],
                length_average=duplication_size_average,
                length_stdev=duplication_size_stdev
            )
            curr_sv_type = 'dup'

        if df_matched['Havana_Transcript_ID'] == '-':
            curr_havana_transcript_id = '-'
        else:
            curr_havana_transcript_id = df_matched['Havana_Transcript_ID'] + '_' + curr_sv_type

        df_temp = pd.DataFrame({
            'Ensembl_Transcript_ID': [df_matched['Ensembl_Transcript_ID'] + '_' + curr_sv_type],
            'Ensembl_Gene_ID': [df_matched['Ensembl_Gene_ID']],
            'Havana_Transcript_ID': [curr_havana_transcript_id],
            'Havana_Gene_ID': [df_matched['Havana_Gene_ID']],
            'Gene_Symbol_Versioned': [df_matched['Gene_Symbol_Versioned']],
            'Gene_Symbol': [df_matched['Gene_Symbol']],
            'Transcript_Length': [len(variant_transcript_sequence)],
            'Transcript_Type': ['protein_coding_variant'],
            'Transcript_Sequence': [variant_transcript_sequence],
            'Variant_Start': [variant_start],
            'Variant_End': [variant_end]
        })
        df_gencode_variants = pd.concat([df_gencode_variants, df_temp])
    print(len(df_gencode_variants), 'variant transcript sequences generated')
    return df_gencode_variants


def generate_variant_transcript_with_insertion(transcript_sequence,
                                               length_average,
                                               length_stdev):
    # Step 1. Generation an insertion sequence
    atcg = ['A', 'T', 'C', 'G']
    if length_average > len(transcript_sequence):
        length_average = len(transcript_sequence)
    insertion_length = math.ceil(np.random.normal(length_average, length_stdev))
    insertion_sequence = ''.join(random.choices(atcg, k=insertion_length))

    # Step 2. Choose an insertion site
    insertion_site = random.randint(1, len(transcript_sequence))

    # Step 3. Generate variant transcript with the insertion
    variant_transcript_sequence = \
        transcript_sequence[0:insertion_site] + \
        insertion_sequence + \
        transcript_sequence[insertion_site:]

    return variant_transcript_sequence, insertion_site, insertion_site


def generate_variant_transcript_with_deletion(transcript_sequence,
                                              length_average,
                                              length_stdev):
    # Step 1. Choose a deletion site
    deletion_start = random.randint(1, len(transcript_sequence))

    # Step 2. Choose a deletion length
    if length_average > len(transcript_sequence):
        length_average = len(transcript_sequence)
    deletion_length = math.ceil(np.random.normal(length_average, length_stdev))
    deletion_end = deletion_start + deletion_length

    # Step 3. Generate variant transcript with the deletion
    variant_transcript_sequence = \
        transcript_sequence[0:deletion_start] + \
        transcript_sequence[deletion_start + deletion_end:]

    return variant_transcript_sequence, deletion_start, deletion_end


def generate_variant_transcript_with_duplication(transcript_sequence,
                                                 length_average,
                                                 length_stdev):
    # Step 1. Choose a duplication site
    duplication_start = random.randint(1, len(transcript_sequence))

    # Step 2. Choose a duplication length
    if length_average > len(transcript_sequence):
        length_average = len(transcript_sequence)
    duplication_end = math.ceil(np.random.normal(length_average, length_stdev)) + duplication_start

    # Step 3. Generate variant transcript with the duplication
    variant_transcript_sequence = \
        transcript_sequence[0:duplication_start] + \
        transcript_sequence[duplication_start:duplication_end] + \
        transcript_sequence[duplication_start:duplication_end] + \
        transcript_sequence[duplication_end:]

    return variant_transcript_sequence, duplication_start, duplication_end


def write_transcript_sequences_to_fasta_file(df, output_fasta_file):
    with open(output_fasta_file, 'w') as f:
        for index, row in df.iterrows():
            curr_ensembl_transcript_id = row['Ensembl_Transcript_ID']
            curr_ensembl_gene_id = row['Ensembl_Gene_ID']
            curr_havana_transcript_id = row['Havana_Transcript_ID']
            curr_havana_gene_id = row['Havana_Gene_ID']
            curr_gene_symbol_versioned = row['Gene_Symbol_Versioned']
            curr_gene_symbol = row['Gene_Symbol']
            curr_transcript_length = row['Transcript_Length']
            curr_transcript_type = row['Transcript_Type']
            curr_transcript_sequence = row['Transcript_Sequence']
            f.write('>' + \
                    curr_ensembl_transcript_id + '|' + \
                    curr_ensembl_gene_id + '|' + \
                    curr_havana_gene_id + '|' + \
                    curr_havana_transcript_id + '|' + \
                    curr_gene_symbol_versioned + '|' + \
                    curr_gene_symbol + '|' + \
                    str(curr_transcript_length) + '|' + \
                    curr_transcript_type + '|' + '\n')
            f.write(curr_transcript_sequence + '\n')


def randomly_generate_weights(n):
    weights = []
    for i in range(0, n):
        random_max = random.randint(0, 101)
        random_number = random.randint(0, random_max)
        weights.append(random_number)
    weights = np.array(weights) / np.sum(np.array(weights))
    return sorted(weights, reverse=True)


def generate_fastq_file_worker(df,
                               read_id_suffix,
                               read_length_avg,
                               read_length_stdev,
                               base_phred_score_avg,
                               base_phred_score_stdev,
                               shared_list):
    print(len(df), 'items received in dataframe (worker)')
    for index, row in df.iterrows():
        # Simulate reads
        sequence = str(row['Transcript_Sequence'])
        num_reads = int(row['Transcript_Reads'])

        if num_reads == 0:
            continue

        if len(sequence) <= read_length_avg:
            max_start = 0
        else:
            max_start = len(sequence) - read_length_avg

        simulated_read_ids = []
        simulated_reads = []
        simulated_base_scores = []

        for i in range(0, num_reads):
            # Randomly generate a read ID
            # Example: @m64220e_210621_231818/66840/ccs
            random_read_id = read_id_suffix + str(i + 1) + '/ccs'
            simulated_read_ids.append(random_read_id)

            # Randomly select a start site to read
            random_start = random.randint(0, max_start)

            # Randomly generate a read length
            random_read_length = math.ceil(np.random.normal(read_length_avg, read_length_stdev))
            random_end = random_start + random_read_length - 1

            # Randomly generate read
            random_read = sequence[random_start:random_end]
            simulated_reads.append(random_read)

            # Randomly generate base phred scores
            random_base_phred_scores = np.array([int(np.random.normal(base_phred_score_avg, base_phred_score_stdev)) for j in range(0, len(random_read))])
            random_base_phred_scores[random_base_phred_scores > 126] = 126
            random_base_phred_scores = list(random_base_phred_scores)
            random_base_phred_scores_chars = ''.join([chr(j) for j in random_base_phred_scores])
            simulated_base_scores.append(random_base_phred_scores_chars)

        for i in range(0, len(simulated_read_ids)):
            shared_list.append((simulated_read_ids[i],
                                simulated_reads[i],
                                simulated_base_scores[i]))


def generate_fastq_file(df,
                        candidate_genes,
                        total_num_reads,
                        read_length_avg,
                        read_length_stdev,
                        base_phred_score_avg,
                        base_phred_score_stdev,
                        output_fastq_file,
                        output_transcript_proportions_file,
                        num_cores):
    df = df.loc[df['Gene_Symbol'].isin(candidate_genes),:]
    total_transcript_bases = df['Transcript_Length'].sum()
    print(total_transcript_bases, 'total bases in', len(df), 'transcripts')

    df_transcript_proportions = pd.DataFrame()
    counter = 1
    for curr_candidate_gene in candidate_genes:
        # Fetch all transcripts
        df_matched = df.loc[
            df['Gene_Symbol'] == curr_candidate_gene,:
        ]

        # Calculate number of reads to allocate to current gene
        curr_transcript_bases = df_matched['Transcript_Length'].sum()
        curr_gene_num_reads = int((curr_transcript_bases / total_transcript_bases) * total_num_reads)
        print(str(counter) + '/' + str(len(candidate_genes)), curr_candidate_gene, ':', curr_gene_num_reads, 'reads allocated')
        counter += 1

        # Calculate number of reads to allocate to each transcript
        df_matched_ref = df_matched.loc[df_matched['Transcript_Type'] == 'protein_coding',:]
        df_matched_ref = df_matched_ref.sort_values(by=['Transcript_Length'], ascending=False)
        df_matched_var = df_matched.loc[df_matched['Transcript_Type'] == 'protein_coding_variant',:]
        random_weights = randomly_generate_weights(n=len(df_matched_ref) + 1)
        transcript_num_reads = [int(curr_gene_num_reads * i) for i in random_weights]

        # Append variant reads
        df_matched_var_temp = df_matched_var.copy()
        df_matched_var_temp.reset_index(inplace=True)
        df_matched_var_temp.loc[0, 'Transcript_Reads'] = transcript_num_reads[0]
        df_transcript_proportions = pd.concat([df_transcript_proportions, df_matched_var_temp.iloc[[0]]])

        # Append reference reads
        for i in range(0, len(df_matched_ref)):
            df_matched_ref_temp = df_matched_ref.iloc[[i]].copy()
            df_matched_ref_temp.reset_index(inplace=True)
            df_matched_ref_temp.loc[0, 'Transcript_Reads'] = transcript_num_reads[i + 1]
            df_transcript_proportions = pd.concat([df_transcript_proportions, df_matched_ref_temp])
    df_transcript_proportions.to_csv(output_transcript_proportions_file,
                                     sep='\t', index=False)

    # Parallelize simulation
    num_1 = str(random.randint(10000, 99999)) + 'e'
    num_2 = random.randint(100000, 200000)
    num_3 = num_2 + random.randint(10000, 99999)
    read_id_suffix = '@m' + num_1 + '_' + str(num_2) + '_' + str(num_3) + '/'

    list_df = np.array_split(df_transcript_proportions, num_cores)
    print(len(list_df), 'dataframes to parallelize')
    pool = mp.Pool(processes=num_cores)
    manager = mp.Manager()
    shared_list = manager.list()
    for df_curr in list_df:
        pool.apply_async(generate_fastq_file_worker,
                         args=[
                             df_curr,
                             read_id_suffix,
                             read_length_avg,
                             read_length_stdev,
                             base_phred_score_avg,
                             base_phred_score_stdev,
                             shared_list
                         ])
    pool.close()
    pool.join()

    # Write reads to FASTQ file
    print(len(shared_list), 'items in shared list')
    with open(output_fastq_file, 'w') as f:
        for i in shared_list:
            curr_read_id = i[0]
            curr_read_sequence = i[1]
            curr_base_quality = i[2]
            f.write(curr_read_id + '\n')
            f.write(curr_read_sequence + '\n')
            f.write('+\n')
            f.write(curr_base_quality + '\n')


if __name__ == "__main__":
    args = parse_args()

    # Step 1. Read GENCODE transcript sequences
    df_gencode = read_gencode_fasta_file(fasta_file=args.gencode_transcripts_fasta_file)

    # Step 3. Read list of oncogenes and tumor suppressor genes
    df_bailey_s1 = read_bailey_oncogenes_and_tsgs_list_file(excel_file=args.bailey_et_al_cell_2018_s1_excel_file)

    # Step 4. Identify oncogenes and tumor suppressors in GENCODE
    candidate_genes = set.intersection(set(df_gencode['Gene_Symbol'].unique()),
                                       set(df_bailey_s1['Gene'].unique()))
    print(len(candidate_genes), 'oncogenes and TSGs found in GENCODE')

    # Step 5. Generate variant transcript sequences
    df_gencode_variants = generate_variant_transcript_sequences(
        candidate_genes=candidate_genes,
        df_gencode=df_gencode,
        deletion_size_average=args.deletion_size_average,
        deletion_size_stdev=args.deletion_size_stdev,
        insertion_size_average=args.insertion_size_average,
        insertion_size_stdev=args.insertion_size_stdev,
        duplication_size_average=args.duplication_size_average,
        duplication_size_stdev=args.duplication_size_stdev
    )
    df_gencode_variants.to_csv(
        args.output_dir + '/' + args.sample_id + '.tsv',
        sep='\t', index=False
    )

    # Step 6. Write to FASTA file
    df = pd.concat([df_gencode, df_gencode_variants])
    write_transcript_sequences_to_fasta_file(
        df=df,
        output_fasta_file=args.output_dir + '/' + args.sample_id + '.fa'
    )

    # Step 7. Generate FASTQ file
    generate_fastq_file(
        df=df,
        candidate_genes=candidate_genes,
        total_num_reads=args.total_num_reads,
        read_length_avg=args.read_length_average,
        read_length_stdev=args.read_length_stdev,
        base_phred_score_avg=args.base_phred_score_average,
        base_phred_score_stdev=args.base_phred_score_stdev,
        output_fastq_file=args.output_dir + '/' + args.sample_id + '.fastq',
        output_transcript_proportions_file=args.output_dir + '/' + args.sample_id + '_ground_truth.tsv',
        num_cores=args.num_cores
    )
