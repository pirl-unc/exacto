#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
reading files.

Author: Jin Seok (Andy) Lee

Last updated date: July 20, 2022
"""


import pandas as pd
import numpy as np
import logging


def read_gencode_gtf_file(gencode_gtf_file: str) -> pd.DataFrame:
    """
    Reads a GENCODE GTF file and returns a DataFrame.

    Args
    ----
    gencode_gtf_file    :   GTF file path.

    Returns
    -------
    df_gtf              :   DataFrame with the following columns:
                            'gene_id',
                            'gene_name',
                            'gene_type',
                            'transcript_id',
                            'transcript_name',
                            'transcript_type',
                            'transcript_support_level',
                            'exon_id',
                            'exon_number',
                            'chrom',
                            'start',
                            'end',
                            'strand',
                            'level'
    """
    logging.info('Started reading GENCODE GTF file.')

    data = {
        'gene_id': [],
        'gene_name': [],
        'gene_type': [],
        'transcript_id': [],
        'transcript_name': [],
        'transcript_type': [],
        'transcript_support_level': [],
        'exon_id': [],
        'exon_number': [],
        'chrom': [],
        'start': [],
        'end': [],
        'strand': [],
        'level': []
    }

    with open(gencode_gtf_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0:2] == '##':
                continue
            elements = line.split('\t')
            if elements[2] == 'exon':
                curr_chrom = elements[0]
                curr_start = int(elements[3])
                curr_end = int(elements[4])
                curr_strand = str(elements[6])
                curr_metadata = str(elements[8]).split(';')
                curr_metadata_dict = {}
                for curr_metadata_elements in curr_metadata:
                    if curr_metadata_elements == '':
                        continue
                    if curr_metadata_elements[0] == ' ':
                        curr_metadata_elements = curr_metadata_elements[1:]
                    curr_metadata_elements_ = curr_metadata_elements.split(' ')
                    curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')

                data['gene_id'].append(curr_metadata_dict['gene_id'])
                data['gene_name'].append(curr_metadata_dict['gene_name'])
                data['gene_type'].append(curr_metadata_dict['gene_type'])
                data['transcript_id'].append(curr_metadata_dict['transcript_id'])
                data['transcript_name'].append(curr_metadata_dict['transcript_name'])
                data['transcript_type'].append(curr_metadata_dict['transcript_type'])
                data['transcript_support_level'].append(curr_metadata_dict.get('transcript_support_level', np.nan))
                data['exon_id'].append(curr_metadata_dict['exon_id'])
                data['exon_number'].append(int(curr_metadata_dict['exon_number']))
                data['level'].append(curr_metadata_dict.get('level', np.nan))
                data['chrom'].append(curr_chrom)
                data['start'].append(curr_start)
                data['end'].append(curr_end)
                data['strand'].append(curr_strand)

    df = pd.DataFrame(data)
    logging.info('Finished reading GENCODE GTF file.')
    return df


def read_gencode_refseq_file(gencode_refseq_metadata_file: str) -> pd.DataFrame:
    """
    Reads a GENCODE RefSeq metadata file and returns a DataFrame

    Args
    ----
    ref
    :param refseq_metadata_file:
    :return:
    """
    logging.info('Started reading GENCODE RefSeq file.')
    df_refseq = pd.read_csv(gencode_refseq_metadata_file,
                            sep='\t',
                            header=None)
    df_refseq.columns = ['ensembl_transcript_id',
                         'refseq_transcript_id',
                         'refseq_protein_id']
    logging.info('Finished reading GENCODE RefSeq file.')
    return df_refseq


def read_gencode_transcripts_fasta_file(gencode_fasta_file: str) -> pd.DataFrame:
    logging.info('Started reading GENCODE transcripts FASTA file.')

    data = {
        'ensembl_transcript_id': [],
        'ensembl_gene_id': [],
        'havana_gene_id': [],
        'havana_transcript_id': [],
        'gene_name_versioned': [],
        'gene_name': [],
        'transcript_length': [],
        'utr5_start': [],
        'utr5_end': [],
        'cds_start': [],
        'cds_end': [],
        'utr3_start': [],
        'utr3_end': [],
        'sequence': []
    }
    first = True
    curr_sequence = ''

    with open(gencode_fasta_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0] == '>':
                # Append previously stored information
                if first:
                    first = False
                else:
                    data['ensembl_transcript_id'].append(curr_ensembl_transcript_id)
                    data['ensembl_gene_id'].append(curr_ensembl_gene_id)
                    data['havana_gene_id'].append(curr_havana_gene_id)
                    data['havana_transcript_id'].append(curr_havana_transcript_id)
                    data['gene_name_versioned'].append(curr_gene_name_versioned)
                    data['gene_name'].append(curr_gene_name)
                    data['transcript_length'].append(curr_transcript_length)
                    data['utr5_start'].append(curr_utr5_start)
                    data['utr5_end'].append(curr_utr5_end)
                    data['cds_start'].append(curr_cds_start)
                    data['cds_end'].append(curr_cds_end)
                    data['utr3_start'].append(curr_utr3_start)
                    data['utr3_end'].append(curr_utr3_end)
                    data['sequence'].append(curr_sequence)
                    curr_sequence = ''

                line = line[1:]
                line_elements = line.split('|')
                curr_ensembl_transcript_id = line_elements[0]
                curr_ensembl_gene_id = line_elements[1]
                curr_havana_gene_id = line_elements[2]
                curr_havana_transcript_id = line_elements[3]
                curr_gene_name_versioned = line_elements[4]
                curr_gene_name = line_elements[5]
                curr_transcript_length = int(line_elements[6])
                curr_utr5_start = -1
                curr_utr5_end = -1
                curr_cds_start = -1
                curr_cds_end = -1
                curr_utr3_start = -1
                curr_utr3_end = -1
                for element in line_elements:
                    if 'UTR5:' in element:
                        curr_utr5_start = int(element.split(':')[1].split('-')[0])
                        curr_utr5_end = int(element.split(':')[1].split('-')[1])
                    if 'CDS:' in element:
                        curr_cds_start = int(element.split(':')[1].split('-')[0])
                        curr_cds_end = int(element.split(':')[1].split('-')[1])
                    if 'UTR3:' in element:
                        curr_utr3_start = int(element.split(':')[1].split('-')[0])
                        curr_utr3_end = int(element.split(':')[1].split('-')[1])
            else:
                curr_sequence += line

        # Store last transcript
        data['ensembl_transcript_id'].append(curr_ensembl_transcript_id)
        data['ensembl_gene_id'].append(curr_ensembl_gene_id)
        data['havana_gene_id'].append(curr_havana_gene_id)
        data['havana_transcript_id'].append(curr_havana_transcript_id)
        data['gene_name_versioned'].append(curr_gene_name_versioned)
        data['gene_name'].append(curr_gene_name)
        data['transcript_length'].append(curr_transcript_length)
        data['utr5_start'].append(curr_utr5_start)
        data['utr5_end'].append(curr_utr5_end)
        data['cds_start'].append(curr_cds_start)
        data['cds_end'].append(curr_cds_end)
        data['utr3_start'].append(curr_utr3_start)
        data['utr3_end'].append(curr_utr3_end)
        data['sequence'].append(curr_sequence)

    df = pd.DataFrame(data)
    logging.info('Finished reading GENCODE transcripts FASTA file.')
    return df

