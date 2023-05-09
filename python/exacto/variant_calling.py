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


import multiprocessing as mp
import numpy as np
import pandas as pd
import pysam
import re
from typing import Dict
from .logging import get_logger
from .constants import VariantTypes
from .variant_call import VariantCall


logger = get_logger(__name__)


def identify_rna_variant_in_cs_tag(
        chromosome: str,
        start: int,
        cs_tag: str,
        read_id: str
) -> Dict:
    """
    Calls RNA variants based on start position and CS tag.

    Parameters
    ----------
    start_pos               :   Start position of alignment.
    cs_tag                  :   CS tag.

    Returns
    -------
    rna_variant_calls_dict  :   Dictionary with the following keys:
                                'chromosome_1',
                                'position_1',
                                'chromosome_2',
                                'position_2',
                                'variant_type',
                                'reference_allele',
                                'alternate_allele',
                                'variant_sequence',
                                'variant_size'
                                'alternate_allele_read_ids'
    """
    rna_variant_calls_dict = {
        'chromosome_1': [],
        'position_1': [],
        'chromosome_2': [],
        'position_2': [],
        'variant_type': [],
        'reference_allele': [],
        'alternate_allele': [],
        'variant_sequence': [],
        'variant_size': [],
        'alternate_allele_read_ids': []
    }
    cs_tag_delimited = re.split(pattern='([=:*+~-])', string=cs_tag)
    curr_pos = start - 1
    curr_tag = ''
    for curr_element in cs_tag_delimited:
        if curr_element == '':
            continue
        elif curr_element in ['=', ':', '*', '+', '~', '-']:
            curr_tag = curr_element
            continue
        else:
            if curr_tag == ':':
                curr_pos += int(curr_element)
            elif curr_tag == '*':
                curr_pos += 1
                rna_variant_calls_dict['chromosome_1'].append(chromosome)
                rna_variant_calls_dict['position_1'].append(curr_pos)
                rna_variant_calls_dict['chromosome_2'].append(chromosome)
                rna_variant_calls_dict['position_2'].append(curr_pos)
                rna_variant_calls_dict['variant_type'].append(VariantTypes.SINGLE_NUCLEOTIDE_VARIANT)
                rna_variant_calls_dict['reference_allele'].append(curr_element[0])
                rna_variant_calls_dict['alternate_allele'].append(curr_element[1])
                rna_variant_calls_dict['variant_sequence'].append(curr_element[1])
                rna_variant_calls_dict['variant_size'].append(1)
                rna_variant_calls_dict['alternate_allele_read_ids'].append(read_id)
            elif curr_tag == '-':
                rna_variant_calls_dict['chromosome_1'].append(chromosome)
                rna_variant_calls_dict['position_1'].append(curr_pos + 1)
                rna_variant_calls_dict['chromosome_2'].append(chromosome)
                rna_variant_calls_dict['position_2'].append(curr_pos + len(curr_element))
                rna_variant_calls_dict['variant_type'].append(VariantTypes.DELETION)
                rna_variant_calls_dict['reference_allele'].append(curr_element)
                rna_variant_calls_dict['alternate_allele'].append('')
                rna_variant_calls_dict['variant_sequence'].append(curr_element)
                rna_variant_calls_dict['variant_size'].append(len(curr_element))
                rna_variant_calls_dict['alternate_allele_read_ids'].append(read_id)
                curr_pos += len(curr_element)
            elif curr_tag == '+':
                rna_variant_calls_dict['chromosome_1'].append(chromosome)
                rna_variant_calls_dict['position_1'].append(curr_pos)
                rna_variant_calls_dict['chromosome_2'].append(chromosome)
                rna_variant_calls_dict['position_2'].append(curr_pos)
                rna_variant_calls_dict['variant_type'].append(VariantTypes.INSERTION)
                rna_variant_calls_dict['reference_allele'].append('')
                rna_variant_calls_dict['alternate_allele'].append(curr_element)
                rna_variant_calls_dict['variant_sequence'].append(curr_element)
                rna_variant_calls_dict['variant_size'].append(len(curr_element))
                rna_variant_calls_dict['alternate_allele_read_ids'].append(read_id)
            elif curr_tag == '~':
                intron_length = int(re.findall(r'[0-9]+', curr_element)[0])
                curr_pos += intron_length
            else:
                logger.warning("Unknown CS element: %s" % curr_tag)
    return rna_variant_calls_dict

def call_rna_variants_worker(df_reads: pd.DataFrame) -> Dict:
    """
    Worker function for calling variants based on CS tags.

    Parameters
    ----------
    df_reads                :   DataFrame with the following columns:
                                'read_id'
                                'chrom'
                                'start'
                                'cs_tag'

    rna_variant_calls_dict  :   Dictionary with the following keys:
                                'chromosome_1'
                                'position_1'
                                'chromosome_2'
                                'position_2'
                                'variant_type'
                                'reference_allele'
                                'alternate_allele'
                                'variant_sequence'
                                'variant_size'
                                'alternate_allele_read_ids'
    """
    rna_variant_calls_dict = {
        'chromosome_1': [],
        'position_1': [],
        'chromosome_2': [],
        'position_2': [],
        'variant_type': [],
        'reference_allele': [],
        'alternate_allele': [],
        'variant_sequence': [],
        'variant_size': [],
        'alternate_allele_read_ids': []
    }
    for row in df_reads.to_dict('records'):
        calls = identify_rna_variant_in_cs_tag(
            chromosome=row['chromm'],
            start=row['start'],
            cs_tag=row['cs_tag'],
            read_id=row['read_id']
        )
        for key, value in calls.items():
            rna_variant_calls_dict[key] += value
    return rna_variant_calls_dict

def call_rna_variants(
        bam_file: pysam.AlignmentFile,
        num_processes
    ) -> pd.DataFrame:
    """
    Calls RNA variants in a BAM file.

    Parameters
    ----------
    bam_file            :   pysam.AlignmentFile object.
    num_processes       :   Number of processes.

    Returns
    -------
    df_variant_calls    :   DataFrame with the following columns:
                            'chromosome_1'
                            'position_1'
                            'chromosome_2'
                            'position_2'
                            'variant_type'
                            'reference_allele'
                            'alternate_allele'
                            'variant_sequence'
                            'variant_size'
                            'alternate_allele_read_ids'
    """
    # Step 1. Get all reads into a DataFrame
    reads_data_dict = {
        'read_id': [],
        'chrom': [],
        'start': [],
        'cs_tag': []
    }
    for read in bam_file.fetch():
        curr_chrom = read.reference_name
        curr_read_id = read.qname
        curr_start_pos = read.pos
        curr_cs_tag = read.get_tag(tag='cs')
        reads_data_dict['read_id'].append(curr_read_id)
        reads_data_dict['chrom'].append(curr_chrom)
        reads_data_dict['start'].append(curr_start_pos)
        reads_data_dict['cs_tag'].append(curr_cs_tag)
    df_reads = pd.DataFrame(reads_data_dict)
    logger.info("%i reads in total." % len(df_reads))

    # Step 2. Multiprocess calling variants from reads
    list_df_reads = np.array_split(df_reads, num_processes)
    pool = mp.Pool(processes=num_processes)
    async_results = [pool.apply_async(call_rna_variants_worker, args=(df_curr_reads,)) for df_curr_reads in list_df_reads]
    pool.close()
    pool.join()
    calls_list = [ar.get() for ar in async_results]

    # Step 3. Convert tuples into a DataFrame
    all_calls = {
        'chromosome_1': [],
        'position_1': [],
        'chromosome_2': [],
        'position_2': [],
        'variant_type': [],
        'reference_allele': [],
        'alternate_allele': [],
        'variant_sequence': [],
        'variant_size': [],
        'alternate_allele_read_ids': []
    }
    for calls in calls_list:
        for key, value in calls.items():
            all_calls[key] += value
    df_variants = pd.DataFrame(all_calls)
    return df_variants
