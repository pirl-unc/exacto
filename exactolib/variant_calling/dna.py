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


import pandas as pd
import numpy as np
import re
import pysam
import multiprocessing as mp
from multiprocessing import Process, Manager
from typing import Tuple, List


from ..logging import get_logger


logger = get_logger(__name__)


def call_dna_variants_from_cs_tag(start_pos: int, cs_tag: str) -> pd.DataFrame:
    """
    Calls variants based on start position and CS tag.

    Parameters
    ----------
    start_pos   :   Start position.
    cs_tag      :   CS tag.

    Returns
    -------
    df_variants :   DataFrame with the following columns:
                    'pos'
                    'variant_type'
                    'ref'
                    'alt'
                    'sequence'
    """
    cs_tag_delimited = re.split(pattern='([:*+-])', string=cs_tag)
    curr_pos = start_pos
    curr_tag = ''
    data = {
        'pos': [],
        'variant_type': [],
        'ref': [],
        'alt': [],
        'sequence': [],
        'variant_size': []
    }
    for curr_element in cs_tag_delimited:
        if curr_element == '':
            continue
        elif curr_element == ':':
            curr_tag = curr_element
            continue
        elif curr_element == '*':
            curr_tag = curr_element
            continue
        elif curr_element == '-':
            curr_tag = curr_element
            continue
        elif curr_element == '+':
            curr_tag = curr_element
            continue
        else:
            if curr_tag == ':':
                curr_pos += int(curr_element)
            elif curr_tag == '*':
                data['pos'].append(curr_pos)
                data['variant_type'].append('snv')
                data['ref'].append(curr_element[0])
                data['alt'].append(curr_element[1])
                data['sequence'].append(curr_element[1])
                data['variant_size'].append(1)
                curr_pos += 1
            elif curr_tag == '-':
                data['pos'].append(curr_pos)
                data['variant_type'].append('deletion')
                data['ref'].append(curr_element)
                data['alt'].append('')
                data['sequence'].append(curr_element)
                data['variant_size'].append(len(curr_element))
                curr_pos += len(curr_element)
            elif curr_tag == '+':
                data['pos'].append(curr_pos)
                data['variant_type'].append('insertion')
                data['ref'].append('')
                data['alt'].append(curr_element)
                data['sequence'].append(curr_element)
                data['variant_size'].append(len(curr_element))
            else:
                logger.warning("Unknown CS element: %s" % curr_tag)
    df_variants = pd.DataFrame(data)
    return df_variants


def call_dna_variants_worker(df_reads: pd.DataFrame,
                             shared_list: list):
    """
    Worker function for calling variants based on CS tags.

    Parameters
    ----------
    df_reads    :   DataFrame with the following columns:
                    'read_id'
                    'chrom'
                    'start'
                    'cs_tag'

    shared_list :   Shared list (mp.Manager().list() object) to which
                    called variants will be appended.
    """
    for row in df_reads.to_dict('records'):
        curr_read_id = row['read_id']
        curr_chrom = row['chrom']
        curr_start = row['start']
        curr_cs_tag = row['cs_tag']

        cs_tag_delimited = re.split(pattern='([:*+-])', string=curr_cs_tag)
        curr_pos = curr_start
        curr_tag = ''
        for curr_element in cs_tag_delimited:
            if curr_element == '':
                continue
            elif curr_element == ':':
                curr_tag = curr_element
                continue
            elif curr_element == '*':
                curr_tag = curr_element
                continue
            elif curr_element == '-':
                curr_tag = curr_element
                continue
            elif curr_element == '+':
                curr_tag = curr_element
                continue
            else:
                if curr_tag == ':':
                    curr_pos += int(curr_element)
                elif curr_tag == '*':
                    shared_list.append([curr_read_id,       # read_id
                                        curr_chrom,         # chrom
                                        curr_pos,           # pos
                                        'snv',              # variant_type
                                        curr_element[0],    # ref
                                        curr_element[1],    # alt
                                        curr_element[1],    # sequence
                                        1]                  # variant_size
                    )
                    curr_pos += 1
                elif curr_tag == '-':
                    shared_list.append([curr_read_id,       # read_id
                                        curr_chrom,         # chrom
                                        curr_pos,           # pos
                                        'deletion',         # variant_type
                                        curr_element,       # ref
                                        '',                 # alt
                                        curr_element,       # sequence
                                        len(curr_element)]  # variant_size
                    )
                    curr_pos += len(curr_element)
                elif curr_tag == '+':
                    shared_list.append([curr_read_id,       # read_id
                                        curr_chrom,         # chrom
                                        curr_pos,           # pos
                                        'insertion',        # variant_type
                                        '',                 # ref
                                        curr_element,       # alt
                                        curr_element,       # sequence
                                        len(curr_element)]  # variant_size
                    )
                else:
                    logger.warning("Unknown CS element: %s" % curr_tag)


def call_dna_variants(bam_file: pysam.AlignmentFile,
                      target_chromosomes: List[str],
                      num_processes) -> pd.DataFrame:
    """
    Calls DNA variants in a BAM file.

    Parameters
    ----------
    bam_file            :   pysam.AlignmentFile object.
    target_chromosomes  :   Target chromosomes.
    num_processes       :   Number of processes.

    Returns
    -------
    df_variants         :   DataFrame with the following columns:
                            'read_id'
                            'chrom'
                            'pos'
                            'variant_type'
                            'ref'
                            'alt'
                            'sequence'
                            'variant_size'
    """
    # Step 1. Get all reads into a DataFrame
    reads_data = {
        'read_id': [],
        'chrom': [],
        'start': [],
        'cs_tag': []
    }
    for read in bam_file.fetch():
        curr_chrom = read.reference_name

        # Skip if not a target chromosome
        if curr_chrom not in target_chromosomes:
            continue

        curr_read_id = read.qname
        curr_start_pos = read.pos + 1 # BAM file and pysam use 0-based coordinate
        curr_cs_tag = read.get_tag(tag='cs')
        reads_data['read_id'].append(curr_read_id)
        reads_data['chrom'].append(curr_chrom)
        reads_data['start'].append(curr_start_pos)
        reads_data['cs_tag'].append(curr_cs_tag)
    df_reads = pd.DataFrame(reads_data)
    logger.info("%i reads in total." % len(df_reads))

    # Step 2. Multiprocess calling variants from reads
    list_df_reads = np.array_split(df_reads, num_processes)
    pool = mp.Pool(processes=num_processes)
    manager = mp.Manager()
    shared_list = manager.list()
    for df_curr_reads in list_df_reads:
        pool.apply_async(call_dna_variants_worker, args=[df_curr_reads, shared_list])
        logger.info("Launched a variant calling worker process.")
    pool.close()
    pool.join()
    logger.info("All variant calling worker processes have completed.")

    # Step 3. Convert tuples into a DataFrame
    variants_data = {
        'read_id': [],
        'chrom': [],
        'pos': [],
        'variant_type': [],
        'ref': [],
        'alt': [],
        'sequence': [],
        'variant_size': []
    }
    for curr_element in shared_list:
        variants_data['read_id'].append(curr_element[0])
        variants_data['chrom'].append(curr_element[1])
        variants_data['pos'].append(curr_element[2])
        variants_data['variant_type'].append(curr_element[3])
        variants_data['ref'].append(curr_element[4])
        variants_data['alt'].append(curr_element[5])
        variants_data['sequence'].append(curr_element[6])
        variants_data['variant_size'].append(curr_element[7])
    df_variants = pd.DataFrame(variants_data)
    return df_variants



