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
The purpose of this python3 script is to implement the AlignmentMap dataclass.
"""


import multiprocessing as mp
import numpy as np
import pandas as pd
import pysam
import re
from dataclasses import dataclass
from typing import Dict
from .aligned_read import AlignedRead
from .constants import VariantTypes, NucleicAcidTypes
from .logging import get_logger


logger = get_logger(__name__)


@dataclass(frozen=True)
class AlignmentMap:
    bam: pysam.AlignmentFile
    nucleic_acid_type: str

    @property
    def bam_file(self):
        return str(self.bam.filename.decode())

    @property
    def chromosome_sizes(self):
        chrom_sizes = {}  # key = chromosome, value = length
        for curr_chrom in self.bam.references:
            curr_size = self.bam.get_reference_length(curr_chrom)
            chrom_sizes[curr_chrom] = curr_size
        return chrom_sizes

    @property
    def num_reads(self):
        num_reads = 0
        for _ in self.bam.fetch():
            num_reads += 1
        return num_reads

    @staticmethod
    def identify_rna_variants_in_cs_tag(
            chromosome: str,
            start: int,
            cs_tag: str,
            read_id: str
    ) -> Dict:
        """
        Identifies RNA variants based on start position and CS tag.

        Parameters
        ----------
        chromosome              :   Chromosome.
        start                   :   Start position of alignment.
        cs_tag                  :   CS tag.
        read_id                 :   Read ID.

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
                    raise Exception('Unknown CS tag element: %s' % curr_tag)
        return rna_variant_calls_dict

    @staticmethod
    def call_rna_variants_worker(reads):
        """
        Worker function for calling RNA variants.

        Parameters
        ----------
        reads               :   List of instances of 'Read' class.

        Returns
        -------
        variant_calls       :   List of VariantCall instances.
        """
        calls = {
            'chromosome': [],
            'start': [],
            'variant_type': [],
            'ref': [],
            'alt': [],
            'variant_sequence': [],
            'variant_size': [],
            'read_id': []
        }
        idx = 1
        for read in reads:
            # Parse CS tag
            curr_cs_tag_delimited = re.split(pattern='([:*+~-])', string=read.cs_tag)
            curr_pos = read.start
            curr_tag = ''
            for curr_element in curr_cs_tag_delimited:
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
                elif curr_element == '~':
                    curr_tag = curr_element
                    continue
                else:
                    if curr_tag == ':':
                        curr_pos += int(curr_element)
                    elif curr_tag == '*':
                        calls['chromosome'].append(read.chromosome)
                        calls['start'].append(curr_pos)
                        calls['variant_type'].append(VariantTypes.SINGLE_NUCLEOTIDE_VARIANT)
                        calls['ref'].append(curr_element[0])
                        calls['alt'].append(curr_element[1])
                        calls['variant_sequence'].append(curr_element[1])
                        calls['variant_size'].append(1)
                        calls['read_id'].append(read.id)
                        idx += 1
                        curr_pos += 1
                    elif curr_tag == '-':
                        calls['chromosome'].append(read.chromosome)
                        calls['start'].append(curr_pos)
                        calls['variant_type'].append(VariantTypes.DELETION)
                        calls['ref'].append(curr_element)
                        calls['alt'].append('')
                        calls['variant_sequence'].append(curr_element)
                        calls['variant_size'].append(len(curr_element))
                        calls['read_id'].append(read.id)
                        idx += 1
                        curr_pos += len(curr_element)
                    elif curr_tag == '+':
                        calls['chromosome'].append(read.chromosome)
                        calls['start'].append(curr_pos)
                        calls['variant_type'].append(VariantTypes.INSERTION)
                        calls['ref'].append('')
                        calls['alt'].append(curr_element)
                        calls['variant_sequence'].append(curr_element)
                        calls['variant_size'].append(len(curr_element))
                        calls['read_id'].append(read.id)
                        idx += 1
                    elif curr_tag == '~':
                        continue
                    else:
                        logger.warning("Unknown CS tag: %s" % curr_tag)
        return calls

    def call_variants(self, num_processes: int) -> pd.DataFrame:
        """
        Calls variants in self.bam.

        Parameters
        ----------
        num_processes       :   Number of processes.

        Returns
        -------
        df_calls            :   pd.DataFrame
        """
        # Step 1. Fetch all reads
        logger.info('Started fetching all reads in BAM file.')
        reads = []
        for read in self.bam.fetch():
            read = AlignedRead(
                id=read.qname,
                chromosome=read.reference_name,
                start=read.pos + 1, # BAM file and pysam use 0-based coordinate
                cs_tag=read.get_tag(tag='cs')
            )
            reads.append(read)
        logger.info('Finished fetching all reads in BAM file.')

        # Step 2. Multiprocess variant calling
        logger.info('Started calling variants.')
        reads_list = np.array_split(reads, num_processes)
        pool = mp.Pool(processes=num_processes)
        if self.nucleic_acid_type == NucleicAcidTypes.RNA:
            worker_fn = AlignmentMap.call_rna_variants_worker
        async_results = [pool.apply_async(worker_fn, args=(reads_list[i], i)) for i in range(0, len(reads_list))]
        pool.close()
        pool.join()
        calls_list = [ar.get() for ar in async_results]
        logger.info('Finished calling variants.')

        # Step 3. Merge calls
        logger.info('Started merging variant calls into a list.')
        calls_merged = {
            'chromosome': [],
            'start': [],
            'variant_type': [],
            'ref': [],
            'alt': [],
            'variant_sequence': [],
            'variant_size': [],
            'read_id': []
        }
        for calls in calls_list:
            calls_merged = {key: value + calls[key] for key, value in calls_merged.items()}
        df_calls = pd.DataFrame(calls_merged)
        df_calls['id'] = df_calls['chromosome'] + '_' + \
                         df_calls['start'].astype(str) + '_' + \
                         df_calls['variant_type'] + '_' + \
                         df_calls['variant_sequence']
        df_read_ids = df_calls.groupby('id').agg({'read_id': ','.join})
        df_read_ids_count = df_calls.groupby('id').size().to_frame('size')
        df_calls = df_calls.loc[:, ['chromosome',
                                    'start',
                                    'variant_type',
                                    'ref',
                                    'alt',
                                    'variant_sequence',
                                    'variant_size',
                                    'id']]
        df_calls = df_calls.join(df_read_ids, how='outer', on='id')
        df_calls = df_calls.join(df_read_ids_count, how='outer', on='id')
        df_calls['read_ids'] = df_calls['read_id']
        df_calls['read_ids_count'] = df_calls['size']
        df_calls = df_calls.loc[:,['chromosome',
                                   'start',
                                   'variant_type',
                                   'ref',
                                   'alt',
                                   'variant_sequence',
                                   'variant_size',
                                   'read_ids',
                                   'read_ids_count']]
        df_calls.drop_duplicates(inplace=True)
        logger.info('Finished merging variant calls into a list.')
        return df_calls
