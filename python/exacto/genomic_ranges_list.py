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
The purpose of this python3 script is to implement the GenomicRangesList dataclass.
"""


import pandas as pd
from bisect import bisect_left, bisect_right, insort
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Type, Dict
from .genomic_range import GenomicRange
from .logging import get_logger


logger = get_logger(__name__)


@dataclass(frozen=True)
class GenomicRangesList:
    genomic_ranges = defaultdict(list)

    @property
    def size(self):
        size = 0
        for _ in self.genomic_ranges.values():
            size += 1
        return size

    @staticmethod
    def load_dataframe(df) -> Type["GenomicRangesList"]:
        """
        Reads a DataFrame and returns a GenomicRangesList object.

        Parameters
        ----------
        df                  :   DataFrame.

        Returns
        -------
        genomic_ranges_list :   GenomicRangesList object.
        """
        genomic_ranges_list = GenomicRangesList()
        for index, row in df.iterrows():
            chromosome = str(row['chromosome'])
            start = int(row['start'])
            end = int(row['end'])
            genomic_range = GenomicRange(
                chromosome=chromosome,
                start=start,
                end=end
            )
            genomic_ranges_list.add_genomic_range(genomic_range=genomic_range)
        return genomic_ranges_list

    @staticmethod
    def read_tsv_file(tsv_file):
        """
        Reads a TSV file and returns a GenomicRangesList object.

        Parameters
        ----------
        tsv_file            :   TSV file.

        Returns
        -------
        genomic_ranges_list :   GenomicRangesList object.
        """
        df = pd.read_csv(tsv_file, sep='\t', low_memory=False, memory_map=True)
        return GenomicRangesList.load_dataframe(df=df)

    def add_genomic_range(self, genomic_range: GenomicRange):
        """
        Adds a GenomicRange object.

        Parameters
        ----------
        genomic_range   :   GenomicRange object.
        """
        insort(self.genomic_ranges[genomic_range.chromosome], genomic_range)

    def find_overlaps(self, chromosome, start, end) -> List[GenomicRange]:
        """
        Finds GenomicRange objects that overlap with query position.

        Parameters
        ----------
        chromosome      :   Chromosome.
        start           :   Start position.
        end             :   End position.

        Returns
        -------
        genomic_ranges  :   List of GenomicRange objects.
        """
        # Get GenomicRange objects that match the query position
        genomic_ranges = []
        for genomic_range in self.genomic_ranges[chromosome]:
            if genomic_range.overlaps(chromosome=chromosome, start=start, end=end):
                genomic_ranges.append(genomic_range)
        return genomic_ranges

    def remove_genomic_range(self, genomic_range: GenomicRange):
        """
        Removes a GenomicRange object.

        Parameters
        ----------
        genomic_range   :   GenomicRange object.
        """
        for i in range(0, len(self.genomic_ranges[genomic_range.chromosome])):
            if self.genomic_ranges[genomic_range.chromosome][i] == genomic_range:
                self.genomic_ranges[genomic_range.chromosome].remove(i)
                break