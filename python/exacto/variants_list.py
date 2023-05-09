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
The purpose of this python3 script is to implement the VariantsList dataclass.
"""


import pandas as pd
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass, field
from bisect import bisect_left, bisect_right, insort
from typing import List, Type, Dict, Optional, Tuple
from .common import get_typed_value, get_variant_calling_method_attr_types, overlaps_any
from .constants import *
from .default_parameters import *
from .genomic_ranges_list import GenomicRangesList
from .logging import get_logger
from .nucleotide_sequence import NucleotideSequence
from .variant import Variant
from .variant_annotation import VariantAnnotation
from .variant_call import VariantCall
from .variant_filter import VariantFilter


logger = get_logger(__name__)


@dataclass(frozen=True)
class VariantsList:
    variants: List[Variant] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.variants)

    @property
    def variant_call_ids(self) -> List[str]:
        variant_call_ids = []
        for variant in self.variants:
            for variant_call in variant.variant_calls:
                variant_call_ids.append(variant_call.id)
        return variant_call_ids

    @property
    def variant_ids(self) -> List[str]:
        return [i.id for i in self.variants]

    @staticmethod
    def convert_row_value(value, default_value, type):
        """
        Converts a row value.

        Parameters
        ----------
        value           :   Value.
        default_value   :   Default value.
        type            :   Desired value type.

        Returns
        -------
        value           :   Type converted value.
        """
        if pd.isna(value):
            return default_value
        else:
            try:
                if type == str:
                    value = str(value)
                if type == int:
                    value = int(value)
                if type == float:
                    value = float(value)
                if type == bool:
                    value = bool(value)
                return value
            except:
                return default_value

    @staticmethod
    def load_dataframe(df) -> Type["VariantsList"]:
        """
        Reads a DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df              :   DataFrame.

        Returns
        -------
        variants_list   :   An instance of the VariantsList class.
        """
        variants_list = VariantsList()
        for index, row in df.iterrows():
            variant_id = VariantsList.convert_row_value(value=row['variant_id'], default_value=None, type=str)
            variant_call_id = VariantsList.convert_row_value(value=row['variant_call_id'], default_value=None, type=str)
            source_id = VariantsList.convert_row_value(value=row['source_id'], default_value=None, type=str)
            sample_id = VariantsList.convert_row_value(value=row['sample_id'], default_value=None, type=str)
            nucleic_acid = VariantsList.convert_row_value(value=row['nucleic_acid'], default_value=None, type=str)
            variant_calling_method = VariantsList.convert_row_value(value=row['variant_calling_method'], default_value=None, type=str)
            sequencing_platform = VariantsList.convert_row_value(value=row['sequencing_platform'], default_value=None, type=str)
            chromosome_1 = VariantsList.convert_row_value(value=row['chromosome_1'], default_value=None, type=str)
            position_1 = VariantsList.convert_row_value(value=row['position_1'], default_value=None, type=int)
            chromosome_2 = VariantsList.convert_row_value(value=row['chromosome_2'], default_value=None, type=str)
            position_2 = VariantsList.convert_row_value(value=row['position_2'], default_value=None, type=int)
            reference_allele = VariantsList.convert_row_value(value=row['reference_allele'], default_value=None, type=str)
            alternate_allele = VariantsList.convert_row_value(value=row['alternate_allele'], default_value=None, type=str)
            filter = VariantsList.convert_row_value(value=row['filter'], default_value=None, type=str)
            quality_score = VariantsList.convert_row_value(value=row['quality_score'], default_value=None, type=float)
            precise = VariantsList.convert_row_value(value=row['precise'], default_value=None, type=bool)
            variant_type = VariantsList.convert_row_value(value=row['variant_type'], default_value=None, type=str)
            variant_subtype = VariantsList.convert_row_value(value=row['variant_subtype'], default_value=None, type=str)
            variant_size = VariantsList.convert_row_value(value=row['variant_size'], default_value=None, type=int)
            total_read_count = VariantsList.convert_row_value(value=row['total_read_count'], default_value=None, type=int)
            reference_allele_read_count = VariantsList.convert_row_value(value=row['reference_allele_read_count'], default_value=None, type=int)
            alternate_allele_read_count = VariantsList.convert_row_value(value=row['alternate_allele_read_count'], default_value=None, type=int)
            alternate_allele_fraction = VariantsList.convert_row_value(value=row['alternate_allele_fraction'], default_value=None, type=float)
            alternate_allele_softclip_direction = VariantsList.convert_row_value(value=row['alternate_allele_softclip_direction'], default_value=None, type=str)
            variant_sequences_ = VariantsList.convert_row_value(value=row['variant_sequences'], default_value='', type=str)
            alternate_allele_read_ids = VariantsList.convert_row_value(value=row['alternate_allele_read_ids'], default_value='', type=str)

            # Variant sequences
            variant_sequences = []
            if variant_sequences_ != '':
                for seq in variant_sequences_.split(';'):
                    variant_sequences.append(NucleotideSequence(sequence=seq))

            # Alternate allele read IDs
            if alternate_allele_read_ids != '':
                alternate_allele_read_ids = alternate_allele_read_ids.split(';')

            # Tool attributes
            tool_attributes = {}
            if not pd.isna(row['tool_attributes']):
                curr_tool_attr_types = get_variant_calling_method_attr_types(variant_calling_method=variant_calling_method)
                for curr_attr in row['tool_attributes'].split(';'):
                    curr_attr_key = curr_attr.split('=')[0]
                    curr_attr_value = curr_attr.split('=')[1]
                    curr_attr_value = get_typed_value(
                        value=curr_attr_value,
                        default_value=None,
                        type=curr_tool_attr_types[curr_attr_key]
                    )
                    if curr_attr_value is not None:
                        tool_attributes[curr_attr_key] = curr_attr_value

            # Annotations
            # todo handle reading annotations
            # if row['pos_1_annotation_chrom'] != '':
            #     pos_1_annotation_counts = len(row['pos_1_annotation_chrom'].split(';'))
            #     for idx in range(0, pos_1_annotation_counts):
            #         pos_1_annotation_chrom = row['pos_1_annotation_chrom'].split(';')[idx]
            #         pos_1_annotation_chrom = row['pos_1_annotation_chrom'].split(';')[idx]
            #         pos_1_annotation_chrom = row['pos_1_annotation_chrom'].split(';')[idx]
            # variant_annotation = VariantAnnotation()
            # self.__convert_tsv_file_element(value=row['pos_1_annotation_chrom'], nested=True, )

            variant_call = VariantCall(
                id=variant_call_id,
                source_id=source_id,
                sample_id=sample_id,
                nucleic_acid=nucleic_acid,
                variant_calling_method=variant_calling_method,
                sequencing_platform=sequencing_platform,
                chromosome_1=chromosome_1,
                position_1=position_1,
                chromosome_2=chromosome_2,
                position_2=position_2,
                reference_allele=reference_allele,
                alternate_allele=alternate_allele,
                filter=filter,
                quality_score=quality_score,
                precise=precise,
                variant_type=variant_type,
                variant_subtype=variant_subtype,
                variant_size=variant_size,
                variant_sequences=variant_sequences,
                total_read_count=total_read_count,
                reference_allele_read_count=reference_allele_read_count,
                alternate_allele_read_count=alternate_allele_read_count,
                alternate_allele_fraction=alternate_allele_fraction,
                alternate_allele_read_ids=alternate_allele_read_ids,
                alternate_allele_softclip_direction=alternate_allele_softclip_direction,
                tool_attributes=tool_attributes
            )
            variants_list.add_variant_call_to_variant(
                variant_call=variant_call,
                variant_id=variant_id
            )
        logger.info("Loaded %i variants and %i variant calls" % (variants_list.size, len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def merge(variants_lists: List[Type["VariantsList"]],
              max_neighbor_distance: int) -> Type["VariantsList"]:
        """
        Merges a list of VariantsList objects and returns one VariantsList object.

        Parameters
        ----------
        variants_lists                  :   List of VariantsList objects.
        max_neighbor_distance           :   Maximum neighbor distance.
                                            This value is used to decide if
                                            a VariantCall should be appended to an existing Variant.
                                            If there exists a VariantCall in a given Variant
                                            where the distances to both position_1 and position_2
                                            are equal to or less than max_neighbor_distance
                                            to position_1 and position_2 of the specified VariantCall,
                                            respectively, then the specified VariantCall is
                                            appended to the Variant. If such VariantCall
                                            is not identified, then a new Variant is constructed
                                            and added to self.variants.

        Returns
        -------
        variants_list                   :   VariantsList object.
        """
        # Step 1. Create an empty VariantsList
        variants_list_merged = VariantsList()

        # Step 2. Iterate through each VariantsList object and append VariantCall objects.
        for variants_list in variants_lists:
            for variant in variants_list.variants:
                for variant_call in variant.variant_calls:
                    variants_list_merged.add_variant_call(
                        variant_call=variant_call,
                        max_neighbor_distance=max_neighbor_distance
                    )
        return variants_list_merged

    @staticmethod
    def read_tsv_file(tsv_file: str) -> Type["VariantsList"]:
        """
        Reads a TSV file and returns an instance of the VariantsList class.

        Parameters
        ----------
        tsv_file        :   TSV file.

        Returns
        -------
        variants_list   :   An instance of the VariantsList class.
        """
        df = pd.read_csv(tsv_file, sep='\t')
        return VariantsList.load_dataframe(df=df)

    def add_variant(self, variant: Variant):
        """
        Adds a Variant object.

        Parameters
        ----------
        variant                         :   Variant object.
        """
        insort(self.variants, variant)

    def add_variant_call(self, variant_call: VariantCall, max_neighbor_distance: int):
        """
        Adds a VariantCall object.

        Parameters
        ----------
        variant_call                    :   VariantCall object.
        max_neighbor_distance           :   Maximum neighbor distance.
                                            This value is used to decide if
                                            the specified VariantCall should be
                                            appended to an existing Variant.
                                            If there exists a VariantCall in a given Variant
                                            where the distances to both pos_1 and pos_2
                                            are equal to or less than max_neighbor_distance
                                            to pos_1 and pos_2 of the specified VariantCall,
                                            respectively, then the specified VariantCall is
                                            appended to the Variant. If such VariantCall
                                            is not identified, then a new Variant is constructed
                                            and added to self.variants.
        """
        # Step 1. Add variant_call if it can be appended to an existing Variant
        if self.size > 0:
            matched_indices = self.find_variant_indices(
                chromosome_1=variant_call.chromosome_1,
                chromosome_2=variant_call.chromosome_2,
                variant_types=VariantTypes.QueryTypeDictionary[variant_call.variant_type]
            )
            for idx in matched_indices:
                matched_variant_calls = self.variants[idx].find_variant_calls(
                    chromosome_1=variant_call.chromosome_1,
                    chromosome_2=variant_call.chromosome_2,
                    position_1_start=variant_call.position_1 - max_neighbor_distance,
                    position_1_end=variant_call.position_1 + max_neighbor_distance,
                    position_2_start=variant_call.position_2 - max_neighbor_distance,
                    position_2_end=variant_call.position_2 + max_neighbor_distance
                )
                if len(matched_variant_calls) > 0:
                    self.variants[idx].add_variant_call(variant_call=variant_call)
                    return

        # Add a new Variant
        variant = Variant(
            id='variant_%i' % (self.size + 1),
            chromosome_1=variant_call.chromosome_1,
            chromosome_2=variant_call.chromosome_2
        )
        variant.add_variant_call(variant_call=variant_call)
        insort(self.variants, variant)
        return

    def add_variant_call_to_variant(self, variant_call: VariantCall, variant_id: str):
        """
        Adds a VariantCall object to an existing Variant object.

        Parameters
        ----------
        variant_call                    :   VariantCall object.
        variant_id                      :   Variant ID.
        """
        for i in range(0, len(self.variants)):
            if self.variants[i].id == variant_id:
                self.variants[i].add_variant_call(variant_call=variant_call)
                return
        variant = Variant(
            id=variant_id,
            chromosome_1=variant_call.chromosome_1,
            chromosome_2=variant_call.chromosome_2
        )
        variant.add_variant_call(variant_call=variant_call)
        insort(self.variants, variant)
        return

    def find_variants(
            self,
            chromosome_1: str,
            chromosome_2: str,
            variant_types: List[str]
    ) -> List[Variant]:
        """
        Returns all variants that match chromosome_1 and chromosome_2.

        Parameters
        ----------
        chromosome_1    :   Chromosome 1.
        chromosome_2    :   Chromosome 2.
        variant_types   :   List of variant types to match.

        Returns
        -------
        variants        :   List of Variant objects.
        """
        # Find the leftmost index where chromosome_1 and chromosome_2 match
        left_index = bisect_left(
            self.variants,
            Variant(id='',
                    chromosome_1=chromosome_1,
                    chromosome_2=chromosome_2)
        )
        # Find the rightmost index where chromosome_1 and chromosome_2 match
        right_index = bisect_right(
            self.variants,
            Variant(
                id='',
                chromosome_1=chromosome_1,
                chromosome_2=chromosome_2
            )
        )
        variants = []
        for variant in self.variants[left_index:right_index]:
            if variant.chromosome_1 == chromosome_1 and \
                    variant.chromosome_2 == chromosome_2  and \
                    variant.variant_type in variant_types:
                variants.append(variant)
        return variants

    def find_variant_by_id(self, id: str) -> int:
        """
        Returns index of variant with query ID.

        Parameters
        ----------
        id          :   Variant ID.

        Returns
        -------
        index       :   Index to variant ID. Returns None if a matching
                        Variant object could not be found.
        """
        for i in range(0, self.variants):
            if self.variants[i].id == id:
                return i
        return None

    def find_variant_indices(
            self,
            chromosome_1: str,
            chromosome_2: str,
            variant_types: List[str]
    ) -> List[int]:
        """
        Returns indices of all variants that match chromosome_1 and chromosome_2.

        Parameters
        ----------
        chromosome_1    :   Chromosome 1.
        chromosome_2    :   Chromosome 2.
        variant_types   :   List of variant types.

        Returns
        -------
        indices         :   List of index integers.
        """
        # Find the leftmost index where chromosome_1 and chromosome_2 match
        left_index = bisect_left(
            self.variants,
            Variant(id='',
                    chromosome_1=chromosome_1,
                    chromosome_2=chromosome_2)
        )

        # Find the rightmost index where chromosome_1 and chromosome_2 match
        right_index = bisect_right(
            self.variants,
            Variant(id='',
                    chromosome_1=chromosome_1,
                    chromosome_2=chromosome_2)
        )

        # Get index values of Variant objects that match the query variant types
        indices = []
        for i in range(left_index, right_index):
            if self.variants[i].chromosome_1 == chromosome_1 and \
                    self.variants[i].chromosome_2 == chromosome_2  and \
                    self.variants[i].variant_type in variant_types:
                indices.append(i)
        return indices

    def filter_worker(
            self,
            variants: List[Variant],
            variant_filters: List[VariantFilter]
    ) -> List[str]:
        """
        Multiprocessing worker function for identifying variant IDs to remove.

        Parameters
        ----------
        variants                :   List of Variant objects.
        variant_filters         :   List of VariantFilter objects.

        Returns
        -------
        variant_ids_to_keep     :   List of variant IDs to keep.
        """
        variant_ids_to_keep = set()
        for variant in variants:
            remove = False
            for variant_filter in variant_filters:
                if not variant_filter.keep(variant=variant):
                    remove = True
            if not remove:
                variant_ids_to_keep.add(variant.id)
        return list(variant_ids_to_keep)

    def filter(
            self,
            variant_filters: List[VariantFilter],
            num_processes: int
    ) -> List[Variant]:
        """
        Filters variants based on VariantFilter objects.

        Parameters
        ----------
        variant_filters     :   List of VariantFilter objects.
        num_processes       :   Number of processes.

        Returns
        -------
        variants            :   List of Variant objects that satisfy all
                                supplied VariantFilter objects.
        """
        # Split the variants into multiple lists
        variants_list = np.array_split(self.variants, num_processes)

        # Multiprocess identification of which variant IDs to filter out
        pool = mp.Pool(processes=num_processes)
        async_results = [pool.apply_async(self.filter_worker, args=(variants, variant_filters)) for variants in variants_list]
        pool.close()
        pool.join()

        # Merge variant IDs to keep
        variant_ids_to_keep_list = [ar.get() for ar in async_results]
        variant_ids_to_keep = set()
        for curr_list in variant_ids_to_keep_list:
            for id in curr_list:
                variant_ids_to_keep.add(id)

        # Get variants that satisfy all supplied filters
        variants = []
        for variant in self.variants:
            if variant.id in variant_ids_to_keep:
                insort(variants, variant)

        return variants

    def filter_regions_worker(
            self,
            variants: List[Variant],
            genomic_ranges_list: GenomicRangesList,
            padding: int
    ) -> List[str]:
        """
        Multiprocessing worker function for identifying variant IDs
        that are near GenomicRangesList object.

        Parameters
        ----------
        variants                :   List of Variant objects.
        genomic_ranges_list     :   GenomicRangesList object.
        padding                 :   Padding.

        Returns
        -------
        variant_ids             :   List of variant IDs that are near any
                                    GenomicRange object in the queried
                                    GenomicRangesList object.
        """
        variant_ids = set()
        for variant in variants:
            for variant_call in variant.variant_calls:
                matched_genomic_ranges_1 = genomic_ranges_list.find_overlaps(
                    chromosome=variant_call.chromosome_1,
                    start=variant_call.position_1 - padding,
                    end=variant_call.position_1 + padding
                )
                matched_genomic_ranges_2 = genomic_ranges_list.find_overlaps(
                    chromosome=variant_call.chromosome_2,
                    start=variant_call.position_2 - padding,
                    end=variant_call.position_2 + padding
                )
                if len(matched_genomic_ranges_1) > 0 or len(matched_genomic_ranges_2) > 0:
                    variant_ids.add(variant.id)
                    break
        return list(variant_ids)

    def filter_regions(
            self,
            genomic_ranges_list: GenomicRangesList,
            padding: int,
            num_processes: int
    ) -> List[Variant]:
        """
        Filters variants based on GenomicRange objects.

        Parameters
        ----------
        genomic_ranges_list     :   GenomicRangesList object.
        padding                 :   Padding to apply to start and end
                                    positions of each GenomicRange object.
        num_processes           :   Number of processes.

        Returns
        -------
        variants                :   List of Variant objects that are near
                                    the queried GenomicRangeList object.
        """
        # Split the variants into multiple lists
        variants_list = np.array_split(self.variants, num_processes)

        # Multiprocess identification of which variant IDs are near
        # the queried GenomicRangeList object
        pool = mp.Pool(processes=num_processes)
        async_results = [pool.apply_async(self.filter_regions_worker, args=(variants, genomic_ranges_list, padding)) for variants in variants_list]
        pool.close()
        pool.join()

        # Merge variant IDs
        variant_ids_list = [ar.get() for ar in async_results]
        variant_ids = set()
        for curr_list in variant_ids_list:
            for id in curr_list:
                variant_ids.add(id)

        # Get variants that satisfy all supplied filters
        variants = []
        for variant in self.variants:
            if variant.id in variant_ids:
                insort(variants, variant)

        return variants

    def find_nearby_variants_worker(
            self,
            variants: List[Variant],
            padding: int
    ) -> List[str]:
        """
        Multiprocessing worker function for identifying variant IDs that
        are near query variants.

        Parameters
        ----------
        variants            :   List of Variant objects.
        padding             :   Padding.

        Returns
        -------
        variant_ids         :   List of variant IDs that are near query variants.
        """
        nearby_variant_ids = set()
        for variant in variants:
            variant_types = VariantTypes.QueryTypeDictionary[variant.variant_type]
            variants_ = self.find_variants(
                chromosome_1=variant.chromosome_1,
                chromosome_2=variant.chromosome_2,
                variant_types=variant_types
            )
            for variant_ in variants_:
                variant_calls = variant_.find_variant_calls(
                    chromosome_1=variant.chromosome_1,
                    chromosome_2=variant.chromosome_2,
                    position_1_start=min(variant.position_1) - padding,
                    position_1_end=max(variant.position_1) + padding,
                    position_2_start=min(variant.position_2) - padding,
                    position_2_end=max(variant.position_2) + padding
                )
                if len(variant_calls) > 0:
                    nearby_variant_ids.add(variant_.id)
        return list(nearby_variant_ids)

    def find_nearby_variants(
            self,
            variants: List[Variant],
            padding: int,
            num_processes: int,
    ) -> List[Variant]:
        """
        Finds variants that are near query Variant objects.

        Parameters
        ----------
        variants        :   List of Variant objects.
        padding         :   Padding to apply to position_1 and position_2 of each Variant.
        num_processes   :   Number of processes.

        Returns
        -------
        nearby_variants :   List of Variant objects.
        """
        # Split the variants into multiple lists
        variants_list = np.array_split(variants, num_processes)

        # Multiprocess identification of which variant IDs to filter out
        pool = mp.Pool(processes=num_processes)
        async_results = [pool.apply_async(self.find_nearby_variants_worker, args=(variants, padding)) for variants in variants_list]
        pool.close()
        pool.join()

        # Merge nearby variant IDs into one
        nearby_variant_ids_list = [ar.get() for ar in async_results]
        nearby_variant_ids = set()
        for nearby_variant_ids_ in nearby_variant_ids_list:
            for id in nearby_variant_ids_:
                nearby_variant_ids.add(id)

        # Get nearby variants
        nearby_variants = []
        for variant in self.variants:
            if variant.id in nearby_variant_ids:
                insort(nearby_variants, variant)

        return nearby_variants

    def remove(self, index):
        try:
            del self.variants[index]
        except:
            raise Exception('Invalid index: %i' % index)

    def remove_by_id(self, id):
        try:
            for i in range(0, self.size):
                if self.variants[i].id == id:
                    del self.variants[i]
                    return
        except:
            raise Exception('Invalid id: %s' % id)

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())

    def to_dict(self) -> Dict:
        data = defaultdict(list)
        for variant in self.variants:
            for key, values in variant.to_dict().items():
                for value in values:
                    data[key].append(value)
        return data



