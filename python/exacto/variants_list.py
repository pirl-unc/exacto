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


import json
import itertools
import pandas as pd
import multiprocessing as mp
import numpy as np
from collections import defaultdict, OrderedDict
from dataclasses import dataclass, field
from functools import partial
from typing import Dict, List, Tuple, Type
from .constants import *
from .genomic_range import GenomicRange
from .genomic_ranges_list import GenomicRangesList
from .logging import get_logger
from .utilities import get_typed_value, get_variant_calling_method_attr_types
from .variant import Variant
from .variant_annotation import VariantAnnotation
from .variant_call import VariantCall
from .variant_filter import VariantFilter
from exacto import exactors


logger = get_logger(__name__)


@dataclass(frozen=True)
class VariantsList:
    variants: List[Variant] = field(default_factory=list)
    _variants_dict: Dict[str, int] = field(default_factory=dict) # bookkeeping purposes

    def __post_init__(self):
        for i in range(0, len(self.variants)):
            self._variants_dict[self.variants[i].id] = i

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
        variant_ids = []
        for variant in self.variants:
            variant_ids.append(variant.id)
        return variant_ids

    def add_variant(self, variant: Variant):
        """
        Adds a Variant object.

        Parameters
        ----------
        variant     :   Variant object.
        """
        curr_len = len(self.variants)
        self._variants_dict[variant.id] = curr_len
        self.variants.append(variant)

    def get_variant(self, variant_id: str) -> Variant:
        return self.variants[self._variants_dict[variant_id]]

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

    def filter(
            self,
            variant_filters: List[VariantFilter],
            num_threads: int
    ) -> List[Variant]:
        """
        Filter variants by a list of VariantFilter objects
        and returns a list of Variant objects that pass
        all VariantFilter objects.

        Parameters
        ----------
        variant_filters         :   List of VariantFilter objects.
        num_threads             :   Number of threads.

        Returns
        -------
        filtered_variants_list  :   VariantsList object.
        """
        # Step 1. Serialize VariantsList object
        variants_list_serialized = json.dumps(self.to_dict())

        # Step 2. Serialize VariantFilter objects
        variant_filters_serialized = []
        for variant_filter in variant_filters:
            variant_filters_serialized.append(json.dumps(variant_filter.to_dict()))

        # Step 3. Filter VariantsList object
        json_str = exactors.filter_variants_list(
            variants_list_serialized,
            variant_filters_serialized,
            num_threads
        )

        # Step 4. Deserialize filtered VariantsList object
        filtered_variants_list = VariantsList.load_serialized_json(json_str=json_str)

        return filtered_variants_list.variants

    def find_nearby_variants(
            self,
            query_variants_list: Type["VariantsList"],
            num_threads: int,
            max_neighbor_distance: int) -> List[Tuple[Variant, List[Variant]]]:
        """
        Find nearby variants.

        Parameters
        ----------
        query_variants_list     :   Query VariantsList object.
        num_threads             :   Number of threads.
        max_neighbor_distance   :   Maximum neighbor distance.

        Returns
        -------
        nearby_variants         :   List where each tuple is
                                    (target Variant, list of nearby query Variants).
        """
        # Step 1. Serialize VariantsList object
        target_variants_list_serialized = json.dumps(self.to_dict())

        # Step 2. Serialize VariantsList object
        query_variants_list_serialized = json.dumps(query_variants_list.to_dict())

        # Step 3. Find nearby variants
        nearby_variants_dict = exactors.find_variants_near_query_variants(
            target_variants_list_serialized,
            query_variants_list_serialized,
            num_threads,
            max_neighbor_distance
        )

        # Step 4. Prepare returning data structure
        nearby_variants = []
        for key, value in nearby_variants_dict.items():
            target_variant = self.get_variant(variant_id=key)
            query_variants = []
            for query_variant_id in value:
                query_variants.append(query_variants_list.get_variant(query_variant_id))
            nearby_variants.append((target_variant, query_variants))
        return nearby_variants

    def overlap_regions(
            self,
            genomic_ranges_list: GenomicRangesList,
            padding: int,
            num_threads: int
    ) -> List[Tuple[Variant, List[GenomicRange]]]:
        """
        Filters variants based on GenomicRange objects.

        Parameters
        ----------
        genomic_ranges_list         :   GenomicRangesList object.
        padding                     :   Padding to apply to GenomicRange start and end.
        num_threads                 :   Number of threads.

        Returns
        -------
        variants                    :   List of Variant objects that are near
                                        the queried GenomicRangeList object.
        """
        # Step 1. Serialize VariantsList object
        variants_list_serialized = json.dumps(self.to_dict())

        # Step 2. Serialize GenomicRangesList object
        genomic_ranges_list_serialized = json.dumps(genomic_ranges_list.to_dict())

        # Step 3. Identify Variant objects that overlap GenomicRange objects
        nearby_variant_ids = exactors.find_variants_overlapping_regions(
            variants_list_serialized,
            genomic_ranges_list_serialized,
            num_threads,
            padding
        )

        # Step 4. Get Variant and GenomicRange objects
        nearby_variants = []
        for variant_id, genomic_range_ids in nearby_variant_ids.items():
            variant = self.get_variant(variant_id=variant_id)
            genomic_ranges = []
            for genomic_range_id in genomic_range_ids:
                genomic_range = genomic_ranges_list.get_genomic_range(genomic_range_id)
                genomic_ranges.append(genomic_range)
            nearby_variants.append((variant, genomic_ranges))
        return nearby_variants

    @staticmethod
    def load_dataframe(df: pd.DataFrame) -> Type["VariantsList"]:
        """
        Reads a DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df              :   DataFrame.

        Returns
        -------
        variants_list   :   An instance of the VariantsList class.
        """
        variants_dict = {}  # key = variant ID, value = Variant object
        for row in df.to_dict('records'):
            variant_id = VariantsList.convert_row_value(value=row['variant_id'], default_value='', type=str)
            variant_call_id = VariantsList.convert_row_value(value=row['variant_call_id'], default_value='', type=str)
            source_id = VariantsList.convert_row_value(value=row['source_id'], default_value='', type=str)
            sample_id = VariantsList.convert_row_value(value=row['sample_id'], default_value='', type=str)
            phase_block_id = VariantsList.convert_row_value(value=row['phase_block_id'], default_value='', type=str)
            clone_id = VariantsList.convert_row_value(value=row['clone_id'], default_value='', type=str)
            nucleic_acid = VariantsList.convert_row_value(value=row['nucleic_acid'], default_value='', type=str)
            variant_calling_method = VariantsList.convert_row_value(value=row['variant_calling_method'], default_value='', type=str)
            sequencing_platform = VariantsList.convert_row_value(value=row['sequencing_platform'], default_value='', type=str)
            chromosome_1 = VariantsList.convert_row_value(value=row['chromosome_1'], default_value='', type=str)
            position_1 = VariantsList.convert_row_value(value=row['position_1'], default_value=-1, type=int)
            chromosome_2 = VariantsList.convert_row_value(value=row['chromosome_2'], default_value='', type=str)
            position_2 = VariantsList.convert_row_value(value=row['position_2'], default_value=-1, type=int)
            reference_allele = VariantsList.convert_row_value(value=row['reference_allele'], default_value='', type=str)
            alternate_allele = VariantsList.convert_row_value(value=row['alternate_allele'], default_value='', type=str)
            filter = VariantsList.convert_row_value(value=row['filter'], default_value='', type=str)
            quality_score = VariantsList.convert_row_value(value=row['quality_score'], default_value=-1.0, type=float)
            precise = VariantsList.convert_row_value(value=row['precise'], default_value=False, type=bool)
            variant_type = VariantsList.convert_row_value(value=row['variant_type'], default_value='', type=str)
            variant_subtype = VariantsList.convert_row_value(value=row['variant_subtype'], default_value='', type=str)
            variant_size = VariantsList.convert_row_value(value=row['variant_size'], default_value=-1, type=int)
            total_read_count = VariantsList.convert_row_value(value=row['total_read_count'], default_value=-1, type=int)
            reference_allele_read_count = VariantsList.convert_row_value(value=row['reference_allele_read_count'], default_value=-1, type=int)
            alternate_allele_read_count = VariantsList.convert_row_value(value=row['alternate_allele_read_count'], default_value=-1, type=int)
            alternate_allele_fraction = VariantsList.convert_row_value(value=row['alternate_allele_fraction'], default_value=-1.0, type=float)
            variant_sequences_ = VariantsList.convert_row_value(value=row['variant_sequences'], default_value='', type=str)
            tags_ = VariantsList.convert_row_value(value=row['tags'], default_value='', type=str)
            alternate_allele_read_ids_ = VariantsList.convert_row_value(value=row['alternate_allele_read_ids'], default_value='', type=str)

            # Variant sequences
            variant_sequences = []
            if variant_sequences_ != '':
                for seq in variant_sequences_.split(';'):
                    variant_sequences.append(str(seq))

            # Alternate allele read IDs
            alternate_allele_read_ids = []
            if alternate_allele_read_ids_ != '':
                for read_id in alternate_allele_read_ids_.split(';'):
                    alternate_allele_read_ids.append(str(read_id))

            # Tool attributes
            tool_attributes = OrderedDict()
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
                        tool_attributes[curr_attr_key] = str(curr_attr_value)

            # Annotations
            # todo handle reading annotations
            # if row['pos_1_annotation_chrom'] != '':
            #     pos_1_annotation_counts = len(row['pos_1_annotation_chrom'].split(';'))
            #     for idx in range(0, pos_1_annotation_counts):
            #         pos_1_annotation_chrom = row['pos_1_annotation_chrom'].split(';')[idx]
            #         pos_1_annotation_chrom = row['pos_1_annotation_chrom'].split(';')[idx]
            #         pos_1_annotation_chrom = row['pos_1_annotation_chrom'].split(';')[idx]
            # annotation = VariantAnnotation()
            # self.__convert_tsv_file_element(value=row['pos_1_annotation_chrom'], nested=True, )

            # Tags
            tags = []
            if tags_ != '':
                tags = tags_.split(';')

            variant_call = VariantCall(
                id=variant_call_id,
                source_id=source_id,
                sample_id=sample_id,
                phase_block_id=phase_block_id,
                clone_id=clone_id,
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
                tool_attributes=tool_attributes,
                tags=tags
            )

            if variant_id not in variants_dict.keys():
                variants_dict[variant_id] = Variant(id=variant_id)
            variants_dict[variant_id].add_variant_call(variant_call=variant_call)

        variants_list = VariantsList()
        for variant in variants_dict.values():
            variants_list.add_variant(variant=variant)

        logger.info("Loaded %i variants and %i variant calls." % (variants_list.size, len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def load_serialized_json(json_str: str) -> Type["VariantsList"]:
        """
        Loads a VaraintsList object from a serialized JSON string.

        Parameters
        ----------
        json_str        :   JSON string.

        Returns
        -------
        variants_list   :   VariantsList object.
        """
        variants_list = VariantsList()
        variants_list_dict = json.loads(json_str)
        for variant_dict in variants_list_dict['variants']:
            variant = Variant(id=variant_dict['id'])
            for variant_call_dict in variant_dict['variant_calls']:
                position_1_annotations_dict = variant_call_dict['position_1_annotations']
                position_2_annotations_dict = variant_call_dict['position_2_annotations']
                del variant_call_dict['position_1_annotations']
                del variant_call_dict['position_2_annotations']
                variant_call = VariantCall(**variant_call_dict)
                for position_1_annotation_dict in position_1_annotations_dict:
                    variant_call.add_position_1_annotation(
                        variant_annotation=VariantAnnotation(**position_1_annotation_dict)
                    )
                for position_2_annotation_dict in position_2_annotations_dict:
                    variant_call.add_position_2_annotation(
                        variant_annotation=VariantAnnotation(**position_2_annotation_dict)
                    )
                variant.add_variant_call(variant_call=variant_call)
            variants_list.add_variant(variant=variant)
        return variants_list

    @staticmethod
    def merge(variants_lists: List[Type["VariantsList"]],
              num_threads: int,
              max_neighbor_distance: int) -> Type["VariantsList"]:
        """
        Merges a list of VariantsList objects and returns one VariantsList object.

        Parameters
        ----------
        variants_lists                  :   List of VariantsList objects.
        num_threads                     :   Number of threads.
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
        # Step 1. Serialize all VariantsList objects
        variants_lists_serialized = []
        for variants_list in variants_lists:
            variants_lists_serialized.append(json.dumps(variants_list.to_dict()))

        # Step 2. Merge VariantsList objects
        json_str = exactors.merge_variants_lists(
            variants_lists_serialized,
            num_threads,
            max_neighbor_distance
        )
        merged_variants_list = VariantsList.load_serialized_json(json_str=json_str)
        return merged_variants_list

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
        df = pd.read_csv(tsv_file, sep='\t', low_memory=False, memory_map=True)
        return VariantsList.load_dataframe(df=df)

    # def remove(self, variant: Variant):
    #     variant_query_types = ','.join(VariantTypes.QueryTypeDictionary[variant.variant_type])
    #     key = '%s-%s-%s' % (variant.chromosome_1, variant.chromosome_2, variant_query_types)
    #     try:
    #         index = self.variants[key].index(variant)
    #         del self.variants[key][index]
    #     except:
    #         raise Exception('Variant object does not exist: %s' % variant)

    def to_dict(self):
        data = {
            'variants': [variant.to_dict() for variant in self.variants]
        }
        return data

    def to_dataframe_row(self) -> Dict:
        data = defaultdict(list)
        for variant in self.variants:
            for key, values in variant.to_dataframe_row().items():
                for value in values:
                    data[key].append(value)
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())
