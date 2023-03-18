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
from copy import deepcopy
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Type, Dict
from .common import safely_convert_value, get_variant_calling_method_attr_types, overlaps_any
from .constants import *
from .default_parameters import *
from .variant import Variant
from .variant_call import VariantCall
from .variant_filter import VariantFilter
from .variant_annotation import VariantAnnotation
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class VariantsList:
    variants: List[Variant] = field(default_factory=list)

    @property
    def variant_ids(self) -> List[str]:
        return [variant.id for variant in self.variants]

    @property
    def variant_call_ids(self) -> List[str]:
        variant_call_ids = []
        for variant in self.variants:
            for variant_call in variant.variant_calls:
                variant_call_ids.append(variant_call.id)
        return variant_call_ids

    @property
    def size(self):
        return len(self.variant_ids)

    @staticmethod
    def merge(
            variants_lists: List,
            max_neighbor_distance: int,
            enforce_variant_type_matching: bool = True
        ) -> Type["VariantsList"]:
        """
        Merges a list of VariantsList instances into one.

        Parameters
        ----------
        variants_lists                  :   List of VariantsList instances.
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
        enforce_variant_type_matching   :   If true, variant_type must match for 2 variant calls
                                            to be considered to be in the same variant (default: True).

        Returns
        -------
        variants_list   :   An instance of VariantsList.
        """
        if len(variants_lists) > 1:
            variants_list = deepcopy(variants_lists[0])
            for i in range(1, len(variants_lists)):
                for j in range(0, len(variants_lists[i].variants)):
                    for k in range(0, len(variants_lists[i].variants[j].variant_calls)):
                        variants_list.add_variant_call(
                            variant_call=variants_lists[i].variants[j].variant_calls[k],
                            enforce_variant_type_matching=enforce_variant_type_matching,
                            max_neighbor_distance=max_neighbor_distance
                        )
            return variants_list
        else:
            return variants_lists[0]

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
            if type == str:
                value = str(value)
            if type == int:
                value = int(value)
            if type == float:
                value = float(value)
            if type == bool:
                value = bool(value)
            return value

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
        variants_list = VariantsList()
        for index, row in df.iterrows():
            variant_call = VariantCall()
            variant_call.id = VariantsList.convert_row_value(value=row['variant_call_id'], default_value=None, type=str)
            variant_call.source_id = VariantsList.convert_row_value(value=row['source_id'], default_value=None, type=str)
            variant_call.tumor_sample_id = VariantsList.convert_row_value(value=row['tumor_sample_id'], default_value=None, type=str)
            variant_call.normal_sample_id = VariantsList.convert_row_value(value=row['normal_sample_id'], default_value=None, type=str)
            variant_call.nucleic_acid = VariantsList.convert_row_value(value=row['nucleic_acid'], default_value=None, type=str)
            variant_call.variant_calling_method = VariantsList.convert_row_value(value=row['variant_calling_method'], default_value=None, type=str)
            variant_call.sequencing_platform = VariantsList.convert_row_value(value=row['sequencing_platform'], default_value=None, type=str)
            variant_call.chr_1 = VariantsList.convert_row_value(value=row['chr_1'], default_value=None, type=str)
            variant_call.pos_1 = VariantsList.convert_row_value(value=row['pos_1'], default_value=None, type=int)
            variant_call.chr_2 = VariantsList.convert_row_value(value=row['chr_2'], default_value=None, type=str)
            variant_call.pos_2 = VariantsList.convert_row_value(value=row['pos_2'], default_value=None, type=int)
            variant_call.ref = VariantsList.convert_row_value(value=row['ref'], default_value=None, type=str)
            variant_call.alt = VariantsList.convert_row_value(value=row['alt'], default_value=None, type=str)
            variant_call.filter = VariantsList.convert_row_value(value=row['filter'], default_value=None, type=str)
            variant_call.quality_score = VariantsList.convert_row_value(value=row['quality_score'], default_value=None, type=float)
            variant_call.precise = VariantsList.convert_row_value(value=row['precise'], default_value=None, type=bool)
            variant_call.variant_type = VariantsList.convert_row_value(value=row['variant_type'], default_value=None, type=str)
            variant_call.variant_subtype = VariantsList.convert_row_value(value=row['variant_subtype'], default_value=None, type=str)
            variant_call.variant_size = VariantsList.convert_row_value(value=row['variant_size'], default_value=None, type=int)
            variant_call.total_tumor_reads = VariantsList.convert_row_value(value=row['total_tumor_reads'], default_value=None, type=int)
            variant_call.ref_tumor_reads = VariantsList.convert_row_value(value=row['ref_tumor_reads'], default_value=None, type=int)
            variant_call.alt_tumor_reads = VariantsList.convert_row_value(value=row['alt_tumor_reads'], default_value=None, type=int)
            variant_call.other_tumor_reads = VariantsList.convert_row_value(value=row['other_tumor_reads'], default_value=None, type=int)
            variant_call.alt_tumor_fraction = VariantsList.convert_row_value(value=row['alt_tumor_fraction'], default_value=None, type=float)
            variant_call.total_normal_reads = VariantsList.convert_row_value(value=row['total_normal_reads'], default_value=None, type=int)
            variant_call.ref_normal_reads = VariantsList.convert_row_value(value=row['ref_normal_reads'], default_value=None, type=int)
            variant_call.alt_normal_reads = VariantsList.convert_row_value(value=row['alt_normal_reads'], default_value=None, type=int)
            variant_call.other_normal_reads = VariantsList.convert_row_value(value=row['other_normal_reads'], default_value=None, type=int)
            variant_call.alt_normal_fraction = VariantsList.convert_row_value(value=row['alt_normal_fraction'], default_value=None, type=float)
            variant_call.alt_tumor_softclip_direction = VariantsList.convert_row_value(value=row['alt_tumor_softclip_direction'], default_value=None, type=str)
            variant_call.alt_normal_softclip_direction = VariantsList.convert_row_value(value=row['alt_normal_softclip_direction'], default_value=None, type=str)

            variant_sequences = VariantsList.convert_row_value(value=row['variant_sequences'], default_value='', type=str)
            alt_tumor_read_ids = VariantsList.convert_row_value(value=row['alt_tumor_read_ids'], default_value='', type=str)
            alt_normal_read_ids = VariantsList.convert_row_value(value=row['alt_normal_read_ids'], default_value='', type=str)

            if variant_sequences != '':
                variant_call.variant_sequences = variant_sequences.split(';')
            if alt_tumor_read_ids != '':
                variant_call.alt_tumor_read_ids = alt_tumor_read_ids.split(';')
            if alt_normal_read_ids != '':
                variant_call.alt_normal_read_ids = alt_normal_read_ids.split(';')

            # Tool attributes
            if row['tool_attributes'] != '':
                curr_tool_attr_types = get_variant_calling_method_attr_types(variant_calling_method=variant_call.variant_calling_method)
                for curr_attr in row['tool_attributes'].split(';'):
                    curr_attr_key = curr_attr.split('=')[0]
                    curr_attr_key_query = curr_attr_key.replace('_tumor_', '_')
                    curr_attr_key_query = curr_attr_key_query.replace('_normal_', '_')
                    curr_attr_value = curr_attr.split('=')[1]
                    curr_attr_value = safely_convert_value(
                        value=curr_attr_value,
                        default_value=None,
                        type=curr_tool_attr_types[curr_attr_key_query]
                    )
                    if curr_attr_value is not None:
                        variant_call.tool_attributes[curr_attr_key] = curr_attr_value

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

            if row['variant_id'] not in variants_list.variant_ids:
                variant = Variant()
                variant.id = row['variant_id']
                variant.variant_calls.append(variant_call)
                variants_list.variants.append(variant)
            else:
                for i in range(0, len(variants_list.variants)):
                    if variants_list.variants[i].id == row['variant_id']:
                        variants_list.variants[i].variant_calls.append(variant_call)
        return variants_list

    def add_variant_call(
            self,
            variant_call: VariantCall,
            max_neighbor_distance: int,
            enforce_variant_type_matching: bool = True
        ):
        """
        Adds a VariantCall object.

        Parameters
        ----------
        variant_call                    :   An instance of the VariantCall class.
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
        enforce_variant_type_matching   :   If true, variant_type must match for 2 variant calls
                                            to be considered to be in the same variant (default: True).
        """
        # Add variant_call if it can be appended to an existing Variant
        query_variant_type = VariantTypes.QueryTypeDictionary[variant_call.variant_type]
        matched_variant_ids = []
        for i in range(0, len(self.variants)):
            for j in range(0, len(self.variants[i].variant_calls)):
                pos_1_delta = abs(self.variants[i].variant_calls[j].pos_1 - variant_call.pos_1)
                pos_2_delta = abs(self.variants[i].variant_calls[j].pos_2 - variant_call.pos_2)
                if (self.variants[i].variant_calls[j].chr_1 == variant_call.chr_1) and \
                        (self.variants[i].variant_calls[j].chr_2 == variant_call.chr_2) and \
                        (pos_1_delta <= max_neighbor_distance) and \
                        (pos_2_delta <= max_neighbor_distance):
                    if enforce_variant_type_matching:
                        if self.variants[i].variant_calls[j].variant_type in query_variant_type:
                            self.variants[i].variant_calls.append(variant_call)
                            matched_variant_ids.append(self.variants[i].variant_calls[j].id)
                            return
                    else:
                        self.variants[i].variant_calls.append(variant_call)
                        matched_variant_ids.append(self.variants[i].variant_calls[j].id)
                        return

        # Add a new Variant
        variant = Variant(id='variant_%i' % (len(self.variants) + 1))
        variant.variant_calls.append(variant_call)

    def load_annotations(self):
        # todo load Ensembl or Gencode annotations into each VariantAnnotation Gene
        pass

    def filter(self, variant_filters: List[VariantFilter]):
        """
        Applies a filter condition to self.variants

        Parameters
        ----------
        variant_call_filters :   A list of instances of the VariantCallFilter class.
        """
        # Step 1. Identify variants to remove
        variant_ids_to_remove = set()
        for variant_filter in variant_filters:
            for variant in self.variants:
                if not variant_filter.is_predicate(variant=variant):
                    variant_ids_to_remove.add(variant.id)

        # Step 2. Remove variants
        variants = []
        for variant in self.variants:
            if variant.id not in variant_ids_to_remove:
                variants.append(variant)
        self.variants = variants

    def filter_regions(
            self,
            df_excluded_regions: pd.DataFrame,
            excluded_regions_padding: int
        ):
        """
        Filters variant calls that are near excluded regions.

        Parameters
        ----------
        df_excluded_regions         :   DataFrame of regions to exclude.
                                        Expected headers: 'chrom', 'chromStart', 'chromEnd'
        excluded_regions_padding    :   Number of bases to pad each excluded region.
        """
        # Step 1. Apply padding to excluded regions
        df_excluded_regions['pos_1'] = df_excluded_regions['chromStart'] - excluded_regions_padding
        df_excluded_regions['pos_2'] = df_excluded_regions['chromEnd'] + excluded_regions_padding
        df_excluded_regions['chr_1'] = df_excluded_regions['chrom']
        df_excluded_regions['chr_2'] = df_excluded_regions['chrom']

        # Step 2. Filter variant calls
        variants = []
        for variant in self.variants:
            # Identify variant call IDs to remove
            variant_call_ids_to_remove = set()
            for variant_call in variant.variant_calls:
                pos_1_overlaps = overlaps_any(
                    df=df_excluded_regions,
                    chrom=variant_call.chr_1,
                    start=variant_call.pos_1,
                    end=variant_call.pos_1
                )
                pos_2_overlaps = overlaps_any(
                    df=df_excluded_regions,
                    chrom=variant_call.chr_2,
                    start=variant_call.pos_2,
                    end=variant_call.pos_2
                )
                if pos_1_overlaps or pos_2_overlaps:
                    variant_call_ids_to_remove.add(variant_call.id)

            variant_ = Variant()
            variant_.id = variant.id
            for variant_call in variant.variant_calls:
                if variant_call.id not in variant_call_ids_to_remove:
                    variant_.variant_calls.append(variant_call)
            if len(variant_.variant_calls) > 0:
                variants.append(variant_)

        self.variants = variants

    def filter_variants(
            self,
            df_excluded_variants: pd.DataFrame,
            excluded_variant_padding: int,
            enforce_variant_type_checking: bool
        ):
        """
        Filters variant calls that are near excluded variants.

        Parameters
        ----------
        df_excluded_variants            :   DataFrame. Expected columns are:
                                            'chr_1', 'pos_1', 'chr_2', 'pos_2'
        excluded_variant_padding        :   Number of bases to pad each excluded variant.
        enforce_variant_type_checking   :   Enforce variant type checking.
        """
        # Step 1. Apply padding to excluded variants
        df_excluded_variants['pos_1_start'] = df_excluded_variants.apply(
            lambda row: int(row.pos_1 - excluded_variant_padding), axis=1
        )
        df_excluded_variants['pos_1_end'] = df_excluded_variants.apply(
            lambda row: int(row.pos_1 + excluded_variant_padding), axis=1
        )
        df_excluded_variants['pos_2_start'] = df_excluded_variants.apply(
            lambda row: int(row.pos_2 - excluded_variant_padding), axis=1
        )
        df_excluded_variants['pos_2_end'] = df_excluded_variants.apply(
            lambda row: int(row.pos_2 + excluded_variant_padding), axis=1
        )

        # Step 2. Filter variants
        variants = []
        for variant in self.variants:
            # Identify variant call IDs to remove
            variant_call_ids_to_remove = set()
            for variant_call in variant.variant_calls:
                if enforce_variant_type_checking:
                    query_variant_type = VariantTypes.QueryTypeDictionary[variant_call.variant_type]
                    conditions = \
                        ((df_excluded_variants['chr_1'] == variant_call.chr_1) &
                         (df_excluded_variants['pos_1_start'] <= variant_call.pos_1) &
                         (df_excluded_variants['pos_1_end'] >= variant_call.pos_1) &
                         (df_excluded_variants['variant_type'].isin(query_variant_type))) | \
                        ((df_excluded_variants['chr_2'] == variant_call.chr_1) &
                         (df_excluded_variants['pos_2_start'] <= variant_call.pos_1) &
                         (df_excluded_variants['pos_2_end'] >= variant_call.pos_1) &
                         (df_excluded_variants['variant_type'].isin(query_variant_type))) | \
                        ((df_excluded_variants['chr_1'] == variant_call.chr_2) &
                         (df_excluded_variants['pos_1_start'] <= variant_call.pos_2) &
                         (df_excluded_variants['pos_1_end'] >= variant_call.pos_2) &
                         (df_excluded_variants['variant_type'].isin(query_variant_type))) | \
                        ((df_excluded_variants['chr_2'] == variant_call.chr_2) &
                         (df_excluded_variants['pos_2_start'] <= variant_call.pos_2) &
                         (df_excluded_variants['pos_2_end'] >= variant_call.pos_2) &
                         (df_excluded_variants['variant_type'].isin(query_variant_type)))
                    df_matched = df_excluded_variants[conditions]
                else:
                    conditions = \
                        ((df_excluded_variants['chr_1'] == variant_call.chr_1) &
                         (df_excluded_variants['pos_1_start'] <= variant_call.pos_1) &
                         (df_excluded_variants['pos_1_end'] >= variant_call.pos_1)) | \
                        ((df_excluded_variants['chr_2'] == variant_call.chr_1) &
                         (df_excluded_variants['pos_2_start'] <= variant_call.pos_1) &
                         (df_excluded_variants['pos_2_end'] >= variant_call.pos_1)) | \
                        ((df_excluded_variants['chr_1'] == variant_call.chr_2) &
                         (df_excluded_variants['pos_1_start'] <= variant_call.pos_2) &
                         (df_excluded_variants['pos_1_end'] >= variant_call.pos_2)) | \
                        ((df_excluded_variants['chr_2'] == variant_call.chr_2) &
                         (df_excluded_variants['pos_2_start'] <= variant_call.pos_2) &
                         (df_excluded_variants['pos_2_end'] >= variant_call.pos_2))
                    df_matched = df_excluded_variants[conditions]
                if len(df_matched) > 0:
                    variant_call_ids_to_remove.add(variant_call.id)

            variant_ = Variant()
            variant_.id = variant.id
            for variant_call in variant.variant_calls:
                if variant_call.id not in variant_call_ids_to_remove:
                    variant_.variant_calls.append(variant_call)
            if len(variant_.variant_calls) > 0:
                variants.append(variant_)
        self.variants = variants

    def to_dict(self) -> Dict:
        data = defaultdict(list)
        for variant in self.variants:
            for key, value in variant.to_dict().items():
                data[key].append(value[0])
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())

