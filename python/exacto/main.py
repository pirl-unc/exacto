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
The purpose of this python3 script is to implement Exacto's main APIs.
"""


import itertools
import pandas as pd
import pysam
from collections import defaultdict
from dataclasses import field
from typing import List, Tuple
from .annotation_db import AnnotationDb
from .constants import NucleicAcidTypes, VariantCallingMethods, VariantCallTags
from .default import *
from .genomic_ranges_list import GenomicRangesList
from .logging import get_logger
from .variants_list import VariantsList
from .variant_filter import VariantFilter
from .vcf import Vcf


logger = get_logger(__name__)


def run_exacto_vcf2tsv(
        vcf_file: str,
        source_id: str,
        variant_calling_method: str,
        sequencing_platform: str,
        case_id: str = '',
        control_id: str = ''
) -> VariantsList:
    """
    Convert a VCF file to a VariantsList.

    Parameters
    ----------
    vcf_file                :   VCF file.
    source_id               :   Source ID (e.g. patient ID or cell line sample ID).
    variant_calling_method  :   Variant calling method.
    sequencing_platform     :   Sequencing platform.
    case_id                 :   Case ID (only necessary if variant_calling_method is 'strelka2-somatic').
    control_id              :   Control ID (only necessary if variant_calling_method is 'strelka2-somatic').

    Returns
    -------
    variants_list           :   VariantsList object.
    """
    df_vcf = Vcf.read_vcf_file(vcf_file=vcf_file)
    if variant_calling_method == VariantCallingMethods.CUTESV:
        variants_list = Vcf.parse_cutesv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DEEPVARIANT:
        variants_list = Vcf.parse_deepvariant_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DELLY2_SOMATIC:
        variants_list = Vcf.parse_delly2_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.GATK4_MUTECT2:
        variants_list = Vcf.parse_gatk4_mutect2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.LUMPY_SOMATIC:
        variants_list = Vcf.parse_lumpy_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.PBSV:
        variants_list = Vcf.parse_pbsv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.SNIFFLES2:
        variants_list = Vcf.parse_sniffles2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.STRELKA2_SOMATIC:
        variants_list = Vcf.parse_strelka2_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            case_id=case_id,
            control_id=control_id
        )
    elif variant_calling_method == VariantCallingMethods.SVIM:
        variants_list = Vcf.parse_svim_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DBSNP:
        variants_list = Vcf.parse_dbsnp_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    else:
        raise Exception('Unsupported variant calling method: %s' % variant_calling_method)
    return variants_list


def run_exacto_filter_variants(
        variants_list: VariantsList,
        variant_filters: List[VariantFilter] = None,
        excluded_variants_list: VariantsList = None,
        excluded_regions_list: GenomicRangesList = None,
        excluded_variants_padding: int = FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING,
        excluded_regions_padding: int = FILTER_VARIANTS_EXCLUDED_REGION_PADDING,
        num_threads: int = FILTER_VARIANTS_NUM_THREADS
) -> Tuple[VariantsList, VariantsList]:
    """
    Filters a VariantsList object.

    Parameters
    ----------
    variants_list               :   VariantsList object.
    variant_filters             :   List of VariantFilter objects.
    excluded_variants_list      :   VariantsList object of variants to exclude.
    excluded_regions_list       :   GenomicRangesList object of regions to exclude.
    excluded_variants_padding   :   Number of bases to pad each variant's positions 1 and 2.
    excluded_regions_padding    :   Number of bases to pad each region to exclude.
    num_threads                 :   Number of threads.

    Returns
    -------
    variants_list_filtered      :   VariantsList object (of filtered variants).
    variants_list_rejected      :   VariantsList object (of rejected variants).
    """
    logger.info('%i variants in the original list before any filtering.' % variants_list.size)
    logger.info('%i variant calls in the original list before any filtering.' % len(variants_list.variant_call_ids))

    # Step 1. Filter out variants based on VariantFilter
    # key   = variant ID
    # value = reasons why the variant was rejected
    rejected_variant_ids_dict = defaultdict(list)
    if variant_filters is not None:
        filtered_variants = variants_list.filter(
            variant_filters=variant_filters,
            num_threads=num_threads,
        )
        filtered_variants_ids = set([variant.id for variant in filtered_variants])
        for variant_id in variants_list.variant_ids:
            if variant_id not in filtered_variants_ids:
                rejected_variant_ids_dict[variant_id].append(VariantCallTags.FAILED_FILTER)
        logger.info('%i variants satisfy all variant filters.' % len(filtered_variants))

    # Step 2. Filter out variants overlapping the excluded regions
    if excluded_regions_list is not None:
        filtered_variants = variants_list.overlap_regions(
            genomic_ranges_list=excluded_regions_list,
            padding=excluded_regions_padding,
            num_threads=num_threads
        )
        for filtered_variant, genomic_ranges in filtered_variants:
            rejected_variant_ids_dict[filtered_variant.id].append(VariantCallTags.NEARBY_EXCLUDED_REGION)
        logger.info('%i variants are near excluded regions.' % len(filtered_variants))

    # Step 3. Filter out variants near the excluded variants
    if excluded_variants_list is not None:
        filtered_variants = variants_list.find_nearby_variants(
            query_variants_list=excluded_variants_list,
            max_neighbor_distance=excluded_variants_padding,
            num_threads=num_threads
        )
        for filtered_variant, query_variants in filtered_variants:
            rejected_variant_ids_dict[filtered_variant.id].append(VariantCallTags.NEARBY_EXCLUDED_VARIANT)
        logger.info('%i variants are near excluded variants.' % len(filtered_variants))

    # Step 4. Create a filtered VariantsList and a rejected VariantsList
    variants_list_filtered = VariantsList()
    variants_list_rejected = VariantsList()
    for variant in variants_list.variants:
        if variant.id in rejected_variant_ids_dict.keys():
            for i in range(0, variant.num_variant_calls):
                for reason in rejected_variant_ids_dict[variant.id]:
                    variant.variant_calls[i].tags.append(reason)
            variants_list_rejected.add_variant(variant=variant)
        else:
            for i in range(0, variant.num_variant_calls):
                variant.variant_calls[i].tags.append(VariantCallTags.PASSED)
            variants_list_filtered.add_variant(variant=variant)

    logger.info('%i variants in the filtered list after all filtering.' % variants_list_filtered.size)
    logger.info('%i variant calls in the filtered list after all filtering.' % len(variants_list_filtered.variant_call_ids))

    return variants_list_filtered, variants_list_rejected


def run_exacto_merge_variants(
        variants_lists: List[VariantsList],
        num_threads: int = MERGE_VARIANTS_NUM_THREADS,
        max_neighbor_distance: int = MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE,
) -> VariantsList:
    """
    Merges VariantsList objects into one.

    Parameters
    ----------
    variants_lists                  :   List of VariantsList objects.
    num_threads                     :   Number of threads.
    max_neighbor_distance           :   Maximum neighbor distance.

    Returns
    -------
    variants_list                   :   VariantsList object.
    """
    variants_list = VariantsList.merge(
        variants_lists=variants_lists,
        num_threads=num_threads,
        max_neighbor_distance=max_neighbor_distance
    )
    return variants_list


def run_exacto_annotate_variants_list(
        variants_list: VariantsList,
        annotation_db: AnnotationDb
) -> VariantsList:
    """
    Annotates a variants list and returns the annotated variants list.

    Parameters
    ----------
    variants_list       :   VariantsList object.
    annotation_db       :   AnnotationDb object.

    Returns
    -------
    variants_list       :   VariantsList object.
    """
    return annotation_db.annotate_variants_list(variants_list=variants_list)
