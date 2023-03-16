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
The purpose of this python3 script is to implement common functions
related to handling VCF files.
"""


import gzip
import pandas as pd
from .logging import get_logger
from .constants import *
from .default_parameters import *
from .vcf_deepvariant import parse_deepvariant_callset
from .vcf_gatk4_mutect2 import parse_gatk4_mutect2_callset
from .vcf_strelka2 import parse_strelka2_somatic_callset, parse_strelka2_germline_callset
from .vcf_sniffles2 import parse_sniffles2_callset
from .vcf_pbsv import parse_pbsv_callset
from .vcf_cutesv import parse_cutesv_callset
from .vcf_svim import parse_svim_callset
from .variants_list import VariantsList


logger = get_logger(__name__)


def __read_vcf_file(vcf_file: str) -> pd.DataFrame:
    """
    Reads a VCF file and returns a DataFrame.

    Parameters
    ----------
    vcf_file    :   VCF file.

    Returns
    -------
    df_vcf      :   DataFrame.
    """
    vcf_names = []
    is_gzipped = False
    if vcf_file.endswith(".gz"):
        is_gzipped = True
        with gzip.open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_names = line.split('\t')
                    break
    else:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_names = line.split('\t')
                    break

    vcf_names = [i.replace('\n', '') for i in vcf_names]
    vcf_names = ['CHROM' if i == '#CHROM' else i for i in vcf_names]
    if is_gzipped:
        df_vcf = pd.read_csv(vcf_file,
                             compression='gzip',
                             comment='#',
                             delim_whitespace=True,
                             header=None,
                             names=vcf_names)
    else:
        df_vcf = pd.read_csv(vcf_file,
                             comment='#',
                             delim_whitespace=True,
                             header=None,
                             names=vcf_names)
    return df_vcf


def read_vcf_file(
        vcf_file: str,
        variant_calling_method: str,
        sequencing_platform: str,
        source_id: str,
        tumor_sample_id: str,
        normal_sample_id: str = ''
    ) -> VariantsList:
    """
    Reads a VCF file and returns an instance of the VariantsList class.

    Parameters
    ----------
    vcf_file                :   VCF file.
    variant_calling_method  :   Variant calling method.
    sequencing_platform     :   Sequencing platform.
    source_id               :   Source ID (e.g. patient ID or cell-line sample ID).
    tumor_sample_id         :   Tumor sample ID. If variant calling was performed
                                with a case (e.g. tumor) sample and
                                a control (e.g. normal) sample, then this parameter
                                must be specified to indicate which column corresponds
                                to the tumor sample.
    normal_sample_id        :   Normal sample ID. If variant calling was performed
                                with a case (e.g. tumor) sample and
                                a control (e.g. normal) sample, then this parameter
                                must be specified to indicate which column corresponds
                                to the normal sample.

    Returns
    -------
    variants_list           :   An instance of the VariantsList class.
    """
    # Step 1. Check input parameters
    if variant_calling_method not in VariantCallingMethods.ALL:
        raise Exception("The parameter 'variant_calling_method' must be one of the following: '%s'." %
                        ', '.join(VariantCallingMethods.ALL))

    # Step 2. Read the VCF file
    df_vcf = __read_vcf_file(vcf_file=vcf_file)

    # Step 3. Parse the DataFrame based on the variant calling method
    if variant_calling_method == VariantCallingMethods.DEEPVARIANT:
        variants_list = parse_deepvariant_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    if variant_calling_method == VariantCallingMethods.GATK4_MUTECT2:
        variants_list = parse_gatk4_mutect2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            tumor_sample_id=tumor_sample_id,
            normal_sample_id=normal_sample_id
        )
    if variant_calling_method == VariantCallingMethods.STRELKA2_SOMATIC:
        variants_list = parse_strelka2_somatic_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id,
            tumor_sample_id=tumor_sample_id,
            normal_sample_id=normal_sample_id
        )
    if variant_calling_method == VariantCallingMethods.STRELKA2_GERMLINE:
        variants_list = parse_strelka2_germline_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    if variant_calling_method == VariantCallingMethods.SNIFFLES2:
        variants_list = parse_sniffles2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    if variant_calling_method == VariantCallingMethods.PBSV:
        variants_list = parse_pbsv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    if variant_calling_method == VariantCallingMethods.CUTESV:
        variants_list = parse_cutesv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    if variant_calling_method == VariantCallingMethods.SVIM:
        variants_list = parse_svim_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    return variants_list

