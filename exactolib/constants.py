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
The purpose of this python3 script is to define Exacto constants.
"""


class SequencingPlatforms:
    ILLUMINA = 'illumina'
    PACBIO_HIFI_CCS = 'pacbio-hifi-ccs'
    OXFORD_NANOPORE_TECHNOLOGIES = 'ont'
    ALL = [ILLUMINA,
           PACBIO_HIFI_CCS,
           OXFORD_NANOPORE_TECHNOLOGIES]


class NucleicAcidTypes:
    DNA = 'dna'
    RNA = 'rna'
    ALL = [DNA,
           RNA]


class VariantTypes:
    SV = 'sv'
    SNV_INDEL = 'snv_indel'
    SNV = 'snv'
    INSERTION = 'ins'
    DELETION = 'del'
    INVERSION_5TO5 = 'inv_5to5'
    INVERSION_5TO3 = 'inv_5to3'
    INVERSION_3TO3 = 'inv_3to3'
    TANDEM_DUPLICATION = 'tandem_dup'
    SEGMENTAL_DUPLICATION = 'seg_dup'
    TRANSLOCATION = 'tra'
    SPLICE_VARIANT = 'spvar'
    INTRON_RETENTION = 'intron_ret'
    CIRCULAR_RNA = 'circ_rna'
    ENDOGENOUS_RETROVIRUS = 'erv'
    ALL = [SV,
           SNV_INDEL,
           SNV,
           INSERTION,
           DELETION,
           INVERSION_5TO5,
           INVERSION_5TO3,
           INVERSION_3TO3,
           TANDEM_DUPLICATION,
           SEGMENTAL_DUPLICATION,
           TRANSLOCATION,
           SPLICE_VARIANT,
           INTRON_RETENTION,
           CIRCULAR_RNA,
           ENDOGENOUS_RETROVIRUS]


class VariantCallingMethods:

    class StructuralVariantCallingMethods:
        SNIFFLES = 'sniffles'
        SNIFFLES2 = 'sniffles2'
        SVIM = 'svim'
        CUTESV = 'cutesv'
        DELLY2 = 'delly2'
        LUMPY = 'lumpy'
        PBSV = 'pbsv'
        ALL = [SNIFFLES,
               SNIFFLES2,
               SVIM,
               CUTESV,
               DELLY2,
               LUMPY,
               PBSV]

    class SmallVariantCallingMethods:
        GATK4_MUTECT2 = "gatk4_mutect2"
        DEEPVARIANT = "deepvariant"
        STRELKA2 = "strelka2"
        ALL = [GATK4_MUTECT2,
               DEEPVARIANT,
               STRELKA2]

    ALL = StructuralVariantCallingMethods.ALL + \
          SmallVariantCallingMethods.ALL


class AnnotationSources:
    ENSEMBL = 'ensembl'
    GENCODE = 'gencode'
    ANNOVAR = 'annovar'
    ALL = [ENSEMBL, GENCODE, ANNOVAR]


class StructuralVariantTypes:
    INSERTION = 'INS'
    DELETION = 'DEL'
    DUPLICATION = 'DUP'
    INVERSION = 'INV'
    TRANSLOCATION = 'TRA'
    BREAKPOINT =' BND'
    ALL = [
        INSERTION,
        DELETION,
        DUPLICATION,
        INVERSION,
        TRANSLOCATION,
        BREAKPOINT
    ]

class SmallVariantTypes:
    SINGLE_NUCLEOTIDE_VARIANT = 'SNV'
    MULTI_NUCLEOTIDE_VARIANT = 'MNV'
    SMALL_INSERTION = 'INS'
    SMALL_DELETION = 'DEL'
    ALL = [
        SINGLE_NUCLEOTIDE_VARIANT,
        MULTI_NUCLEOTIDE_VARIANT,
        SMALL_INSERTION,
        SMALL_DELETION
    ]
