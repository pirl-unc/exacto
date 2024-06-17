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


class VariantTypes:
    SINGLE_NUCLEOTIDE_VARIANT = 'SNV'   # DNA
    MULTI_NUCLEOTIDE_VARIANT = 'MNV'    # DNA
    INSERTION = 'INS'                   # DNA
    DELETION = 'DEL'                    # DNA
    INVERSION = 'INV'                   # DNA
    DUPLICATION = 'DUP'                 # DNA
    TRANSLOCATION = 'TRA'               # DNA
    BREAKPOINT = 'BND'                  # DNA
    SPLICE_VARIANT = 'SPV'              # RNA
    NOVEL_ISOFORM = 'NIS'               # RNA
    INTRON_RETENTION = 'INR'            # RNA
    FUSION_GENE = 'FUS'                 # RNA
    CIRCULAR_RNA = 'CIR'                # RNA
    JET = 'JET'                         # RNA
    ENDOGENOUS_RETROVIRUS = 'ERV'       # RNA
    REFERENCE = 'REF'                   # REFERENCE
    ALL = [
        SINGLE_NUCLEOTIDE_VARIANT,
        MULTI_NUCLEOTIDE_VARIANT,
        INSERTION,
        DELETION,
        INVERSION,
        DUPLICATION,
        TRANSLOCATION,
        BREAKPOINT,
        SPLICE_VARIANT,
        NOVEL_ISOFORM,
        INTRON_RETENTION,
        FUSION_GENE,
        CIRCULAR_RNA,
        JET,
        ENDOGENOUS_RETROVIRUS
    ]

    FULL_NAMES = {
        REFERENCE: 'Reference',
        SINGLE_NUCLEOTIDE_VARIANT: 'Single-nucleotide Variant',
        INSERTION: 'Insertion',
        DELETION: 'Deletion',
        DUPLICATION: 'Duplication',
        INVERSION: 'Inversion',
        TRANSLOCATION: 'Translocation'
    }

    class DeletionSubtypes:
        EXONIC_DELETION = 'DEL_EXONIC'
        PARTIAL_EXONIC_DELETION = 'DEL_PARTIALEXONIC'
        FIVE_PRIME_SPLICE_SITE_DELETION = 'DEL_5PRIMESPLICESITE'
        THREE_PRIME_SPLICE_SITE_DELETION = 'DEL_3PRIMESPLICESITE'
        ALL = [
            EXONIC_DELETION,
            PARTIAL_EXONIC_DELETION,
            FIVE_PRIME_SPLICE_SITE_DELETION,
            THREE_PRIME_SPLICE_SITE_DELETION
        ]

    class DuplicationSubtypes:
        TANDEM_DUPLICATION = 'DUP_TANDEM'
        SEGMENTAL_DUPLICATION = 'DUP_SEG'
        INTERSPERSED_DUPLICATION = 'DUP_INTERSPERSED'
        ALL = [
            TANDEM_DUPLICATION,
            SEGMENTAL_DUPLICATION,
            INTERSPERSED_DUPLICATION
        ]

    class InsertionSubtypes:
        EXONIC_INSERTION = 'INS_EXONIC'
        INTRONIC_INSERTION = 'INS_INTRONIC'
        FIVE_PRIME_SPLICE_SITE_INSERTION = 'INS_5PRIMESPLICESITE'
        THREE_PRIME_SPLICE_SITE_INSERTION = 'INS_3PRIMESPLICESITE'
        ALL = [
            EXONIC_INSERTION,
            INTRONIC_INSERTION,
            FIVE_PRIME_SPLICE_SITE_INSERTION,
            THREE_PRIME_SPLICE_SITE_INSERTION
        ]

    QueryTypeDictionary = {
        SINGLE_NUCLEOTIDE_VARIANT: [SINGLE_NUCLEOTIDE_VARIANT],
        MULTI_NUCLEOTIDE_VARIANT: [MULTI_NUCLEOTIDE_VARIANT],
        INSERTION: [DUPLICATION, INSERTION],
        DELETION: [DELETION],
        INVERSION: [BREAKPOINT, INVERSION, TRANSLOCATION],
        DUPLICATION: [DUPLICATION, INSERTION],
        TRANSLOCATION: [BREAKPOINT, INVERSION, TRANSLOCATION],
        BREAKPOINT: [BREAKPOINT, INVERSION, TRANSLOCATION],
        SPLICE_VARIANT: [SPLICE_VARIANT],
        NOVEL_ISOFORM: [NOVEL_ISOFORM],
        INTRON_RETENTION: [INTRON_RETENTION],
        FUSION_GENE: [FUSION_GENE],
        CIRCULAR_RNA: [CIRCULAR_RNA],
        JET: [JET],
        ENDOGENOUS_RETROVIRUS: [ENDOGENOUS_RETROVIRUS]
    }


class TranslocationOrientations:
    ORIENTATION_1 = 't[p['  # piece extending to the right of p is joined after t
    ORIENTATION_2 = 't]p]'  # reverse complement piece extending left of p is joined after t
    ORIENTATION_3 = ']p]t'  # piece extending to the left of p is joined before t
    ORIENTATION_4 = '[p[t'  # reverse complement extending right of p is joined before t
    ALL = [
        ORIENTATION_1,
        ORIENTATION_2,
        ORIENTATION_3,
        ORIENTATION_4
    ]


class Strands:
    POSITIVE = '+'
    NEGATIVE = '-'
    BOTH_STRANDS = '+-'

