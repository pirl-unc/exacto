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


DEFAULT_ATTRIBUTE_VALUE = ''


class NucleicAcidTypes:
    DNA = 'dna'
    RNA = 'rna'
    DNA_RNA = 'dna_rna'
    ALL = [
        DNA,
        RNA,
        DNA_RNA
    ]


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
        SINGLE_NUCLEOTIDE_VARIANT: 'SNV',
        INSERTION: 'Insertion',
        DELETION: 'Deletion',
        DUPLICATION: 'Duplication',
        INVERSION: 'Inversion',
        TRANSLOCATION: 'Translocation'
    }

    class DuplicationSubtypes:
        TANDEM_DUPLICATION = 'DUP_TANDEM'
        SEGMENTAL_DUPLICATION = 'DUP_SEG'
        INTERSPERSED_DUPLICATION = 'DUP_INTERSPERSED'
        ALL = [
            TANDEM_DUPLICATION,
            SEGMENTAL_DUPLICATION,
            INTERSPERSED_DUPLICATION
        ]

    QueryTypeDictionary = {
        SINGLE_NUCLEOTIDE_VARIANT: [SINGLE_NUCLEOTIDE_VARIANT],
        MULTI_NUCLEOTIDE_VARIANT: [MULTI_NUCLEOTIDE_VARIANT],
        INSERTION: [INSERTION, DUPLICATION],
        DELETION: [DELETION],
        INVERSION: [INVERSION, BREAKPOINT],
        DUPLICATION: [DUPLICATION, INSERTION],
        TRANSLOCATION: [TRANSLOCATION, BREAKPOINT, INVERSION],
        BREAKPOINT: [BREAKPOINT, INVERSION, TRANSLOCATION],
        SPLICE_VARIANT: [SPLICE_VARIANT],
        NOVEL_ISOFORM: [NOVEL_ISOFORM],
        INTRON_RETENTION: [INTRON_RETENTION],
        FUSION_GENE: [FUSION_GENE],
        CIRCULAR_RNA: [CIRCULAR_RNA],
        JET: [JET],
        ENDOGENOUS_RETROVIRUS: [ENDOGENOUS_RETROVIRUS]
    }


class VariantCallingMethods:
    GATK4_MUTECT2 = 'gatk4-mutect2'
    DEEPVARIANT = 'deepvariant'
    EXACTO = 'exacto'
    STRELKA2 = 'strelka2'
    SNIFFLES2 = 'sniffles2'
    SVIM = 'svim'
    CUTESV = 'cutesv'
    DELLY2 = 'delly2'
    LUMPY = 'lumpy'
    PBSV = 'pbsv'
    ALL = [
        GATK4_MUTECT2,
        DEEPVARIANT,
        STRELKA2,
        SNIFFLES2,
        SVIM,
        CUTESV,
        DELLY2,
        LUMPY,
        PBSV
    ]

    class AttributeTypes:
        CUTESV = {
            'id': str,
            'svtype': str,
            'svlen': int,
            'chr2': str,
            'end': int,
            'cipos': str,
            'cilen': str,
            're': int,
            'strand': str,
            'rnames': str,
            'af': float,
            'precise': bool,
            'gt': str,
            'gq': float,
            'pl': str,
            'dr': int,
            'dv': int
        }
        DEEPVARIANT = {
            'id': str,
            'end': int,
            'gt': str,
            'gq': int,
            'dp': int,
            'min_dp': int,
            'ad': str,
            'vaf': float,
            'pl': str,
            'med_dp': int
        }
        GATK4_MUTECT2 = {
            'id': str,
            'as_filterstatus': str,
            'as_sb_table': str,
            'as_uniq_alt_read_count': int,
            'contq': float,
            'ecnt': int,
            'germq': int,
            'mbq': int,
            'mfrl': int,
            'mmq': int,
            'mpos': int,
            'nalod': float,
            'ncount': int,
            'nlod': float,
            'ocm': int,
            'pon': bool,
            'popaf': float,
            'af': float,
            'roq': float,
            'rpa': int,
            'ru': str,
            'seqq': int,
            'str': bool,
            'strandq': int,
            'strq': int,
            'tlod': float,
            'ad': str,
            'dp': int,
            'f1r2': str,
            'f2r1': str,
            'fad': str,
            'gq': float,
            'gt': str,
            'pgt': str,
            'pid': str,
            'pl': int,
            'ps': int,
            'sb': str
        }
        PBSV = {
            'id': str,
            'svtype': str,
            'end': int,
            'svlen': int,
            'svann': str,
            'cipos': str,
            'mateid': str,
            'matedist': int,
            'precise': bool,
            'gt': str,
            'dp': int,
            'ad': str,
            'sac': str
        }
        SNIFFLES2 = {
            'id': str,
            'svlen': int,
            'svtype': str,
            'chr2': str,
            'support': int,
            'support_inline': int,
            'support_long': int,
            'end': int,
            'stdev_pos': float,
            'stdev_len': float,
            'coverage': str,
            'strand': str,
            'ac': int,
            'supp_vec': str,
            'consensus_support': int,
            'rnames': str,
            'af': float,
            'nm': float,
            'phase': str,
            'gt': str,
            'gq': int,
            'dr': int,
            'dv': int,
            'precise': bool
        }
        STRELKA2 = {
            'id': str,
            'qss': int,
            'tqss': int,
            'nt': str,
            'qss_nt': int,
            'tqss_nt': int,
            'sgt': str,
            'mq': float,
            'mq0': int,
            'readposranksum': float,
            'pnoise': float,
            'pnoise2': float,
            'somaticevs': float,
            'qsi': int,
            'tqsi': int,
            'qsi_nt': int,
            'tqsi_nt': int,
            'ru': str,
            'rc': int,
            'ic': int,
            'ihp': int,
            'somatic': bool,
            'overlap': bool,
            'fdp': int,
            'sdp': int,
            'subdp': int,
            'au': str,
            'cu': str,
            'gu': str,
            'tu': str,
            'dp': int,
            'dp2': int,
            'tar': int,
            'tir': int,
            'tor': int,
            'snvsb': float,
            'dp50': float,
            'fdp50': float,
            'subdp50': float,
            'bcn50': float,
            'end': int,
            'snvhpol': int,
            'cigar': str,
            'refrep': int,
            'idrep': int,
            'blockavg_min30p3a': bool,
            'gt': str,
            'gq': int,
            'gqx': int,
            'dpf': int,
            'min_dp': int,
            'ad': str,
            'adf': str,
            'adr': str,
            'ft': str,
            'dpi': int,
            'pl': int,
            'ps': int,
            'sb': float

        }
        SVIM = {
            'id': str,
            'svtype': str,
            'end': int,
            'svlen': int,
            'support': int,
            'std_span': float,
            'std_pos': float,
            'std_pos1': float,
            'std_pos2': float,
            'zmws': int,
            'seqs': str,
            'reads': str,
            'cutpaste': bool,
            'gt': bool,
            'dp': int,
            'ad': str,
            'cn': int
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


class GenomicRegionTypes:
    EXONIC = 'exonic'
    INTRONIC = 'intronic'
    FIVE_PRIME_UTR = '5prime_utr'
    THREE_PRIME_UTR = '3prime_utr'
    INTERGENIC = 'intergenic'
    ALL = [
        EXONIC,
        INTRONIC,
        FIVE_PRIME_UTR,
        THREE_PRIME_UTR,
        INTERGENIC
    ]


class FunctionalConsequenceTypes:
    FIVE_PRIME_SPLICE_SITE = '5prime_splice_site'
    THREE_PRIME_SPLICE_SITE = '3prime_splice_site'
    ALL = [
        FIVE_PRIME_SPLICE_SITE,
        THREE_PRIME_SPLICE_SITE
    ]


class AnnotationSources:
    ENSEMBL = 'ensembl'
    GENCODE = 'gencode'
    ANNOVAR = 'annovar'
    ALL = [
        ENSEMBL,
        GENCODE,
        ANNOVAR
    ]


class Sexes:
    MALE = 'male'
    FEMALE = 'female'
    ALL = [
        MALE,
        FEMALE
    ]


class Strands:
    POSITIVE = '+'
    NEGATIVE = '-'
    BOTH_STRANDS = '+-'


class SoftclipDirections:
    TOWARDS_THREE_PRIME = '3prime'
    TOWWARDS_FIVE_PRIME = '5prime'


class VariantFilterSampleTypes:
    CASE = 'case'
    CONTROL = 'control'


class VariantFilterQuantifiers:
    ALL = 'all'
    ANY = 'any'
    MEDIAN = 'median'
    AVERAGE = 'average'
    MIN = 'min'
    MAX = 'max'


class VariantFilterOperators:
    LESS_THAN = '<'
    LESS_THAN_OR_EQUAL_TO = '<='
    GREATER_THAN = '>'
    GREATER_THAN_OR_EQUAL_TO = '>='
    EQUALS = '=='
    IN = 'in'


RNA_CODONS = {
    # U
    'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',  # UxU
    'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',  # UxC
    'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',  # UxA
    'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',  # UxG

    # C
    'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',  # CxU
    'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',  # CxC
    'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',  # CxA
    'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',  # CxG

    # A
    'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',  # AxU
    'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',  # AxC
    'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',  # AxA
    'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',  # AxG

    # G
    'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',  # GxU
    'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',  # GxC
    'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',  # GxA
    'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'   # GxG
}

AMINO_ACID_CODES = {
    'Cys': 'C',
    'Asp': 'D',
    'Ser': 'S',
    'Gln': 'Q',
    'Lys': 'K',
    'Trp': 'W',
    'Asn': 'N',
    'Pro': 'P',
    'Thr': 'T',
    'Phe': 'F',
    'Ala': 'A',
    'Gly': 'G',
    'Ile': 'I',
    'Leu': 'L',
    'His': 'H',
    'Arg': 'R',
    'Met': 'M',
    'Val': 'V',
    'Glu': 'E',
    'Tyr': 'Y',
    '---': '*'
}

AMINO_ACID_THREE_LETTER_CODES = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    '*': '---'
}