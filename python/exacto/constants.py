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

    class DuplicationSubtypes:
        TANDEM_DUPLICATION = 'DUP_TANDEM'
        SEGMENTAL_DUPLICATION = 'DUP_SEG'
        ALL = [
            TANDEM_DUPLICATION,
            SEGMENTAL_DUPLICATION
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
    GATK4_MUTECT2 = "gatk4-mutect2"
    DEEPVARIANT = "deepvariant"
    STRELKA2_GERMLINE = "strelka2-germline"
    STRELKA2_SOMATIC = "strelka2-somatic"
    SNIFFLES2 = 'sniffles2'
    SVIM = 'svim'
    CUTESV = 'cutesv'
    DELLY2 = 'delly2'
    LUMPY = 'lumpy'
    PBSV = 'pbsv'
    ALL = [
        GATK4_MUTECT2,
        DEEPVARIANT,
        STRELKA2_GERMLINE,
        STRELKA2_SOMATIC,
        SNIFFLES2,
        SVIM,
        CUTESV,
        DELLY2,
        LUMPY,
        PBSV
    ]

    class AttributeTypes:
        CUTESV = {
            'cutesv_id': str,
            'cutesv_svtype': str,
            'cutesv_svlen': int,
            'cutesv_chr2': str,
            'cutesv_end': int,
            'cutesv_cipos': str,
            'cutesv_cilen': str,
            'cutesv_re': int,
            'cutesv_strand': str,
            'cutesv_rnames': str,
            'cutesv_af': float,
            'cutesv_precise': bool,
            'cutesv_gt': str,
            'cutesv_gq': float,
            'cutesv_pl': str,
            'cutesv_dr': int,
            'cutesv_dv': int
        }
        DEEPVARIANT = {
            'deepvariant_id': str,
            'deepvariant_end': int,
            'deepvariant_gt': str,
            'deepvariant_gq': int,
            'deepvariant_dp': int,
            'deepvariant_min_dp': int,
            'deepvariant_ad': str,
            'deepvariant_vaf': float,
            'deepvariant_pl': str,
            'deepvariant_med_dp': int
        }
        GATK4_MUTECT2 = {
            'gatk4-mutect2_id': str,
            'gatk4-mutect2_as_filterstatus': str,
            'gatk4-mutect2_as_sb_table': str,
            'gatk4-mutect2_as_uniq_alt_read_count': int,
            'gatk4-mutect2_contq': float,
            'gatk4-mutect2_dp': int,
            'gatk4-mutect2_ecnt': int,
            'gatk4-mutect2_germq': int,
            'gatk4-mutect2_mbq': int,
            'gatk4-mutect2_mfrl': int,
            'gatk4-mutect2_mmq': int,
            'gatk4-mutect2_mpos': int,
            'gatk4-mutect2_nalod': float,
            'gatk4-mutect2_ncount': int,
            'gatk4-mutect2_nlod': float,
            'gatk4-mutect2_ocm': int,
            'gatk4-mutect2_pon': bool,
            'gatk4-mutect2_popaf': float,
            'gatk4-mutect2_af': float,
            'gatk4-mutect2_roq': float,
            'gatk4-mutect2_rpa': int,
            'gatk4-mutect2_ru': str,
            'gatk4-mutect2_seqq': int,
            'gatk4-mutect2_str': bool,
            'gatk4-mutect2_strandq': int,
            'gatk4-mutect2_strq': int,
            'gatk4-mutect2_tlod': float,
            'gatk4-mutect2_ad': str,
            'gatk4-mutect2_dp': int,
            'gatk4-mutect2_f1r2': str,
            'gatk4-mutect2_f2r1': str,
            'gatk4-mutect2_fad': str,
            'gatk4-mutect2_gq': float,
            'gatk4-mutect2_gt': str,
            'gatk4-mutect2_pgt': str,
            'gatk4-mutect2_pid': str,
            'gatk4-mutect2_pl': int,
            'gatk4-mutect2_ps': int,
            'gatk4-mutect2_sb': str
        }
        PBSV = {
            'pbsv_id': str,
            'pbsv_svtype': str,
            'pbsv_end': int,
            'pbsv_svlen': int,
            'pbsv_svann': str,
            'pbsv_cipos': str,
            'pbsv_mateid': str,
            'pbsv_matedist': int,
            'pbsv_precise': bool,
            'pbsv_gt': str,
            'pbsv_dp': int,
            'pbsv_ad': str,
            'pbsv_sac': str
        }
        SNIFFLES2 = {
            'sniffles2_id': str,
            'sniffles2_svlen': int,
            'sniffles2_svtype': str,
            'sniffles2_chr2': str,
            'sniffles2_support': int,
            'sniffles2_support_inline': int,
            'sniffles2_support_long': int,
            'sniffles2_end': int,
            'sniffles2_stdev_pos': float,
            'sniffles2_stdev_len': float,
            'sniffles2_coverage': str,
            'sniffles2_strand': str,
            'sniffles2_ac': int,
            'sniffles2_supp_vec': str,
            'sniffles2_consensus_support': int,
            'sniffles2_rnames': str,
            'sniffles2_af': float,
            'sniffles2_nm': float,
            'sniffles2_phase': str,
            'sniffles2_gt': str,
            'sniffles2_gq': int,
            'sniffles2_dr': int,
            'sniffles2_dv': int,
            'sniffles2_precise': bool
        }
        STRELKA2_SOMATIC = {
            'strelka2-somatic_id': str,
            'strelka2-somatic_qss': int,
            'strelka2-somatic_tqss': int,
            'strelka2-somatic_nt': str,
            'strelka2-somatic_qss_nt': int,
            'strelka2-somatic_tqss_nt': int,
            'strelka2-somatic_sgt': str,
            'strelka2-somatic_dp': int,
            'strelka2-somatic_mq': float,
            'strelka2-somatic_mq0': int,
            'strelka2-somatic_readposranksum': float,
            'strelka2-somatic_pnoise': float,
            'strelka2-somatic_pnoise2': float,
            'strelka2-somatic_somaticevs': float,
            'strelka2-somatic_qsi': int,
            'strelka2-somatic_tqsi': int,
            'strelka2-somatic_qsi_nt': int,
            'strelka2-somatic_tqsi_nt': int,
            'strelka2-somatic_ru': str,
            'strelka2-somatic_rc': int,
            'strelka2-somatic_ic': int,
            'strelka2-somatic_ihp': int,
            'strelka2-somatic_somatic': bool,
            'strelka2-somatic_overlap': bool,
            'strelka2-somatic_fdp': int,
            'strelka2-somatic_sdp': int,
            'strelka2-somatic_subdp': int,
            'strelka2-somatic_au': str,
            'strelka2-somatic_cu': str,
            'strelka2-somatic_gu': str,
            'strelka2-somatic_tu': str,
            'strelka2-somatic_dp': int,
            'strelka2-somatic_dp2': int,
            'strelka2-somatic_tar': int,
            'strelka2-somatic_tir': int,
            'strelka2-somatic_tor': int,
            'strelka2-somatic_snvsb': float,
            'strelka2-somatic_dp50': float,
            'strelka2-somatic_fdp50': float,
            'strelka2-somatic_subdp50': float,
            'strelka2-somatic_bcn50': float
        }
        STRELKA2_GERMLINE = {
            'strelka2-germline_id': str,
            'strelka2-germline_end': int,
            'strelka2-germline_snvhpol': int,
            'strelka2-germline_cigar': str,
            'strelka2-germline_ru': str,
            'strelka2-germline_refrep': int,
            'strelka2-germline_idrep': int,
            'strelka2-germline_mq': int,
            'strelka2-germline_blockavg_min30p3a': bool,
            'strelka2-germline_gt': str,
            'strelka2-germline_gq': int,
            'strelka2-germline_gqx': int,
            'strelka2-germline_dp': int,
            'strelka2-germline_dpf': int,
            'strelka2-germline_min_dp': int,
            'strelka2-germline_ad': str,
            'strelka2-germline_adf': str,
            'strelka2-germline_adr': str,
            'strelka2-germline_ft': str,
            'strelka2-germline_dpi': int,
            'strelka2-germline_pl': int,
            'strelka2-germline_ps': int,
            'strelka2-germline_sb': float
        }
        SVIM = {
            'svim_id': str,
            'svim_svtype': str,
            'svim_end': int,
            'svim_svlen': int,
            'svim_support': int,
            'svim_std_span': float,
            'svim_std_pos': float,
            'svim_std_pos1': float,
            'svim_std_pos2': float,
            'svim_seqs': str,
            'svim_reads': str,
            'svim_cutpaste': bool,
            'svim_gt': bool,
            'svim_dp': int,
            'svim_ad': str,
            'svim_cn': int
        }

class GenomicRegionTypes:
    EXONIC = 'exonic'
    INTRONIC = 'intronic'
    INTERGENIC = 'intergenic'
    ALL = [
        EXONIC,
        INTRONIC,
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


class GeneTypes:
    PROTEIN_CODING = 'protein_coding'


class TranscriptTypes:
    PROTEIN_CODING = 'protein_coding'


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


