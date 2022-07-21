#!/usr/bin/python3

"""
The purpose of this python3 script is to define constants.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


class Constants:

    class SequencingPlatforms:
        ILLUMINA = 'illumina'
        PACBIO_HIFI_CCS = 'pacbio-hifi-ccs'
        OXFORD_NANOPORE_TECHNOLOGIES = 'ont'
        ALL = [ILLUMINA,
               PACBIO_HIFI_CCS,
               OXFORD_NANOPORE_TECHNOLOGIES]

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