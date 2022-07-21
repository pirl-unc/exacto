#!/usr/bin/python3

"""
The purpose of this python3 script is to implement classes and functions
related to generating a graph genome.

Last updated date: May 23, 2022

Author: Jin Seok (Andy) Lee
"""


import pysam


class GenomicSequence:

    def __init__(self,
                 base: str,
                 genotype: str,
                 ref_chrom: str,
                 ref_pos: int,
                 is_variant: bool,
                 variant_id: str):
        """
        Args
        ----
        base        : Base (A, C, G, T).
        genotype    : '0/0' (homozygous reference),
                      '0/1' (heterozygous) or
                      '1/1' (homozygous variant).
        ref_chrom   : Reference genome chromosome name.
        ref_pos     : Reference genome
        """


class GenomeGraph:

    def __init__(self,
                 vcf_file: str,
                 reference_genome_fasta_file: str,
                 chromosomes: list):
        self.__vcf_file = vcf_file
        self.__reference_genome_fasta_file = reference_genome_fasta_file
        self.__chromosomes = chromosomes
        self.build()


    def build(self):
        # Step 1. Iterate through each chromosome and apply variants
        fasta = pysam.FastaFile(self.__reference_genome_fasta_file)
        for curr_chr in self.__chromosomes:
            curr_chr_length = fasta.get_reference_length(curr_chr)
            print(curr_chr, ':', curr_chr_length)
