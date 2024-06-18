import os
import pysam
import random


HG38_FASTA_FILE = "/Users/leework/Documents/Research/projects/seqdata/references/hg38.fa"


def reverse_complement(sequence):
    complement_dict = {'A': 'T',
                       'T': 'A',
                       'C': 'G',
                       'G': 'C',
                       'a': 't',
                       't': 'a',
                       'c': 'g',
                       'g': 'c'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence


if __name__ == "__main__":
    fasta = pysam.FastaFile(HG38_FASTA_FILE)

    # Step 1. Get TMPRSS2-ERG translocation sequence
    tmprss2_sequence = reverse_complement(sequence=fasta.fetch("chr21", 41494453, 41508065))
    erg_sequence = reverse_complement(sequence=fasta.fetch("chr21", 38380027, 38391687))
    sequence_1 = tmprss2_sequence + erg_sequence

    # Step 2. Get EGFR sequence
    egfr_sequence_1 = fasta.fetch("chr7", 55019000, 55020000)
    egfr_sequence_2 = fasta.fetch("chr7", 55050000, 55055000)
    egfr_sequence_3 = fasta.fetch("chr7", 55105000, 55110000)
    sequence_2 = egfr_sequence_1 + egfr_sequence_2 + egfr_sequence_3

    # Step 3. Write to FASTQ file
    output_file = '../../../test/data/fastq/hg38_tumor_long_read_dna_1.fastq'
    with open(output_file, 'w') as f:
        for i in range(0, 60):
            read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                random.randint(100000,999999),
                                                i+1)
            base_quality_scores = [chr(96) for _ in range(0,len(sequence_1))]
            base_quality_scores = ''.join(base_quality_scores)
            f.write(read_id + '\n')
            f.write(sequence_1 + '\n')
            f.write('+\n')
            f.write(base_quality_scores + '\n')
        for i in range(0, 60):
            read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                random.randint(100000,999999),
                                                i+1)
            base_quality_scores = [chr(96) for _ in range(0,len(sequence_2))]
            base_quality_scores = ''.join(base_quality_scores)
            f.write(read_id + '\n')
            f.write(sequence_2 + '\n')
            f.write('+\n')
            f.write(base_quality_scores + '\n')
    os.system('gzip %s' % output_file)


