import os
import random
from typing import List


def create_fastq_file(
        sequences: List[str],
        output_fastq_file: str,
        num_reads: List[int],
        stranded: bool = False
):
    """
    Create a FASTQ (gzipped) file.

    Parameters:
        sequences               :   List of sequences.
        output_fastq_file       :   Output FASTQ file.
        num_reads               :   Number of reads.
        stranded                :   If True, then only the original sequence
                                    is generated for the FASTQ. If False, at 0.5 probability,
                                    the original sequence and the reverse complemented sequence is generated.
    """
    with open(output_fastq_file, 'w') as f:
        for i,sequence in enumerate(sequences):
            for j in range(0, num_reads[i]):
                if stranded == False:
                    if random.random() < 0.5:
                        sequence_ = sequence
                    else:
                        sequence_ = reverse_complement(sequence)
                else:
                    sequence_ = sequence

                read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                    random.randint(100000,999999),
                                                    j+1)
                base_quality_scores = [chr(96) for _ in range(0,len(sequence_))]
                base_quality_scores = ''.join(base_quality_scores)
                f.write(read_id + '\n')
                f.write(sequence_ + '\n')
                f.write('+\n')
                f.write(base_quality_scores + '\n')
    os.system('gzip %s' % output_fastq_file)


def reverse_complement(sequence: str) -> str:
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'N': 'N',
        'n': 'n'
    }
    reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
    return reverse_complement_seq

