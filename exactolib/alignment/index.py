from collections import defaultdict
import json


class Genome:

    def __init__(self):
        self.__genome = []


    def add_contig(self, sequence: str, chrom: str):
        self.__genome.append((chrom, sequence))


    def write(self, output_fasta_file: str):
        # todo write
        a = 1


class Index:

    def __init__(self,
                 kmer_length=25):
        self.__kmer_length = kmer_length

        # key = k-mer, value = list of tuples (chrom, pos)
        self.__index = defaultdict(list)


    def add_sequence(self,
                     sequence: str,
                     chrom: str,
                     start_pos: int):
        """
        Adds a sequence to the index.

        Args
        ----
        sequence        : sequence.
        chrom           : chromosome name.
        start_pos       : start position (i.e. offset).
        """
        # Step 1. Check the length of the sequence
        if len(sequence) < self.__kmer_length:
            print(
                "Sequence length of smaller than the k-mer of the index (",
                self.__kmer_length, ")"
            )
            return

        # Step 2. Add each k-mer of the sequence
        for i in range(0, len(sequence) - self.__kmer_length + 1):
            self.add_kmer(kmer=sequence[i:i + self.__kmer_length],
                          chrom=chrom,
                          start_pos=i + start_pos)


    def add_kmer(self,
                 kmer: str,
                 chrom: str,
                 start_pos: int):
        """
        Adds a k-mer to the index.

        Args
        ----
        kmer        : k-mer sequence.
        position    : position of k-mer in the genome (format: "<contig>:<pos>")
        """
        # Step 1. Check the length of the k-mer sequence
        if len(kmer) != self.__kmer_length:
            print(
                "k-mer length does not match the pre-defined k-mer length (",
                self.__kmer_length, ")"
            )
            return

        # Step 2. Add k-mer to the index
        self.__index[kmer].append((chrom, start_pos))

    


    # def query(self, sequence: str):
    #
    #
    #     return self.__index[kmer]


    def write(self,
              output_txt_file: str):
        with open(output_txt_file, 'w') as f:
            f.write(json.dumps(output_txt_file))
