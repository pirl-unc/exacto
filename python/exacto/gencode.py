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
The purpose of this python3 script is to implement the Gencode dataclass.
"""


import pysam
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Tuple
from .annotation_db import AnnotationDb
from .constants import *
from .exon import Exon
from .fasta import Fasta
from .gene import Gene
from .gene_set import GeneSet
from .logging import get_logger
from .nucleotide_sequence import NucleotideSequence
from .transcript import Transcript
from .variant import VariantAnnotation
from .variant_call import VariantCall
from .variants_list import VariantsList


logger = get_logger(__name__)


@dataclass
class Gencode(AnnotationDb):
    gtf_file: str
    genome_fasta_file: str
    version: str
    genome: str
    gene_set: GeneSet = GeneSet()
    genome_fasta: pysam.FastaFile = None

    @property
    def source(self):
        return AnnotationSources.GENCODE

    def __post_init__(self):
        # Step 1. Read FASTA file
        self.genome_fasta = pysam.FastaFile(self.genome_fasta_file)

        # Step 2. Read GTF file
        self.__read_gtf_file_genes()  # add genes
        self.__read_gtf_file_transcripts()  # add transcripts
        self.__read_gtf_file_exons()  # add exons
        self.__read_gtf_file_start_codons()  # update start codon start and end positions
        self.__read_gtf_file_stop_codons()  # update stop codon start and end positions
        self.__read_gtf_file_utr()  # update UTR start and end positions

    @staticmethod
    def get_stable_ensembl_id(id: str):
        """
        Returns the stable Ensembl ID and version.

        Parameters
        ----------
        id          :   Ensembl ID (e.g. 'ENSG00001.1').

        Returns
        -------
        stable_id   :   Stable Ensembl ID (without version number).
        version     :   Version.
        """
        if (('ENSG' in id) or ('ENST' in id or 'ENSE' in id)) and '.' in id:
            stable_id = id.split('.')[0]
            version = id.split('.')[1]
            return stable_id, version
        else:
            return id, None

    def __read_gtf_file_genes(self):
        """
        Reads and loads all 'gene' rows in GENCODE GTF file.
        """
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'gene':
                    curr_gene_chrom = str(elements[0])
                    curr_gene_source = str(elements[1])
                    curr_gene_start = int(elements[3])
                    curr_gene_end = int(elements[4])
                    curr_gene_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_gene_stable_id, curr_gene_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['gene_id']))
                    curr_gene_name = str(curr_metadata_dict['gene_name'])
                    curr_gene_type = str(curr_metadata_dict['gene_type'])
                    curr_gene_level = int(curr_metadata_dict['level'])
                    gene = Gene(
                        id=curr_gene_id,
                        stable_id=curr_gene_stable_id,
                        source=curr_gene_source,
                        source_version=self.version,
                        name=curr_gene_name,
                        chromosome=curr_gene_chrom,
                        start=curr_gene_start,
                        end=curr_gene_end,
                        strand=curr_gene_strand,
                        type=curr_gene_type,
                        level=curr_gene_level,
                        version=curr_gene_version,
                        genome=self.genome
                    )
                    self.gene_set.add_gene(gene=gene)
        logger.info('Loaded %i genes in total.' % len(self.gene_set.gene_ids))

    def __read_gtf_file_transcripts(self):
        """
        Reads and loads all 'transcript' rows in GENCODE GTF file.
        """
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'transcript':
                    curr_transcript_chrom = str(elements[0])
                    curr_transcript_source = str(elements[1])
                    curr_transcript_start = int(elements[3])
                    curr_transcript_end = int(elements[4])
                    curr_transcript_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_transcript_stable_id, curr_transcript_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['transcript_id'][0]))
                    curr_transcript_type = str(curr_metadata_dict['transcript_type'][0])
                    curr_transcript_name = str(curr_metadata_dict['transcript_name'][0])
                    curr_transcript_tags = [str(tag).replace('"', '') for tag in curr_metadata_dict['tag']]
                    try:
                        curr_transcript_level = int(curr_metadata_dict['level'][0])
                    except:
                        curr_transcript_level = None
                    try:
                        curr_transcript_support_level = int(curr_metadata_dict['transcript_support_level'][0])
                    except:
                        curr_transcript_support_level = None
                    transcript = Transcript(
                        id=curr_transcript_id,
                        stable_id=curr_transcript_stable_id,
                        source=curr_transcript_source,
                        source_version=self.version,
                        chromosome=curr_transcript_chrom,
                        start=curr_transcript_start,
                        end=curr_transcript_end,
                        type=curr_transcript_type,
                        strand=curr_transcript_strand,
                        version=curr_transcript_version,
                        name=curr_transcript_name,
                        level=curr_transcript_level,
                        support_level=curr_transcript_support_level,
                        tags=curr_transcript_tags,
                        genome=self.genome
                    )
                    self.gene_set.add_transcript(gene_id=curr_gene_id, transcript=transcript)

        transcripts_count = 0
        for gene in self.gene_set.genes:
            transcripts_count += gene.transcripts_count
        logger.info('Loaded %i transcripts in total.' % transcripts_count)

    def __read_gtf_file_exons(self):
        """
        Reads and loads all 'exon' rows in GENCODE GTF file.
        """
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'exon':
                    curr_exon_chrom = str(elements[0])
                    curr_exon_source = str(elements[1])
                    curr_exon_start = int(elements[3])
                    curr_exon_end = int(elements[4])
                    curr_exon_strand = str(elements[6])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = defaultdict(list)
                    curr_exon_tags = []
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
                    curr_gene_id = str(curr_metadata_dict['gene_id'][0])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
                    curr_exon_id = str(curr_metadata_dict['exon_id'][0])
                    curr_exon_stable_id, curr_exon_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['exon_id'][0]))
                    curr_exon_number = int(curr_metadata_dict['exon_number'][0])
                    curr_exon_tags = [str(tag).replace('"', '') for tag in curr_metadata_dict['tag']]
                    if curr_exon_strand == Strands.POSITIVE:
                        sequence = self.genome_fasta.fetch(curr_exon_chrom, curr_exon_start - 1, curr_exon_end)
                    else:
                        sequence = self.genome_fasta.fetch(curr_exon_chrom, curr_exon_start - 1, curr_exon_end)
                        sequence = NucleotideSequence.reverse_complement(sequence=sequence)
                    exon = Exon(
                        id=curr_exon_id,
                        stable_id=curr_exon_stable_id,
                        source=curr_exon_source,
                        source_version=self.version,
                        chromosome=curr_exon_chrom,
                        start=curr_exon_start,
                        end=curr_exon_end,
                        number=curr_exon_number,
                        strand=curr_exon_strand,
                        sequence=NucleotideSequence(sequence=sequence),
                        version=curr_exon_version,
                        tags=curr_exon_tags
                    )
                    self.gene_set.add_exon(
                        gene_id=curr_gene_id,
                        transcript_id=curr_transcript_id,
                        exon=exon
                    )
        exons_count = 0
        for gene in self.gene_set.genes:
            for transcript in gene.transcripts:
                exons_count += transcript.exons_count
        logger.info('Loaded %i exons in total.' % exons_count)

    def __read_gtf_file_start_codons(self):
        """
        Reads and loads all 'start_codon' rows in GENCODE GTF file.
        """
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'start_codon':
                    curr_start_codon_start = int(elements[3])
                    curr_start_codon_end = int(elements[4])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].start_codon_start = curr_start_codon_start
                    self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].start_codon_end = curr_start_codon_end

    def __read_gtf_file_stop_codons(self):
        """
        Reads and loads all 'stop_codon' rows in GENCODE GTF file.
        """
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'stop_codon':
                    curr_stop_codon_start = int(elements[3])
                    curr_stop_codon_end = int(elements[4])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')

                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].stop_codon_start = curr_stop_codon_start
                    self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].stop_codon_end = curr_stop_codon_end

    def __read_gtf_file_utr(self):
        """
        Reads and loads all 'UTR' rows in GENCODE GTF file.
        """
        with open(self.gtf_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0:2] == '##':
                    continue
                elements = line.split('\t')
                if elements[2] == 'UTR':
                    curr_utr_start = int(elements[3])
                    curr_utr_end = int(elements[4])
                    curr_metadata = str(elements[8]).split(';')
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')

                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].utr_start = curr_utr_start
                    self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].utr_end = curr_utr_end

    def annotate_position(
            self,
            chromosome: str,
            position: int
    ) -> List[VariantAnnotation]:
        """
        Annotates a position and returns a list of VariantAnnotation objects.

        Parameters
        ----------
        chromosome              :   Chromosome.
        position                :   Position.

        Returns
        -------
        variant_annotations     :   List of VariantAnnotation objects.
        """
        variant_annotations = []
        genes = self.gene_set.find_genes(chromosome=chromosome, position=position)
        if len(genes) == 0:
            variant_annotation = VariantAnnotation(
                region=GenomicRegionTypes.INTERGENIC,
                source=AnnotationSources.ENSEMBL,
                source_version=self.version
            )
            variant_annotations.append(variant_annotation)
        else:
            for gene in genes:
                region = GenomicRegionTypes.INTRONIC
                for transcript in gene.transcripts:
                    for exon in transcript.exons:
                        if exon.start <= position <= exon.end:
                            region = GenomicRegionTypes.EXONIC
                variant_annotation = VariantAnnotation(
                    region=region,
                    source=AnnotationSources.GENCODE,
                    source_version=self.version,
                    gene_id=gene.id,
                    gene_stable_id=gene.stable_id,
                    gene_version=gene.version,
                    gene_name=gene.name,
                    gene_type=gene.type,
                    gene_strand=gene.strand,
                    genome=gene.genome
                )
                variant_annotations.append(variant_annotation)
        return variant_annotations

    def annotate_variant_call(
            self,
            variant_call: VariantCall
    ) -> Tuple[List[VariantAnnotation], List[VariantAnnotation]]:
        """
        Annotates a VariantCall object and returns two lists of VariantAnnotation objects.

        Parameters
        ----------
        variant_call        :   VariantCall object.

        Returns
        -------
        pos_1_annotations   :   List of VariantAnnotation objects for position 1.
        pos_2_annotations   :   List of VariantAnnotation objects for position 2.
        """
        pos_1_annotations = self.annotate_position(
            chromosome=variant_call.chromosome_1,
            position=variant_call.position_1
        )
        pos_2_annotations = self.annotate_position(
            chromosome=variant_call.chromosome_2,
            position=variant_call.position_2
        )
        return pos_1_annotations, pos_2_annotations

    def annotate_variants_list(self, variants_list) -> VariantsList:
        """
        Annotates a VariantsList object.

        Parameters
        ----------
        variants_list   :   VariantsList object.

        Returns
        -------
        variants_list   :   VariantsList object with each VariantCall object annotated.
        """
        if self.gene_set is None:
            raise Exception('Please read a GENCODE comprehensive gene annotation '
                            'GTF file before annotating a VariantsList object.')
        for i in range(0, variants_list.size):
            for j in range(0, variants_list.variants[i].size):
                position_1_annotations, position_2_annotations = self.annotate_variant_call(
                    variants_list.variants[i].variant_calls[j]
                )
                for annotation in position_1_annotations:
                    variants_list.variants[i].variant_calls[j].position_1_annotations.append(annotation)
                for annotation in position_2_annotations:
                    variants_list.variants[i].variant_calls[j].position_2_annotations.append(annotation)
        return variants_list
