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
The purpose of this python3 script is to implement the Gencode class.
"""


import pandas as pd
from copy import deepcopy
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .constants import *
from .exon import Exon
from .transcript import Transcript
from .gene import Gene
from .gene_set import GeneSet
from .annotation import Annotation
from .variant import Variant
from .variant_call import VariantCall
from .variants_list import VariantsList
from .variant_annotation import VariantAnnotation
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class Gencode(Annotation):
    gene_set: GeneSet = None
    comprehensive_gene_annotation_gtf_file: str = None
    transcripts_fasta_file: str = None

    @staticmethod
    def get_stable_ensembl_id(id):
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
        with open(self.comprehensive_gene_annotation_gtf_file, 'r') as f:
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
                    if 'ENSG' in curr_gene_id and '.' in curr_gene_id:
                        curr_gene_version = curr_gene_id.split('.')[1]
                        curr_gene_id = curr_gene_id.split('.')[0]
                    curr_gene_name = str(curr_metadata_dict['gene_name'])
                    curr_gene_type = str(curr_metadata_dict['gene_type'])
                    curr_gene_level = int(curr_metadata_dict['level'])
                    gene = Gene(
                        id=curr_gene_id,
                        source=curr_gene_source,
                        name=curr_gene_name,
                        chromosome=curr_gene_chrom,
                        start=curr_gene_start,
                        end=curr_gene_end,
                        strand=curr_gene_strand,
                        type=curr_gene_type,
                        level=curr_gene_level,
                        version=curr_gene_version
                    )
                    self.gene_set.add_gene(gene=gene)
        logger.info('Loaded %i genes in total.' % len(self.gene_set.gene_ids))

    def __read_gtf_file_transcripts(self):
        """
        Reads and loads all 'transcript' rows in GENCODE GTF file.
        """
        with open(self.comprehensive_gene_annotation_gtf_file, 'r') as f:
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
                    curr_metadata_dict = {}
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    if 'ENSG' in curr_gene_id and '.' in curr_gene_id:
                        curr_gene_version = curr_gene_id.split('.')[1]
                        curr_gene_id = curr_gene_id.split('.')[0]
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    if 'ENST' in curr_transcript_id and '.' in curr_transcript_id:
                        curr_transcript_version = curr_transcript_id.split('.')[1]
                        curr_transcript_id = curr_transcript_id.split('.')[0]
                    curr_transcript_type = str(curr_metadata_dict['transcript_type'])
                    try:
                        curr_transcript_level = int(curr_metadata_dict['transcript_support_level'])
                    except:
                        curr_transcript_level = None
                    curr_transcript_name = str(curr_metadata_dict['transcript_name'])
                    transcript = Transcript(
                        id=curr_transcript_id,
                        chromosome=curr_transcript_chrom,
                        start=curr_transcript_start,
                        end=curr_transcript_end,
                        type=curr_transcript_type,
                        strand=curr_transcript_strand,
                        source=curr_transcript_source,
                        version=curr_transcript_version,
                        name=curr_transcript_name,
                        level=curr_transcript_level,

                    )
                    self.gene_set.add_transcript(gene_id=curr_gene_id, transcript=transcript)

        transcripts_count = 0
        for gene in self.gene_set.genes.values():
            transcripts_count += len(gene.transcript_ids)
        logger.info('Loaded %i transcripts in total.' % transcripts_count)

    def __read_gtf_file_exons(self):
        """
        Reads and loads all 'exon' rows in GENCODE GTF file.
        """
        with open(self.comprehensive_gene_annotation_gtf_file, 'r') as f:
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
                    curr_metadata_dict = {}
                    curr_exon_tags = []
                    for curr_metadata_elements in curr_metadata:
                        if curr_metadata_elements == '':
                            continue
                        if curr_metadata_elements[0] == ' ':
                            curr_metadata_elements = curr_metadata_elements[1:]
                        curr_metadata_elements_ = curr_metadata_elements.split(' ')
                        if curr_metadata_elements_[0] == 'tag':
                            curr_exon_tags.append(str(curr_metadata_elements_[1].replace('"', '')))
                        else:
                            curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
                    curr_gene_id = str(curr_metadata_dict['gene_id'])
                    if 'ENSG' in curr_gene_id and '.' in curr_gene_id:
                        curr_gene_version = curr_gene_id.split('.')[1]
                        curr_gene_id = curr_gene_id.split('.')[0]
                    curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                    if 'ENST' in curr_transcript_id and '.' in curr_transcript_id:
                        curr_transcript_version = curr_transcript_id.split('.')[1]
                        curr_transcript_id = curr_transcript_id.split('.')[0]
                    curr_exon_id = str(curr_metadata_dict['exon_id'])
                    if 'ENSE' in curr_exon_id and '.' in curr_exon_id:
                        curr_exon_version = curr_exon_id.split('.')[1]
                        curr_exon_id = curr_exon_id.split('.')[0]
                    curr_exon_number = int(curr_metadata_dict['exon_number'])
                    exon = Exon(
                        id=curr_exon_id,
                        chromosome=curr_exon_chrom,
                        start=curr_exon_start,
                        end=curr_exon_end,
                        sequence='',
                        number=curr_exon_number,
                        strand=curr_exon_strand,
                        version=curr_exon_version,
                        source=curr_exon_source,
                        tags=curr_exon_tags
                    )
                    self.gene_set.add_exon(
                        gene_id=curr_gene_id,
                        transcript_id=curr_transcript_id,
                        exon=exon
                    )
        exons_count = 0
        for gene in self.gene_set.genes.values():
            for transcript in gene.transcripts.values():
                exons_count += len(transcript.exon_ids)
        logger.info('Loaded %i exons in total.' % exons_count)

    def __read_gtf_file_start_codons(self):
        """
        Reads and loads all 'start_codon' rows in GENCODE GTF file.
        """
        with open(self.comprehensive_gene_annotation_gtf_file, 'r') as f:
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
        with open(self.comprehensive_gene_annotation_gtf_file, 'r') as f:
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
        with open(self.comprehensive_gene_annotation_gtf_file, 'r') as f:
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

    def read_comprehensive_gene_annotation_gtf_file(self, gtf_file: str):
        """
        Reads a GENCODE 'comprehensive gene annotation' GTF file and loads the data into self.gene_set.

        Parameters
        ----------
        gtf_file    :   GTF file.
        """
        logger.info('Started reading GTF file: %s' % gtf_file)
        self.gene_set = GeneSet()
        self.comprehensive_gene_annotation_gtf_file = gtf_file
        self.__read_gtf_file_genes()            # add genes
        self.__read_gtf_file_transcripts()      # add transcripts
        self.__read_gtf_file_exons()            # add exons
        self.__read_gtf_file_start_codons()     # update start codon start and end positions
        self.__read_gtf_file_stop_codons()      # update stop codon start and end positions
        self.__read_gtf_file_utr()              # update UTR start and end positions
        logger.info('Finished reading GTF file: %s' % gtf_file)

    def read_transcript_fasta_file(self, fasta_file: str):
        """
        Reads a GENCODE transcripts FASTA file and loads the data into self.gene_set.

        Parameters
        ----------
        fasta_file  :   GENCODE transcripts FASTA file.
        """
        logger.info('Started reading transcripts FASTA file: %s' % fasta_file)
        self.transcripts_fasta_file = fasta_file
        first = True
        curr_sequence = ''
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line[0] == '>':
                    if not first:
                        # Store sequence information
                        self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].sequence = curr_sequence
                        curr_sequence = ''
                    else:
                        first = False
                    line = line[1:]
                    line_elements = line.split('|')
                    curr_transcript_id = line_elements[0]
                    curr_gene_id = line_elements[1]
                else:
                    curr_sequence += line

    def annotate_variant_call(
            self,
            variant_call: VariantCall) -> Tuple[List[VariantAnnotation], List[VariantAnnotation]]:
        """
        Annotates an instance of the VariantCall class and returns two lists of VariantAnnotation instances.

        Parameters
        ----------
        variant_call        :   An instance of the VariantCall class.

        Returns
        -------
        pos_1_annotations   :   List of instances of the VariantAnnotation class for position 1.
        pos_2_annotations   :   List of instances of the VariantAnnotation class for position 2.
        """
        # Position 1
        pos_1_annotations = []
        for gene in self.gene_set.genes.values():
            if (gene.chromosome == variant_call.chr_1) and (gene.start <= variant_call.pos_1 <= gene.end):
                # Check if position 1 is exonic or intronic
                is_exonic = False
                for transcript in gene.transcripts:
                    for exon in transcript.exons:
                        if exon.start <= variant_call.chr_1 <= exon.end:
                            is_exonic = True
                if is_exonic:
                    region = GenomicRegionTypes.EXONIC
                else:
                    region = GenomicRegionTypes.INTRONIC
                variant_annotation = VariantAnnotation(
                    chrom=variant_call.chr_1,
                    pos=variant_call.pos_1,
                    region=region,
                    gene=deepcopy(gene)
                )
                pos_1_annotations.append(variant_annotation)
        if len(pos_1_annotations) == 0:
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chr_1,
                pos=variant_call.pos_1,
                region=GenomicRegionTypes.INTERGENIC
            )
            pos_1_annotations.append(variant_annotation)

        # Position 2
        pos_2_annotations = []
        for gene in self.gene_set.genes.values():
            if (gene.chromosome == variant_call.chr_2) and (gene.start <= variant_call.pos_2 <= gene.end):
                # Check if position 2 is exonic or intronic
                is_exonic = False
                for transcript in gene.transcripts:
                    for exon in transcript.exons:
                        if exon.start <= variant_call.chr_2 <= exon.end:
                            is_exonic = True
                if is_exonic:
                    region = GenomicRegionTypes.EXONIC
                else:
                    region = GenomicRegionTypes.INTRONIC
                variant_annotation = VariantAnnotation(
                    chrom=variant_call.chr_2,
                    pos=variant_call.pos_2,
                    region=region,
                    gene=deepcopy(gene)
                )
                pos_2_annotations.append(variant_annotation)
        if len(pos_2_annotations) == 0:
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chr_2,
                pos=variant_call.pos_2,
                region=GenomicRegionTypes.INTERGENIC
            )
            pos_2_annotations.append(variant_annotation)
        return pos_1_annotations, pos_2_annotations

    def annotate_variants(self, variants_list) -> VariantsList:
        """
        Annotates variants using self.gene_set.

        Parameters
        ----------
        variants_list   :   An instance of the VariantsList class.

        Returns
        -------
        variants_list   :   An instance of the VariantsList class with
                            each variant_call annotated.
        """
        if self.gene_set is None:
            logger.error('Please read a GENCODE comprehensive gene annotation '
                         'GTF file first before annotating variants.')
            exit(1)
        for i in range(0, len(variants_list.variant_ids)):
            for j in range(0, len(variants_list.variants[i].variant_calls)):
                pos_1_annotations, pos_2_annotations = self.annotate_variant_call(variant_call=variants_list.variants[i].variant_calls[j])
                variants_list.variants[i].variant_calls[j].pos_1_annotations = pos_1_annotations
                variants_list.variants[i].variant_calls[j].pos_2_annotations = pos_2_annotations
        return variants_list
