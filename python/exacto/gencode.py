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


# import pandas as pd
# import pysam
# from collections import defaultdict
# from dataclasses import dataclass
# from typing import List, Tuple
# from .annotation_db import AnnotationDb
# from .exon import Exon
# from .gene import Gene
# from .gene_set import GeneSet
# from .transcript import Transcript
# from ..common.nucleotide_sequence import NucleotideSequence
# from ..constants import *
# from ..logging import get_logger
# from ..vcf.variant import VariantAnnotation
# from ..vcf.variant_call import VariantCall
# from ..vcf.variants_list import VariantsList
#
#
# logger = get_logger(__name__)
#
#
# @dataclass
# class Gencode(AnnotationDb):
#     gtf_file: str
#     genome_fasta_file: str
#     version: str
#     genome: str
#     gene_set: GeneSet = GeneSet()
#     genome_fasta: pysam.FastaFile = None
#
#     def __post_init__(self):
#         # Step 1. Read FASTA file
#         self.genome_fasta = pysam.FastaFile(self.genome_fasta_file)
#
#         # Step 2. Read GTF file
#         self.__read_gtf_file_genes()  # add genes
#         self.__read_gtf_file_transcripts()  # add transcripts
#         self.__read_gtf_file_exons()  # add exons
#         self.__read_gtf_file_start_codons()  # update start codon start and end positions
#         self.__read_gtf_file_stop_codons()  # update stop codon start and end positions
#         self.__read_gtf_file_utr()  # update UTR start and end positions
#
#     @staticmethod
#     def read_gtf_file(gtf_file: str) -> pd.DataFrame:
#         """
#         Reads a GENCODE GTF file and returns a Pandas DataFrame.
#
#         Parameters
#         ----------
#         gtf_file    :   GTF file.
#
#         Returns
#         -------
#         df          :   pd.DataFrame with the following columns:
#                         'chromosome'
#                         'source'
#                         'feature'
#                         'start'
#                         'end'
#                         'score'
#                         'strand'
#                         'frame'
#                         'gene_id'
#                         'gene_type'
#                         'gene_name'
#                         'transcript_type'
#                         'transcript_name'
#                         'transcript_support_level'
#                         'exon_id'
#                         'exon_number'
#                         'protein_id'
#                         'level'
#                         'ccds_id'
#                         'hgnc_id'
#                         'havana_gene'
#                         'havana_transcript'
#                         'tag'
#                         'ont'
#         """
#         data = {
#             'chromosome': [],
#             'source': [],
#             'feature': [],
#             'start': [],
#             'end': [],
#             'score': [],
#             'strand': [],
#             'frame': [],
#             'gene_id': [],
#             'gene_type': [],
#             'gene_name': [],
#             'transcript_type': [],
#             'transcript_name': [],
#             'transcript_support_level': [],
#             'exon_id': [],
#             'exon_number': [],
#             'protein_id': [],
#             'level': [],
#             'ccds_id': [],
#             'hgnc_id': [],
#             'havana_gene': [],
#             'havana_transcript': [],
#             'tag': [],
#             'ont': []
#         }
#         with open(gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 curr_chromosome = str(elements[0])
#                 curr_source = str(elements[1])
#                 curr_feature = str(elements[2])
#                 curr_start = int(elements[3])
#                 curr_end = int(elements[4])
#                 curr_score = str(elements[5])
#                 curr_strand = str(elements[6])
#                 curr_frame = int(elements[7])
#                 curr_metadata = str(elements[8]).split(';')
#                 curr_metadata_dict = {}
#                 curr_tags = []
#                 curr_onts = []
#                 for curr_metadata_elements in curr_metadata:
#                     if curr_metadata_elements == '':
#                         continue
#                     if curr_metadata_elements[0] == ' ':
#                         curr_metadata_elements = curr_metadata_elements[1:]
#                     curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                     curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
#                     if curr_metadata_elements_[0] == 'tag':
#                         curr_tags.append(str(curr_metadata_elements_[1].replace('"', '')))
#                     if curr_metadata_elements_[0] == 'ont':
#                         curr_onts.append(str(curr_metadata_elements_[1].replace('"', '')))
#                 curr_tag_str = ';'.join(curr_tags)
#                 curr_ont_str = ';'.join(curr_onts)
#                 curr_gene_id = str(curr_metadata_dict['gene_id']) if 'gene_id' in curr_metadata_dict.keys() else ''
#                 curr_gene_type = str(curr_metadata_dict['gene_type']) if 'gene_type' in curr_metadata_dict.keys() else ''
#                 curr_gene_name = str(curr_metadata_dict['gene_name']) if 'gene_name' in curr_metadata_dict.keys() else ''
#                 curr_transcript_type = str(curr_metadata_dict['transcript_type']) if 'transcript_type' in curr_metadata_dict.keys() else ''
#                 curr_transcript_name = str(curr_metadata_dict['transcript_name']) if 'transcript_name' in curr_metadata_dict.keys() else ''
#                 curr_transcript_support_level = str(curr_metadata_dict['transcript_support_level']) if 'transcript_support_level' in curr_metadata_dict.keys() else ''
#                 curr_exon_id = str(curr_metadata_dict['exon_id']) if 'exon_id' in curr_metadata_dict.keys() else ''
#                 curr_exon_number = int(curr_metadata_dict['exon_number']) if 'exon_number' in curr_metadata_dict.keys() else -1
#                 curr_protein_id = str(curr_metadata_dict['protein_id']) if 'protein_id' in curr_metadata_dict.keys() else ''
#
#                 if 'level' in curr_metadata_dict.keys():
#                     curr_level = int(curr_metadata_dict['level'])
#                 else:
#                     curr_level = -1
#
#                 if 'ccdsid' in curr_metadata_dict.keys():
#                     curr_ccds_id = str(curr_metadata_dict['ccdsid'])
#                 else:
#                     curr_ccds_id = ''
#
#                 if 'hgnc_id' in curr_metadata_dict.keys():
#                     curr_hgnc_id = str(curr_metadata_dict['hgnc_id'])
#                 else:
#                     curr_hgnc_id = ''
#
#                 if 'havana_gene' in curr_metadata_dict.keys():
#                     curr_havana_gene = str(curr_metadata_dict['havana_gene'])
#                 else:
#                     curr_havana_gene = ''
#
#                 if 'havana_transcript' in curr_metadata_dict.keys():
#                     curr_havana_transcript = str(curr_metadata_dict)
#                 else:
#                     curr_havana_transcript = ''
#
#                 if curr_score == '.':
#                     curr_score = -1.0
#                 else:
#                     curr_score = float(curr_score)
#
#                 data['chromosome'].append(curr_chromosome)
#                 data['source'].append(curr_source)
#                 data['feature'].append(curr_feature)
#                 data['start'].append(curr_start)
#                 data['end'].append(curr_end)
#                 data['score'].append(curr_score)
#                 data['strand'].append(curr_strand)
#                 data['frame'].append(curr_frame)
#                 data['gene_id'].append(curr_gene_id)
#                 data['gene_type'].append(curr_gene_type)
#                 data['transcript_type'].append(curr_transcript_type)
#                 data['transcript_name'].append(curr_transcript_name)
#                 data['transcript_support_level'].append(curr_transcript_support_level)
#                 data['exon_id'].append(curr_exon_id)
#                 data['exon_number'].append(curr_exon_number)
#                 data['protein_id'].append(curr_protein_id)
#                 data['level'].append(curr_level)
#                 data['ccds_id'].append(curr_ccds_id)
#                 data['hgnc_id'].append(curr_hgnc_id)
#                 data['havana_gene'].append(curr_havana_gene)
#                 data['havana_transcript'].append(curr_havana_transcript)
#                 data['tag'].append(curr_tag_str)
#                 data['ont'].append(curr_ont_str)
#         return pd.DataFrame(data)
#
#     @staticmethod
#     def get_stable_ensembl_id(id: str):
#         """
#         Returns the stable Ensembl ID and version.
#
#         Parameters
#         ----------
#         id          :   Ensembl ID (e.g. 'ENSG00001.1').
#
#         Returns
#         -------
#         stable_id   :   Stable Ensembl ID (without version number).
#         version     :   Version.
#         """
#         if (('ENSG' in id) or ('ENST' in id or 'ENSE' in id)) and '.' in id:
#             stable_id = id.split('.')[0]
#             version = id.split('.')[1]
#             return stable_id, version
#         else:
#             return id, None
#
#     def __read_gtf_file_genes(self):
#         """
#         Reads and loads all 'gene' rows in GENCODE GTF file.
#         """
#         with open(self.gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 if elements[2] == 'gene':
#                     curr_gene_chrom = str(elements[0])
#                     curr_gene_source = str(elements[1])
#                     curr_gene_start = int(elements[3])
#                     curr_gene_end = int(elements[4])
#                     curr_gene_strand = str(elements[6])
#                     curr_metadata = str(elements[8]).split(';')
#                     curr_metadata_dict = {}
#                     for curr_metadata_elements in curr_metadata:
#                         if curr_metadata_elements == '':
#                             continue
#                         if curr_metadata_elements[0] == ' ':
#                             curr_metadata_elements = curr_metadata_elements[1:]
#                         curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                         curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
#                     curr_gene_id = str(curr_metadata_dict['gene_id'])
#                     curr_gene_stable_id, curr_gene_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['gene_id']))
#                     curr_gene_name = str(curr_metadata_dict['gene_name'])
#                     curr_gene_type = str(curr_metadata_dict['gene_type'])
#                     curr_gene_level = int(curr_metadata_dict['level'])
#                     gene = Gene(
#                         id=curr_gene_id,
#                         stable_id=curr_gene_stable_id,
#                         source=curr_gene_source,
#                         source_version=self.version,
#                         name=curr_gene_name,
#                         chromosome=curr_gene_chrom,
#                         start=curr_gene_start,
#                         end=curr_gene_end,
#                         strand=curr_gene_strand,
#                         type=curr_gene_type,
#                         level=curr_gene_level,
#                         version=curr_gene_version,
#                         genome=self.genome
#                     )
#                     self.gene_set.add_gene(gene=gene)
#         logger.info('Loaded %i genes in total.' % len(self.gene_set.gene_ids))
#
#     def __read_gtf_file_transcripts(self):
#         """
#         Reads and loads all 'transcript' rows in GENCODE GTF file.
#         """
#         with open(self.gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 if elements[2] == 'transcript':
#                     curr_transcript_chrom = str(elements[0])
#                     curr_transcript_source = str(elements[1])
#                     curr_transcript_start = int(elements[3])
#                     curr_transcript_end = int(elements[4])
#                     curr_transcript_strand = str(elements[6])
#                     curr_metadata = str(elements[8]).split(';')
#                     curr_metadata_dict = defaultdict(list)
#                     for curr_metadata_elements in curr_metadata:
#                         if curr_metadata_elements == '':
#                             continue
#                         if curr_metadata_elements[0] == ' ':
#                             curr_metadata_elements = curr_metadata_elements[1:]
#                         curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                         curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
#                     curr_gene_id = str(curr_metadata_dict['gene_id'][0])
#                     curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
#                     curr_transcript_stable_id, curr_transcript_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['transcript_id'][0]))
#                     curr_transcript_type = str(curr_metadata_dict['transcript_type'][0])
#                     curr_transcript_name = str(curr_metadata_dict['transcript_name'][0])
#                     curr_transcript_tags = [str(tag).replace('"', '') for tag in curr_metadata_dict['tag']]
#                     try:
#                         curr_transcript_level = int(curr_metadata_dict['level'][0])
#                     except:
#                         curr_transcript_level = None
#                     try:
#                         curr_transcript_support_level = int(curr_metadata_dict['transcript_support_level'][0])
#                     except:
#                         curr_transcript_support_level = None
#                     transcript = Transcript(
#                         id=curr_transcript_id,
#                         stable_id=curr_transcript_stable_id,
#                         source=curr_transcript_source,
#                         source_version=self.version,
#                         chromosome=curr_transcript_chrom,
#                         start=curr_transcript_start,
#                         end=curr_transcript_end,
#                         type=curr_transcript_type,
#                         strand=curr_transcript_strand,
#                         version=curr_transcript_version,
#                         name=curr_transcript_name,
#                         level=curr_transcript_level,
#                         support_level=curr_transcript_support_level,
#                         tags=curr_transcript_tags,
#                         genome=self.genome
#                     )
#                     self.gene_set.add_transcript(gene_id=curr_gene_id, transcript=transcript)
#
#         transcripts_count = 0
#         for gene in self.gene_set.genes:
#             transcripts_count += gene.transcripts_count
#         logger.info('Loaded %i transcripts in total.' % transcripts_count)
#
#     def __read_gtf_file_exons(self):
#         """
#         Reads and loads all 'exon' rows in GENCODE GTF file.
#         """
#         with open(self.gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 if elements[2] == 'exon':
#                     curr_exon_chrom = str(elements[0])
#                     curr_exon_source = str(elements[1])
#                     curr_exon_start = int(elements[3])
#                     curr_exon_end = int(elements[4])
#                     curr_exon_strand = str(elements[6])
#                     curr_metadata = str(elements[8]).split(';')
#                     curr_metadata_dict = defaultdict(list)
#                     curr_exon_tags = []
#                     for curr_metadata_elements in curr_metadata:
#                         if curr_metadata_elements == '':
#                             continue
#                         if curr_metadata_elements[0] == ' ':
#                             curr_metadata_elements = curr_metadata_elements[1:]
#                         curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                         curr_metadata_dict[curr_metadata_elements_[0]].append(curr_metadata_elements_[1].replace('"', ''))
#                     curr_gene_id = str(curr_metadata_dict['gene_id'][0])
#                     curr_transcript_id = str(curr_metadata_dict['transcript_id'][0])
#                     curr_exon_id = str(curr_metadata_dict['exon_id'][0])
#                     curr_exon_stable_id, curr_exon_version = Gencode.get_stable_ensembl_id(id=str(curr_metadata_dict['exon_id'][0]))
#                     curr_exon_number = int(curr_metadata_dict['exon_number'][0])
#                     curr_exon_tags = [str(tag).replace('"', '') for tag in curr_metadata_dict['tag']]
#                     if curr_exon_strand == Strands.POSITIVE:
#                         sequence = self.genome_fasta.fetch(curr_exon_chrom, curr_exon_start - 1, curr_exon_end)
#                     else:
#                         sequence = self.genome_fasta.fetch(curr_exon_chrom, curr_exon_start - 1, curr_exon_end)
#                         sequence = NucleotideSequence.reverse_complement(sequence=sequence)
#                     exon = Exon(
#                         id=curr_exon_id,
#                         stable_id=curr_exon_stable_id,
#                         source=curr_exon_source,
#                         source_version=self.version,
#                         chromosome=curr_exon_chrom,
#                         start=curr_exon_start,
#                         end=curr_exon_end,
#                         number=curr_exon_number,
#                         strand=curr_exon_strand,
#                         sequence=NucleotideSequence(sequence=sequence),
#                         version=curr_exon_version,
#                         tags=curr_exon_tags
#                     )
#                     self.gene_set.add_exon(
#                         gene_id=curr_gene_id,
#                         transcript_id=curr_transcript_id,
#                         exon=exon
#                     )
#         exons_count = 0
#         for gene in self.gene_set.genes:
#             for transcript in gene.transcripts:
#                 exons_count += transcript.exons_count
#         logger.info('Loaded %i exons in total.' % exons_count)
#
#     def __read_gtf_file_start_codons(self):
#         """
#         Reads and loads all 'start_codon' rows in GENCODE GTF file.
#         """
#         with open(self.gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 if elements[2] == 'start_codon':
#                     curr_start_codon_start = int(elements[3])
#                     curr_start_codon_end = int(elements[4])
#                     curr_metadata = str(elements[8]).split(';')
#                     curr_metadata_dict = {}
#                     for curr_metadata_elements in curr_metadata:
#                         if curr_metadata_elements == '':
#                             continue
#                         if curr_metadata_elements[0] == ' ':
#                             curr_metadata_elements = curr_metadata_elements[1:]
#                         curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                         curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
#                     curr_gene_id = str(curr_metadata_dict['gene_id'])
#                     curr_transcript_id = str(curr_metadata_dict['transcript_id'])
#                     self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].start_codon_start = curr_start_codon_start
#                     self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].start_codon_end = curr_start_codon_end
#
#     def __read_gtf_file_stop_codons(self):
#         """
#         Reads and loads all 'stop_codon' rows in GENCODE GTF file.
#         """
#         with open(self.gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 if elements[2] == 'stop_codon':
#                     curr_stop_codon_start = int(elements[3])
#                     curr_stop_codon_end = int(elements[4])
#                     curr_metadata = str(elements[8]).split(';')
#                     curr_metadata_dict = {}
#                     for curr_metadata_elements in curr_metadata:
#                         if curr_metadata_elements == '':
#                             continue
#                         if curr_metadata_elements[0] == ' ':
#                             curr_metadata_elements = curr_metadata_elements[1:]
#                         curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                         curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
#
#                     curr_gene_id = str(curr_metadata_dict['gene_id'])
#                     curr_transcript_id = str(curr_metadata_dict['transcript_id'])
#                     self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].stop_codon_start = curr_stop_codon_start
#                     self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].stop_codon_end = curr_stop_codon_end
#
#     def __read_gtf_file_utr(self):
#         """
#         Reads and loads all 'UTR' rows in GENCODE GTF file.
#         """
#         with open(self.gtf_file, 'r') as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line[0:2] == '##':
#                     continue
#                 elements = line.split('\t')
#                 if elements[2] == 'UTR':
#                     curr_utr_start = int(elements[3])
#                     curr_utr_end = int(elements[4])
#                     curr_metadata = str(elements[8]).split(';')
#                     curr_metadata_dict = {}
#                     for curr_metadata_elements in curr_metadata:
#                         if curr_metadata_elements == '':
#                             continue
#                         if curr_metadata_elements[0] == ' ':
#                             curr_metadata_elements = curr_metadata_elements[1:]
#                         curr_metadata_elements_ = curr_metadata_elements.split(' ')
#                         curr_metadata_dict[curr_metadata_elements_[0]] = curr_metadata_elements_[1].replace('"', '')
#
#                     curr_gene_id = str(curr_metadata_dict['gene_id'])
#                     curr_transcript_id = str(curr_metadata_dict['transcript_id'])
#                     self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].utr_start = curr_utr_start
#                     self.gene_set.genes[curr_gene_id].transcripts[curr_transcript_id].utr_end = curr_utr_end
#
#     def annotate_position(
#             self,
#             chromosome: str,
#             position: int
#     ) -> List[VariantAnnotation]:
#         """
#         Annotates a position and returns a list of VariantAnnotation objects.
#
#         Parameters
#         ----------
#         chromosome              :   Chromosome.
#         position                :   Position.
#
#         Returns
#         -------
#         variant_annotations     :   List of VariantAnnotation objects.
#         """
#         variant_annotations = []
#         genes = self.gene_set.find_genes(chromosome=chromosome, position=position)
#         if len(genes) == 0:
#             variant_annotation = VariantAnnotation(
#                 region=GenomicRegionTypes.INTERGENIC,
#                 source=AnnotationSources.ENSEMBL,
#                 source_version=self.version
#             )
#             variant_annotations.append(variant_annotation)
#         else:
#             for gene in genes:
#                 region = GenomicRegionTypes.INTRONIC
#                 for transcript in gene.transcripts:
#                     for exon in transcript.exons:
#                         if exon.start <= position <= exon.end:
#                             region = GenomicRegionTypes.EXONIC
#                 variant_annotation = VariantAnnotation(
#                     region=region,
#                     source=AnnotationSources.GENCODE,
#                     source_version=self.version,
#                     gene_id=gene.id,
#                     gene_stable_id=gene.stable_id,
#                     gene_version=gene.version,
#                     gene_name=gene.name,
#                     gene_type=gene.type,
#                     gene_strand=gene.strand,
#                     genome=gene.genome
#                 )
#                 variant_annotations.append(variant_annotation)
#         return variant_annotations
#
#     def annotate_variant_call(
#             self,
#             variant_call: VariantCall
#     ) -> Tuple[List[VariantAnnotation], List[VariantAnnotation]]:
#         """
#         Annotates a VariantCall object and returns two lists of VariantAnnotation objects.
#
#         Parameters
#         ----------
#         variant_call        :   VariantCall object.
#
#         Returns
#         -------
#         pos_1_annotations   :   List of VariantAnnotation objects for position 1.
#         pos_2_annotations   :   List of VariantAnnotation objects for position 2.
#         """
#         pos_1_annotations = self.annotate_position(
#             chromosome=variant_call.chromosome_1,
#             position=variant_call.position_1
#         )
#         pos_2_annotations = self.annotate_position(
#             chromosome=variant_call.chromosome_2,
#             position=variant_call.position_2
#         )
#         return pos_1_annotations, pos_2_annotations
#
#     def annotate_variants_list(self, variants_list) -> VariantsList:
#         """
#         Annotates a VariantsList object.
#
#         Parameters
#         ----------
#         variants_list   :   VariantsList object.
#
#         Returns
#         -------
#         variants_list   :   VariantsList object with each VariantCall object annotated.
#         """
#         if self.gene_set is None:
#             raise Exception('Please read a GENCODE comprehensive gene annotation '
#                             'GTF file before annotating a VariantsList object.')
#         for key in variants_list.variants.keys():
#             for i in range(0, len(variants_list.variants[key])):
#                 for j in range(0, len(variants_list.variants[key][i].variant_calls)):
#                     position_1_annotations, position_2_annotations = self.annotate_variant_call(
#                         variants_list.variants[key][i].variant_calls[j]
#                     )
#                     for annotation in position_1_annotations:
#                         variants_list.variants[key][i].variant_calls[j].position_1_annotations.append(annotation)
#                     for annotation in position_2_annotations:
#                         variants_list.variants[key][i].variant_calls[j].position_2_annotations.append(annotation)
#         return variants_list
