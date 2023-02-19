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
The purpose of this python3 script is to implement functions related to
annotating variants using GENCODE.
"""


import pandas as pd
from typing import Tuple, List
from ...logging import get_logger


logger = get_logger(__name__)


class GencodeExon:

    def __init__(
            self,
            exon_id: str,
            exon_number: int,
            exon_chrom: str,
            exon_start: int,
            exon_end: int
        ):
        self.__exon_id = exon_id
        self.__exon_number = exon_number
        self.__exon_chrom = exon_chrom
        self.__exon_start = exon_start
        self.__exon_end = exon_end

    def get_exon_id(self) -> str:
        return self.__exon_id

    def get_exon_number(self) -> int:
        return self.__exon_number

    def get_exon_chrom(self) -> str:
        return self.__exon_chrom

    def get_exon_start(self) -> int:
        return self.__exon_start

    def get_exon_end(self) -> int:
        return self.__exon_end


class GencodeTranscript:

    def __init__(
            self,
            gene_id: str,
            transcript_id: str,
            transcript_name: str,
            transcript_type: str,
            transcript_chrom: str,
            transcript_start: int,
            transcript_end: int,
            transcript_strand:str,
            level: int
        ):
        self.__gene_id = gene_id
        self.__transcript_id = transcript_id
        self.__transcript_name = transcript_name
        self.__transcript_type = transcript_type
        self.__transcript_chrom = transcript_chrom
        self.__transcript_start = transcript_start
        self.__transcript_end = transcript_end
        self.__transcript_strand = transcript_strand
        self.__level = level
        self.__start_codon_start = -1
        self.__start_codon_end = -1
        self.__stop_codon_start = -1
        self.__stop_codon_end = -1
        self.__utr_start = -1
        self.__utr_end = -1
        self.__exons = []

    def set_start_codon(
            self,
            start: int,
            end: int
        ) -> None:
        self.__start_codon_start = start
        self.__start_codon_end = end

    def set_stop_codon(
            self,
            start: int,
            end: int
        ) -> None:
        self.__stop_codon_start = start
        self.__stop_codon_end = end

    def set_utr(
            self,
            start: int,
            end: int
        ) -> None:
        self.__utr_start = start
        self.__utr_end = end

    def add_exon(
            self,
            exon: GencodeExon
        ) -> None:
        self.__exons.append(exon)

    def get_gene_id(self) -> str:
        return self.__gene_id

    def get_transcript_id(self) -> str:
        return self.__transcript_id

    def get_transcript_name(self) -> str:
        return self.__transcript_name

    def get_transcript_type(self) -> str:
        return self.__transcript_type

    def get_transcript_chrom(self) -> str:
        return self.__transcript_chrom

    def get_transcript_start(self) -> int:
        return self.__transcript_start

    def get_transcript_end(self) -> int:
        return self.__transcript_end

    def get_transcript_strand(self) -> str:
        return self.__transcript_strand

    def get_level(self) -> int:
        return self.__level

    def get_start_codon(self) -> int:
        return self.__start_codon_start, self.__start_codon_end

    def get_stop_codon(self) -> int:
        return self.__stop_codon_start, self.__stop_codon_end

    def get_utr(self) -> Tuple[int, int]:
        return self.__utr_start, self.__utr_end

    def get_exons(self) -> List[GencodeExon]:
        return self.__exons

    def get_exon_ids(self) -> List[str]:
        exon_ids = [i.get_exon_id() for i in self.__exons]
        return exon_ids


class GencodeGene:

    def __init__(
            self,
            gene_id: str,
            gene_name: str,
            gene_type: str,
            gene_chrom: str,
            gene_start: int,
            gene_end: int,
            gene_strand: str,
            level: int
        ):
        self.__gene_id = gene_id
        self.__gene_name = gene_name
        self.__gene_type = gene_type
        self.__gene_chrom = gene_chrom
        self.__gene_start = gene_start
        self.__gene_end = gene_end
        self.__gene_strand = gene_strand
        self.__level = level
        self.__transcripts = []

    def set_start_codon(
            self,
            transcript_id: str,
            start_codon_start: int,
            start_codon_end: int
        ) -> None:
        for i in range(0, len(self.__transcripts)):
            if self.__transcripts[i].get_transcript_id() == transcript_id:
                self.__transcripts[i].set_start_codon(start=start_codon_start,
                                                      end=start_codon_end)
                return

    def set_stop_codon(
            self,
            transcript_id: str,
            stop_codon_start: int,
            stop_codon_end: int
        ) -> None:
        for i in range(0, len(self.__transcripts)):
            if self.__transcripts[i].get_transcript_id() == transcript_id:
                self.__transcripts[i].set_stop_codon(start=stop_codon_start,
                                                     end=stop_codon_end)
                return

    def set_utr(
            self,
            transcript_id: str,
            utr_start: int,
            utr_end: int
        ) -> None:
        for i in range(0, len(self.__transcripts)):
            if self.__transcripts[i].get_transcript_id() == transcript_id:
                self.__transcripts[i].set_utr(start=utr_start,
                                              end=utr_end)
                return

    def get_gene_id(self) -> str:
        return self.__gene_id

    def get_gene_name(self) -> str:
        return self.__gene_name

    def get_gene_type(self) -> str:
        return self.__gene_type

    def get_gene_chrom(self) -> str:
        return self.__gene_chrom

    def get_gene_start(self) -> int:
        return self.__gene_start

    def get_gene_end(self) -> int:
        return self.__gene_end

    def get_gene_strand(self) -> str:
        return self.__gene_strand

    def get_level(self) -> int:
        return self.__level

    def add_transcript(
            self,
            transcript: GencodeTranscript
        ) -> None:
        self.__transcripts.append(transcript)

    def get_transcripts(self) -> List[GencodeTranscript]:
        return self.__transcripts

    def get_transcript_ids(self) -> List[str]:
        transcript_ids = [i.get_transcript_id() for i in self.__transcripts]
        return transcript_ids

    def add_exon(
            self,
            transcript_id: str,
            exon: GencodeExon
        ) -> None:
        for i in range(0, len(self.__transcripts)):
            if self.__transcripts[i].get_transcript_id() == transcript_id:
                if exon.get_exon_id() not in self.__transcripts[i].get_exon_ids():
                    self.__transcripts[i].add_exon(exon=exon)
                    return


def read_gencode_gtf_file(
        gencode_gtf_file: str
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads a GENCODE GTF file and returns a DataFrame.

    Parameters
    ----------
    gencode_gtf_file    :   GTF file path.

    Returns
    -------
    df_genes            :   DataFrame with the following columns:
                            'gene_id'
                            'gene_name'
                            'gene_type'
                            'gene_chrom'
                            'gene_start'
                            'gene_end'
                            'gene_strand'
                            'level'
                            'transcripts_count'
    df_transcripts      :   DataFrame with the following columns:
                            'gene_id'
                            'transcript_id'
                            'transcript_name'
                            'transcript_type'
                            'transcript_chrom'
                            'transcript_start'
                            'transcript_end'
                            'transcript_strand'
                            'level'
                            'start_codon_start'
                            'start_codon_end'
                            'stop_codon_start'
                            'stop_codon_end'
                            'utr_start'
                            'utr_end'
                            'exons_count'
    df_exons            :   DataFrame with the following columns:
                            'gene_id'
                            'transcript_id'
                            'exon_id'
                            'exon_number'
                            'exon_chrom'
                            'exon_start'
                            'exon_end'
    """
    logger.info('Started reading GENCODE GTF file.')

    # Step 1. Add all genes
    gencode_genes_dict = {} # key = gene_id, value = an instance of GENCODE_Gene
    with open(gencode_gtf_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0:2] == '##':
                continue
            elements = line.split('\t')
            if elements[2] == 'gene':
                curr_gene_chrom = str(elements[0])
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
                if curr_gene_id not in gencode_genes_dict.keys():
                    gencode_gene = GencodeGene(
                        gene_id=str(curr_metadata_dict['gene_id']),
                        gene_name=str(curr_metadata_dict['gene_name']),
                        gene_type=str(curr_metadata_dict['gene_type']),
                        gene_chrom=curr_gene_chrom,
                        gene_start=curr_gene_start,
                        gene_end=curr_gene_end,
                        gene_strand=curr_gene_strand,
                        level=int(curr_metadata_dict['level'])
                    )
                    gencode_genes_dict[curr_gene_id] = gencode_gene
    logger.info('%i genes in total.' % len(gencode_genes_dict.keys()))

    # Step 2. Add all transcripts
    with open(gencode_gtf_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0:2] == '##':
                continue
            elements = line.split('\t')
            if elements[2] == 'transcript':
                curr_transcript_chrom = str(elements[0])
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
                curr_transcript_id = str(curr_metadata_dict['transcript_id'])
                try:
                    transcript_level = int(curr_metadata_dict['transcript_support_level'])
                except:
                    transcript_level = -1
                if curr_transcript_id not in gencode_genes_dict[curr_gene_id].get_transcript_ids():
                    gencode_transcript = GencodeTranscript(
                        gene_id=str(curr_gene_id),
                        transcript_id=str(curr_metadata_dict['transcript_id']),
                        transcript_name=str(curr_metadata_dict['transcript_name']),
                        transcript_type=str(curr_metadata_dict['transcript_type']),
                        transcript_chrom=curr_transcript_chrom,
                        transcript_start=curr_transcript_start,
                        transcript_end=curr_transcript_end,
                        transcript_strand=curr_transcript_strand,
                        level=transcript_level
                    )
                    gencode_genes_dict[curr_gene_id].add_transcript(transcript=gencode_transcript)

    transcripts_count = 0
    for key,val in gencode_genes_dict.items():
        transcripts_count += len(val.get_transcript_ids())
    logger.info('%i transcripts in total.' % transcripts_count)

    # Step 3. Add all exons
    with open(gencode_gtf_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0:2] == '##':
                continue
            elements = line.split('\t')
            if elements[2] == 'exon':
                curr_exon_chrom = str(elements[0])
                curr_exon_start = int(elements[3])
                curr_exon_end = int(elements[4])
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
                gencode_exon = GencodeExon(
                    exon_id=str(curr_metadata_dict['exon_id']),
                    exon_number=int(curr_metadata_dict['exon_number']),
                    exon_chrom=curr_exon_chrom,
                    exon_start=curr_exon_start,
                    exon_end=curr_exon_end
                )
                gencode_genes_dict[curr_gene_id].add_exon(transcript_id=curr_transcript_id, exon=gencode_exon)

    exons_count = 0
    for key, gene in gencode_genes_dict.items():
        for transcript in gene.get_transcripts():
            exons_count += len(transcript.get_exon_ids())
    logger.info('%i exons in total.' % exons_count)

    # Step 4. Add start codons
    with open(gencode_gtf_file, 'r') as f:
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
                gencode_genes_dict[curr_gene_id].set_start_codon(
                    transcript_id=curr_transcript_id,
                    start_codon_start=curr_start_codon_start,
                    start_codon_end=curr_start_codon_end
                )

    # Step 4. Add stop codons
    with open(gencode_gtf_file, 'r') as f:
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
                gencode_genes_dict[curr_gene_id].set_stop_codon(
                    transcript_id=curr_transcript_id,
                    stop_codon_start=curr_stop_codon_start,
                    stop_codon_end=curr_stop_codon_end
                )

    # Step 4. Add UTR
    with open(gencode_gtf_file, 'r') as f:
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
                gencode_genes_dict[curr_gene_id].set_utr(
                    transcript_id=curr_transcript_id,
                    utr_start=curr_utr_start,
                    utr_end=curr_utr_end
                )

    # Step 5. Store genes, transcripts, and exons
    genes_data = {
        'gene_id': [],
        'gene_name': [],
        'gene_type': [],
        'gene_chrom': [],
        'gene_start': [],
        'gene_end': [],
        'gene_strand': [],
        'level': [],
        'transcripts_count': []
    }
    transcripts_data = {
        'gene_id': [],
        'transcript_id': [],
        'transcript_name': [],
        'transcript_type': [],
        'transcript_chrom': [],
        'transcript_start': [],
        'transcript_end': [],
        'transcript_strand': [],
        'level': [],
        'start_codon_start': [],
        'start_codon_end': [],
        'stop_codon_start': [],
        'stop_codon_end': [],
        'utr_start': [],
        'utr_end': [],
        'exons_count': []
    }
    exons_data = {
        'gene_id': [],
        'transcript_id': [],
        'exon_id': [],
        'exon_number': [],
        'exon_chrom': [],
        'exon_start': [],
        'exon_end': []
    }
    for curr_gene_id, gene in gencode_genes_dict.items():
        genes_data['gene_id'].append(gene.get_gene_id())
        genes_data['gene_name'].append(gene.get_gene_name())
        genes_data['gene_type'].append(gene.get_gene_type())
        genes_data['gene_chrom'].append(gene.get_gene_chrom())
        genes_data['gene_start'].append(gene.get_gene_start())
        genes_data['gene_end'].append(gene.get_gene_end())
        genes_data['gene_strand'].append(gene.get_gene_strand())
        genes_data['level'].append(gene.get_level())
        genes_data['transcripts_count'].append(len(gene.get_transcript_ids()))
        for curr_transcript in gene.get_transcripts():
            transcripts_data['gene_id'].append(curr_gene_id)
            transcripts_data['transcript_id'].append(curr_transcript.get_transcript_id())
            transcripts_data['transcript_name'].append(curr_transcript.get_transcript_name())
            transcripts_data['transcript_type'].append(curr_transcript.get_transcript_type())
            transcripts_data['transcript_chrom'].append(curr_transcript.get_transcript_chrom())
            transcripts_data['transcript_start'].append(curr_transcript.get_transcript_start())
            transcripts_data['transcript_end'].append(curr_transcript.get_transcript_end())
            transcripts_data['transcript_strand'].append(curr_transcript.get_transcript_strand())
            transcripts_data['level'].append(curr_transcript.get_level())
            transcripts_data['start_codon_start'].append(curr_transcript.get_start_codon()[0])
            transcripts_data['start_codon_end'].append(curr_transcript.get_start_codon()[1])
            transcripts_data['stop_codon_start'].append(curr_transcript.get_stop_codon()[0])
            transcripts_data['stop_codon_end'].append(curr_transcript.get_stop_codon()[1])
            transcripts_data['utr_start'].append(curr_transcript.get_utr()[0])
            transcripts_data['utr_end'].append(curr_transcript.get_utr()[1])
            transcripts_data['exons_count'].append(len(curr_transcript.get_exon_ids()))
            for curr_exon in curr_transcript.get_exons():
                exons_data['gene_id'].append(curr_gene_id)
                exons_data['transcript_id'].append(curr_transcript.get_transcript_id())
                exons_data['exon_id'].append(curr_exon.get_exon_id())
                exons_data['exon_number'].append(curr_exon.get_exon_number())
                exons_data['exon_chrom'].append(curr_exon.get_exon_chrom())
                exons_data['exon_start'].append(curr_exon.get_exon_start())
                exons_data['exon_end'].append(curr_exon.get_exon_end())
    df_genes = pd.DataFrame(genes_data)
    df_transcripts = pd.DataFrame(transcripts_data)
    df_exons = pd.DataFrame(exons_data)
    logger.info('Finished reading GENCODE GTF file.')
    return df_genes, df_transcripts, df_exons


def read_gencode_refseq_file(
        gencode_refseq_metadata_file: str
    ) -> pd.DataFrame:
    """
    Reads a GENCODE RefSeq metadata file and returns a DataFrame

    Parameters
    ----------
    gencode_refseq_metadata_file    :   GENCODE RefSeq metadata file.

    Returns
    -------
    DataFrame with the following columns:
    'ensembl_transcript_id'
    'refseq_transcript_id'
    'refseq_protein_id'
    """
    logger.info('Started reading GENCODE RefSeq file.')
    df_refseq = pd.read_csv(gencode_refseq_metadata_file,
                            sep='\t',
                            header=None)
    df_refseq.columns = ['ensembl_transcript_id',
                         'refseq_transcript_id',
                         'refseq_protein_id']
    logger.info('Finished reading GENCODE RefSeq file.')
    return df_refseq


def read_gencode_transcripts_fasta_file(
        gencode_fasta_file: str
    ) -> pd.DataFrame:
    """

    Parameters
    ----------
    gencode_fasta_file  :   GENCODE FASTA file.

    Returns
    -------
    DataFrame with the following columns:
    'ensembl_transcript_id'
    'ensembl_gene_id'
    'havana_gene_id'
    'havana_transcript_id'
    'gene_name_versioned'
    'gene_name'
    'transcript_length'
    'utr5_start'
    'utr5_end'
    'cds_start'
    'cds_end'
    'utr3_start'
    'utr3_end'
    'sequence'
    """
    logger.info('Started reading GENCODE transcripts FASTA file.')

    data = {
        'ensembl_transcript_id': [],
        'ensembl_gene_id': [],
        'havana_gene_id': [],
        'havana_transcript_id': [],
        'gene_name_versioned': [],
        'gene_name': [],
        'transcript_length': [],
        'utr5_start': [],
        'utr5_end': [],
        'cds_start': [],
        'cds_end': [],
        'utr3_start': [],
        'utr3_end': [],
        'sequence': []
    }
    first = True
    curr_sequence = ''

    with open(gencode_fasta_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0] == '>':
                # Append previously stored information
                if first:
                    first = False
                else:
                    data['ensembl_transcript_id'].append(curr_ensembl_transcript_id)
                    data['ensembl_gene_id'].append(curr_ensembl_gene_id)
                    data['havana_gene_id'].append(curr_havana_gene_id)
                    data['havana_transcript_id'].append(curr_havana_transcript_id)
                    data['gene_name_versioned'].append(curr_gene_name_versioned)
                    data['gene_name'].append(curr_gene_name)
                    data['transcript_length'].append(curr_transcript_length)
                    data['utr5_start'].append(curr_utr5_start)
                    data['utr5_end'].append(curr_utr5_end)
                    data['cds_start'].append(curr_cds_start)
                    data['cds_end'].append(curr_cds_end)
                    data['utr3_start'].append(curr_utr3_start)
                    data['utr3_end'].append(curr_utr3_end)
                    data['sequence'].append(curr_sequence)
                    curr_sequence = ''

                line = line[1:]
                line_elements = line.split('|')
                curr_ensembl_transcript_id = line_elements[0]
                curr_ensembl_gene_id = line_elements[1]
                curr_havana_gene_id = line_elements[2]
                curr_havana_transcript_id = line_elements[3]
                curr_gene_name_versioned = line_elements[4]
                curr_gene_name = line_elements[5]
                curr_transcript_length = int(line_elements[6])
                curr_utr5_start = -1
                curr_utr5_end = -1
                curr_cds_start = -1
                curr_cds_end = -1
                curr_utr3_start = -1
                curr_utr3_end = -1
                for element in line_elements:
                    if 'UTR5:' in element:
                        curr_utr5_start = int(element.split(':')[1].split('-')[0])
                        curr_utr5_end = int(element.split(':')[1].split('-')[1])
                    if 'CDS:' in element:
                        curr_cds_start = int(element.split(':')[1].split('-')[0])
                        curr_cds_end = int(element.split(':')[1].split('-')[1])
                    if 'UTR3:' in element:
                        curr_utr3_start = int(element.split(':')[1].split('-')[0])
                        curr_utr3_end = int(element.split(':')[1].split('-')[1])
            else:
                curr_sequence += line

        # Store last transcript
        data['ensembl_transcript_id'].append(curr_ensembl_transcript_id)
        data['ensembl_gene_id'].append(curr_ensembl_gene_id)
        data['havana_gene_id'].append(curr_havana_gene_id)
        data['havana_transcript_id'].append(curr_havana_transcript_id)
        data['gene_name_versioned'].append(curr_gene_name_versioned)
        data['gene_name'].append(curr_gene_name)
        data['transcript_length'].append(curr_transcript_length)
        data['utr5_start'].append(curr_utr5_start)
        data['utr5_end'].append(curr_utr5_end)
        data['cds_start'].append(curr_cds_start)
        data['cds_end'].append(curr_cds_end)
        data['utr3_start'].append(curr_utr3_start)
        data['utr3_end'].append(curr_utr3_end)
        data['sequence'].append(curr_sequence)

    df = pd.DataFrame(data)
    logger.info('Finished reading GENCODE transcripts FASTA file.')
    return df


def subset_gencode_dataframes(
        df_target_regions: pd.DataFrame,
        df_genes: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Subsets GENCODE DataFrames by target regions.

    Parameters
    ----------
    df_target_regions   :   DataFrame with the following columns:
    df_genes            :   DataFrame of genes.
    df_transcripts      :   DataFrame of transcripts.
    df_exons            :   DataFrame of exons.

    Returns
    -------
    df_genes            :   DataFrame of genes.
    df_transcripts      :   DataFrame of transcripts.
    df_exons            :   DataFrame of exons.
    """
    # Filter genes
    df_genes_filtered = pd.DataFrame()
    for _, row in df_target_regions.iterrows():
        df_genes_matched = df_genes.loc[
            (df_genes['gene_chrom'] == row['chrom']) &
            (df_genes['gene_start'] >= row['start']) &
            (df_genes['gene_end'] <= row['end']),:
        ]
        if len(df_genes_matched) > 0:
            df_genes_filtered = pd.concat([df_genes_filtered, df_genes_matched])
    df_genes = df_genes_filtered.drop_duplicates()

    # Check if any gene falls within the specified target regions
    if len(df_genes) == 0:
        logger.error("No gene falls within specified target regions.")
        exit()

    # Filter transcripts
    df_transcripts = df_transcripts.loc[
        df_transcripts['gene_id'].isin(df_genes['gene_id'].unique()),:
    ]

    # Filter exons
    df_exons = df_exons.loc[
        df_exons['transcript_id'].isin(df_transcripts['transcript_id'].unique()),:
    ]
    return df_genes, df_transcripts, df_exons


def annotate_small_variants_using_gencode(
        df_small_variants: pd.DataFrame,
        df_gencode_genes: pd.DataFrame,
        df_gencode_exons: pd.DataFrame
    ) -> pd.DataFrame:
    """
    Annotates small variants using GENCODE.

    Parameters
    ----------
    df_small_variants       :   DataFrame of small variants.
                                Expected columns:
                                'variant_id'
                                'chrom'
                                'pos'
    df_gencode_genes        :   DataFrame of GENCODE genes.
                                Expected columns:
                                'gene_id'
                                'gene_name'
                                'gene_type'
                                'gene_chrom'
                                'gene_start'
                                'gene_end'
                                'gene_strand'
                                'level'
                                'transcripts_count'
    df_gencode_exons        :   DataFrame of GENCODE exons.
                                Expected columns:
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'exon_number'
                                'exon_chrom'
                                'exon_start'
                                'exon_end'

    Returns
    -------
    DataFrame with the following columns appended:
    'ensembl_pos_region'
    'ensembl_pos_gene_id'
    'ensembl_pos_gene_name'
    'ensembl_pos_gene_type'
    'ensembl_pos_gene_strand'
    'ensembl_pos_gene_start'
    'ensembl_pos_gene_end'
    """
    # Step 1. Annotate each variant
    data = {
        'variant_id': [],
        'gencode_pos_region': [],
        'gencode_pos_gene_id': [],
        'gencode_pos_gene_name': [],
        'gencode_pos_gene_type': [],
        'gencode_pos_gene_strand': [],
        'gencode_pos_gene_start': [],
        'gencode_pos_gene_end': [],
        'gencode_pos_exon_id': [],
        'gencode_pos_exon_number': []
    }
    for index, row in df_small_variants.iterrows():
        data['variant_id'].append(row['variant_id'])

        # Position 1 annotations
        curr_chrom = row['chrom']
        curr_pos = row['pos']
        curr_pos_region = ''
        curr_pos_gene_id = ''
        curr_pos_gene_name = ''
        curr_pos_gene_type = ''
        curr_pos_gene_strand = ''
        curr_pos_gene_start = ''
        curr_pos_gene_end = ''
        curr_pos_exon_id = ''
        curr_pos_exon_number = ''
        df_gencode_genes_matched = df_gencode_genes.loc[
            (df_gencode_genes['gene_chrom'] == curr_chrom) &
            (df_gencode_genes['gene_start'] <= curr_pos) &
            (df_gencode_genes['gene_end'] >= curr_pos),:
        ]
        if len(df_gencode_genes_matched) > 0:
             df_gencode_exons_matched = df_gencode_exons.loc[
                (df_gencode_exons['exon_chrom'] == curr_chrom) &
                (df_gencode_exons['exon_start'] <= curr_pos) &
                (df_gencode_exons['exon_end'] >= curr_pos), :
             ]
             if len(df_gencode_exons_matched) > 0:
                 curr_pos_region = 'exonic'
                 curr_pos_exon_id = ','.join(df_gencode_exons_matched['exon_id'].values.tolist())
                 curr_pos_exon_number = ','.join(map(str, df_gencode_exons_matched['exon_number'].values.tolist()))
             else:
                 curr_pos_region = 'intronic'
             curr_pos_gene_id = ','.join(df_gencode_genes_matched['gene_id'].values.tolist())
             curr_pos_gene_name = ','.join(df_gencode_genes_matched['gene_name'].values.tolist())
             curr_pos_gene_type = ','.join(df_gencode_genes_matched['gene_type'].values.tolist())
             curr_pos_gene_strand = ','.join(df_gencode_genes_matched['gene_strand'].values.tolist())
             curr_pos_gene_start = ','.join(map(str, df_gencode_genes_matched['gene_start'].values.tolist()))
             curr_pos_gene_end = ','.join(map(str, df_gencode_genes_matched['gene_end'].values.tolist()))

        else:
            curr_pos_region = 'intergenic'
        data['gencode_pos_region'].append(curr_pos_region)
        data['gencode_pos_gene_id'].append(curr_pos_gene_id)
        data['gencode_pos_gene_name'].append(curr_pos_gene_name)
        data['gencode_pos_gene_type'].append(curr_pos_gene_type)
        data['gencode_pos_gene_strand'].append(curr_pos_gene_strand)
        data['gencode_pos_gene_start'].append(curr_pos_gene_start)
        data['gencode_pos_gene_end'].append(curr_pos_gene_end)
        data['gencode_pos_exon_id'].append(curr_pos_exon_id)
        data['gencode_pos_exon_number'].append(curr_pos_exon_number)

    df_annotations = pd.DataFrame(data)
    df_small_variants = pd.merge(df_small_variants, df_annotations, on='variant_id')
    return df_small_variants


def annotate_structural_variants_using_gencode(
        df_structural_variants: pd.DataFrame,
        df_gencode_genes: pd.DataFrame,
        df_gencode_exons: pd.DataFrame
    ) -> pd.DataFrame:
    """
    Annotates structural variants using GENCODE.

    Parameters
    ----------
    df_structural_variants  :   DataFrame of structural variants.
                                Expected columns:
                                'variant_id'
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
                                'sv_type' (DEL, INS, INV, DUP, BND or TRA)
    df_gencode_genes        :   DataFrame of GENCODE genes.
                                Expected columns:
                                'gene_id'
                                'gene_name'
                                'gene_type'
                                'gene_chrom'
                                'gene_start'
                                'gene_end'
                                'gene_strand'
                                'level'
                                'transcripts_count'
    df_gencode_exons        :   DataFrame of GENCODE exons.
                                Expected columns:
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'exon_number'
                                'exon_chrom'
                                'exon_start'
                                'exon_end'

    Returns
    -------
    DataFrame with the following columns appended:
    'ensembl_pos_1_region'
    'ensembl_pos_1_gene_id'
    'ensembl_pos_1_gene_name'
    'ensembl_pos_1_gene_type'
    'ensembl_pos_1_gene_strand'
    'ensembl_pos_1_gene_start'
    'ensembl_pos_1_gene_end'
    'ensembl_pos_2_region'
    'ensembl_pos_2_gene_id'
    'ensembl_pos_2_gene_name'
    'ensembl_pos_2_gene_type'
    'ensembl_pos_2_gene_strand'
    'ensembl_pos_2_gene_start'
    'ensembl_pos_2_gene_end'
    """
    # Step 1. Annotate each variant
    data = {
        'variant_id': [],
        'gencode_pos_1_region': [],
        'gencode_pos_1_gene_id': [],
        'gencode_pos_1_gene_name': [],
        'gencode_pos_1_gene_type': [],
        'gencode_pos_1_gene_strand': [],
        'gencode_pos_1_gene_start': [],
        'gencode_pos_1_gene_end': [],
        'gencode_pos_1_exon_id': [],
        'gencode_pos_1_exon_number': [],
        'gencode_pos_2_region': [],
        'gencode_pos_2_gene_id': [],
        'gencode_pos_2_gene_name': [],
        'gencode_pos_2_gene_type': [],
        'gencode_pos_2_gene_strand': [],
        'gencode_pos_2_gene_start': [],
        'gencode_pos_2_gene_end': [],
        'gencode_pos_2_exon_id': [],
        'gencode_pos_2_exon_number': []
    }
    for index, row in df_structural_variants.iterrows():
        data['variant_id'].append(row['variant_id'])

        # Position 1 annotations
        curr_chr_1 = row['chr_1']
        curr_pos_1 = row['pos_1']
        curr_pos_1_region = ''
        curr_pos_1_exon_id = ''
        curr_pos_1_exon_number = ''
        curr_pos_1_gene_id = ''
        curr_pos_1_gene_name = ''
        curr_pos_1_gene_type = ''
        curr_pos_1_gene_strand = ''
        curr_pos_1_gene_start = ''
        curr_pos_1_gene_end = ''
        df_gencode_genes_matched = df_gencode_genes.loc[
            (df_gencode_genes['gene_chrom'] == curr_chr_1) &
            (df_gencode_genes['gene_start'] <= curr_pos_1) &
            (df_gencode_genes['gene_end'] >= curr_pos_1),:
        ]
        if len(df_gencode_genes_matched) > 0:
             df_gencode_exons_matched = df_gencode_exons.loc[
                (df_gencode_exons['exon_chrom'] == curr_chr_1) &
                (df_gencode_exons['exon_start'] <= curr_pos_1) &
                (df_gencode_exons['exon_end'] >= curr_pos_1), :
            ]
             if len(df_gencode_exons_matched) > 0:
                 curr_pos_1_region = 'exonic'
                 curr_pos_1_exon_id = ','.join(df_gencode_exons_matched['exon_id'].values.tolist())
                 curr_pos_1_exon_number = ','.join(map(str, df_gencode_exons_matched['exon_number'].values.tolist()))
             else:
                 curr_pos_1_region = 'intronic'
             curr_pos_1_gene_id = ','.join(df_gencode_genes_matched['gene_id'].values.tolist())
             curr_pos_1_gene_name = ','.join(df_gencode_genes_matched['gene_name'].values.tolist())
             curr_pos_1_gene_type = ','.join(df_gencode_genes_matched['gene_type'].values.tolist())
             curr_pos_1_gene_strand = ','.join(df_gencode_genes_matched['gene_strand'].values.tolist())
             curr_pos_1_gene_start = ','.join(map(str, df_gencode_genes_matched['gene_start'].values.tolist()))
             curr_pos_1_gene_end = ','.join(map(str, df_gencode_genes_matched['gene_end'].values.tolist()))

        else:
            curr_pos_1_region = 'intergenic'
        data['gencode_pos_1_region'].append(curr_pos_1_region)
        data['gencode_pos_1_gene_id'].append(curr_pos_1_gene_id)
        data['gencode_pos_1_gene_name'].append(curr_pos_1_gene_name)
        data['gencode_pos_1_gene_type'].append(curr_pos_1_gene_type)
        data['gencode_pos_1_gene_strand'].append(curr_pos_1_gene_strand)
        data['gencode_pos_1_gene_start'].append(curr_pos_1_gene_start)
        data['gencode_pos_1_gene_end'].append(curr_pos_1_gene_end)
        data['gencode_pos_1_exon_id'].append(curr_pos_1_exon_id)
        data['gencode_pos_1_exon_number'].append(curr_pos_1_exon_number)

        # Position 2 annotations
        curr_chr_2 = row['chr_2']
        curr_pos_2 = row['pos_2']
        curr_pos_2_region = ''
        curr_pos_2_exon_id = ''
        curr_pos_2_exon_number = ''
        curr_pos_2_gene_id = ''
        curr_pos_2_gene_name = ''
        curr_pos_2_gene_type = ''
        curr_pos_2_gene_strand = ''
        curr_pos_2_gene_start = ''
        curr_pos_2_gene_end = ''
        df_gencode_genes_matched = df_gencode_genes.loc[
            (df_gencode_genes['gene_chrom'] == curr_chr_2) &
            (df_gencode_genes['gene_start'] <= curr_pos_2) &
            (df_gencode_genes['gene_end'] >= curr_pos_2), :
        ]
        if len(df_gencode_genes_matched) > 0:
            df_gencode_exons_matched = df_gencode_exons.loc[
                (df_gencode_exons['exon_chrom'] == curr_chr_2) &
                (df_gencode_exons['exon_start'] <= curr_pos_2) &
                (df_gencode_exons['exon_end'] >= curr_pos_2), :
            ]
            if len(df_gencode_exons_matched) > 0:
                curr_pos_2_region = 'exonic'
                curr_pos_2_exon_id = ','.join(df_gencode_exons_matched['exon_id'].values.tolist())
                curr_pos_2_exon_number = ','.join(map(str, df_gencode_exons_matched['exon_number'].values.tolist()))
            else:
                curr_pos_2_region = 'intronic'
            curr_pos_2_gene_id = ','.join(df_gencode_genes_matched['gene_id'].values.tolist())
            curr_pos_2_gene_name = ','.join(df_gencode_genes_matched['gene_name'].values.tolist())
            curr_pos_2_gene_type = ','.join(df_gencode_genes_matched['gene_type'].values.tolist())
            curr_pos_2_gene_strand = ','.join(df_gencode_genes_matched['gene_strand'].values.tolist())
            curr_pos_2_gene_start = ','.join(map(str, df_gencode_genes_matched['gene_start'].values.tolist()))
            curr_pos_2_gene_end = ','.join(map(str, df_gencode_genes_matched['gene_end'].values.tolist()))
        else:
            curr_pos_2_region = 'intergenic'
        data['gencode_pos_2_region'].append(curr_pos_2_region)
        data['gencode_pos_2_gene_id'].append(curr_pos_2_gene_id)
        data['gencode_pos_2_gene_name'].append(curr_pos_2_gene_name)
        data['gencode_pos_2_gene_type'].append(curr_pos_2_gene_type)
        data['gencode_pos_2_gene_strand'].append(curr_pos_2_gene_strand)
        data['gencode_pos_2_gene_start'].append(curr_pos_2_gene_start)
        data['gencode_pos_2_gene_end'].append(curr_pos_2_gene_end)
        data['gencode_pos_2_exon_id'].append(curr_pos_2_exon_id)
        data['gencode_pos_2_exon_number'].append(curr_pos_2_exon_number)

    df_annotations = pd.DataFrame(data)
    df_structural_variants = pd.merge(df_structural_variants, df_annotations, on='variant_id')
    return df_structural_variants
