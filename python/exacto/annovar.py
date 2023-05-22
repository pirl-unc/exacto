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
The purpose of this python3 script is to implement the Annovar dataclass.
"""


import pandas as pd
import subprocess as sp
from dataclasses import dataclass
from .annotation_db import AnnotationDb
from .constants import AnnotationSources
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class Annovar(AnnotationDb):
    perl_path: str
    annovar_path: str
    humandb_path: str
    genome_assembly: str
    protocol: str
    operation: str

    @property
    def source(self):
        return AnnotationSources.ANNOVAR

    def write_annovar_input_file(
            self,
            tsv_file: str,
            annovar_input_file: str
    ):
        """
        Writes an ANNOVAR input file from a TSV file.

        Parameters
        ----------
        tsv_file            :   TSV file. Expected columns:
                                'chrom'
                                'pos'
                                'ref'
                                'alt'
                                'genotype'
        annovar_input_file  :   ANNOVAR input file.
        """
        df = pd.read_csv(tsv_file, sep='\t')
        data = {
            'chrom': [],
            'start': [],
            'end': [],
            'ref': [],
            'alt': [],
            'genotype': []
        }
        for row in df.itertuples(index=False):
            if len(row.ref) == 1 and len(row.alt) == 1:  # SNV
                data['chrom'].append(row.chrom.replace('chr', ''))
                data['start'].append(row.pos)
                data['end'].append(row.pos)
                data['ref'].append(row.ref)
                data['alt'].append(row.alt)
            elif len(row.ref) > 1 and len(row.alt) == 1:  # Deletion
                data['chrom'].append(row.chrom.replace('chr', ''))
                data['start'].append(row.pos + 1)
                data['end'].append(row.pos + len(row.ref) - 1)
                data['ref'].append(row.ref[1:])
                data['alt'].append('-')
            elif len(row.ref) == 1 and len(row.alt) > 1:  # Insertion
                data['chrom'].append(row.chrom.replace('chr', ''))
                data['start'].append(row.pos)
                data['end'].append(row.pos)
                data['ref'].append('-')
                data['alt'].append(row.alt[1:])
            else:
                data['chrom'].append(row.chrom.replace('chr', ''))
                data['start'].append(row.pos)
                data['end'].append(row.pos)
                data['ref'].append(row.ref)
                data['alt'].append(row.alt)

            if row.genotype == '0/1' or row.genotype == '0|1' or row.genotype == '1|0' or row.genotype == '0|1':
                data['genotype'].append('het')
            elif row.genotype == '0/0' or row.genotype == '1/1' or row.genotype == '1|1' or row.genotype == '0|0':
                data['genotype'].append('hom')
            else:
                data['genotype'].append('-')

        df_annovar_input = pd.DataFrame(data)
        df_annovar_input.drop_duplicates(inplace=True)
        df_annovar_input.to_csv(annovar_input_file, sep="\t", header=False, index=False)

    def annotate_small_variants_using_annovar(
            self,
            annovar_input_file: str,
            output_file: str
    ) -> pd.DataFrame:
        """
        Annotates small variants using ANNOVAR.

        Parameters
        ----------
        perl_path               :   perl path.
        annovar_path            :   ANNOVAR path.
        humandb_path            :   humandb path.
        annovar_input_file      :   ANNOVAR input file.
        genome_assembly         :   Genome assembly (e.g. 'hg38').
        protocol                :   Protocol (e.g. 'refGene,exac03,1000g2015aug_eur,1000g2015aug_eas,1000g2015aug_sas,clinvar_20210501,cosmic95_coding,avsnp150,dbnsfp42c').
        operation               :   Operation (e.g. 'g,f,f,f,f,f,f,f,f').
        output_file             :   Output file.

        Returns
        -------
        DataFrame of the output ANNOVAR .txt file
        """
        # Step 1. Run ANNOVAR
        # perl table_annovar.pl
        # --buildver <genome_assembly>
        # --protocol <protocol>
        # --operation <operation>
        # --outfile <output_file>
        # --remove
        # --otherinfo
        # <annovar_input_file>
        # <humandb_path>
        cmd = [self.perl_path, self.annovar_path + 'table_annovar.pl',
               '--buildver', self.genome_assembly,
               '--protocol', self.protocol,
               '--operation', self.operation,
               '--outfile', output_file,
               '--remove',
               '--otherinfo',
               annovar_input_file,
               self.humandb_path]
        cmd = ' '.join(cmd)
        sp.call(cmd, shell=True)

        # Step 2. Load output file
        df_annotations = pd.read_csv(output_file + '.' + self.genome_assembly + '_multianno.txt', sep='\t')
        return df_annotations

    def annotate_variants_list(self):
        # todo implement
        pass
