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
The purpose of this python3 script is to implement the Ensembl class.
"""


import pyensembl
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .constants import *
from .variant_call import VariantCall
from .variant_annotation import VariantAnnotation
from .variants_list import VariantsList
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class Ensembl:
    release: int = None
    species: str = None
    ensembl_txt_file: str = None
    use_pyensembl: bool = True

    @property
    def ensembl(self):
        return pyensembl.EnsemblRelease(release=self.release, species=self.species)

    def annotate_variant_call_using_pyensembl(
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
        chr_1 = variant_call.chr_1
        chr_1 = chr_1.replace('chr', '')
        pos_1_genes = self.ensembl.genes_at_locus(contig=chr_1, position=variant_call.pos_1)
        if len(pos_1_genes) == 0:
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chr_1,
                pos=variant_call.pos_1,
                region=GenomicRegionTypes.INTERGENIC
            )
            pos_1_annotations.append(variant_annotation)
        else:
            is_exonic = False
            for gene in pos_1_genes:
                exons = self.ensembl.exon_ids_of_gene_id(gene.gene_id)
                for exon in exons:
                    if exon.start <= variant_call.pos_1 <= exon.end:
                        is_exonic = True
            if is_exonic:
                region = GenomicRegionTypes.EXONIC
            else:
                region = GenomicRegionTypes.INTRONIC
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chr_1,
                pos=variant_call.pos_1,
                region=region
            )
            pos_1_annotations.append(variant_annotation)

        # Position 2
        pos_2_annotations = []
        chr_2 = variant_call.chr_2
        chr_2 = chr_2.replace('chr', '')
        pos_2_genes = self.ensembl.genes_at_locus(contig=chr_2, position=variant_call.pos_2)
        if len(pos_2_genes) == 0:
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chr_2,
                pos=variant_call.pos_2,
                region=GenomicRegionTypes.INTERGENIC
            )
            pos_2_annotations.append(variant_annotation)
        else:
            is_exonic = False
            for gene in pos_2_genes:
                exons = self.ensembl.exon_ids_of_gene_id(gene.gene_id)
                for exon in exons:
                    if exon.start <= variant_call.pos_2 <= exon.end:
                        is_exonic = True
            if is_exonic:
                region = GenomicRegionTypes.EXONIC
            else:
                region = GenomicRegionTypes.INTRONIC
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chr_2,
                pos=variant_call.pos_2,
                region=region
            )
            pos_2_annotations.append(variant_annotation)

        return pos_1_annotations, pos_2_annotations

    def annotate_variants(self, variants_list) -> VariantsList:
        """
        Annotates variants using Ensembl data.

        Parameters
        ----------
        variants_list   :   An instance of the VariantsList class.

        Returns
        -------
        variants_list   :   An instance of the VariantsList class with
                            each variant_call annotated.
        """
        # Check if member variables for pyensembl have been set
        if self.use_pyensembl:
            if self.release is None or self.species is None:
                logger.error('release and species must be set to query pyensembl.')
                exit(1)

        if self.use_pyensembl:
            for i in range(0, len(variants_list.variant_ids)):
                for j in range(0, len(variants_list.variants[i].variant_calls)):
                    pos_1_annotations, pos_2_annotations = self.annotate_variant_call_using_pyensembl(variants_list.variants[i].variant_calls[j])
                    variants_list.variants[i].variant_calls[j].pos_1_annotations = pos_1_annotations
                    variants_list.variants[i].variant_calls[j].pos_2_annotations = pos_2_annotations
        else:
            # todo: implement variant annotation using Ensembl TXT file
            a = 1
        return variants_list
