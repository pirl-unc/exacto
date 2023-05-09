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
from .annotation import Annotation
from .constants import *
from .gene import Gene
from .variant_call import VariantCall
from .variant_annotation import VariantAnnotation
from .variants_list import VariantsList
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class Ensembl(Annotation):
    release: int = None
    species: str = None
    ensembl_txt_file: str = None
    use_pyensembl: bool = True
    _ensembl = None

    @property
    def ensembl(self):
        if self._ensembl is None:
            self._ensembl = pyensembl.EnsemblRelease(release=self.release, species=self.species)
        return self._ensembl

    def annotate_variant_call_using_pyensembl(self, variant_call: VariantCall) -> Tuple[List[VariantAnnotation], List[VariantAnnotation]]:
        """
        Annotates an instance of the VariantCall class and returns two lists of VariantAnnotation instances.

        Parameters
        ----------
        variant_call            :   VariantCall object.

        Returns
        -------
        position_1_annotations  :   List of VariantAnnotation objects for position 1.
        position_2_annotations  :   List of VariantAnnotation objects for position 2.
        """
        # Position 1
        position_1_annotations = []
        chromosome_1 = variant_call.chromosome_1
        chromosome_1 = chromosome_1.replace('chr', '')
        position_1_genes = self.ensembl.genes_at_locus(contig=chromosome_1, position=variant_call.position_1)
        if len(position_1_genes) == 0:
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chromosome_1,
                pos=variant_call.position_1,
                region=GenomicRegionTypes.INTERGENIC
            )
            position_1_annotations.append(variant_annotation)
        else:
            for position_1_gene in position_1_genes:
                gene = Gene(
                    id=position_1_gene.gene_id,
                    source=AnnotationSources.ENSEMBL,
                    name=position_1_gene.gene_name,
                    chromosome=chromosome_1,
                    start=position_1_gene.start,
                    end=position_1_gene.end,
                    strand=position_1_gene.strand,
                    type=position_1_gene.biotype
                )
                is_exonic = False
                exon_ids = self.ensembl.exon_ids_of_gene_id(position_1_gene.gene_id)
                for exon_id in exon_ids:
                    exon = self.ensembl.exon_by_id(exon_id=exon_id)
                    if exon.start <= variant_call.position_1 <= exon.end:
                        is_exonic = True
                if is_exonic:
                    region = GenomicRegionTypes.EXONIC
                else:
                    region = GenomicRegionTypes.INTRONIC
                variant_annotation = VariantAnnotation(
                    chrom=variant_call.chromosome_1,
                    pos=variant_call.position_1,
                    region=region,
                    gene=gene
                )
                position_1_annotations.append(variant_annotation)

        # Position 2
        position_2_annotations = []
        chromosome_2 = variant_call.chromosome_2
        chromosome_2 = chromosome_2.replace('chr', '')
        position_2_genes = self.ensembl.genes_at_locus(contig=chromosome_2, position=variant_call.position_2)
        if len(position_2_genes) == 0:
            variant_annotation = VariantAnnotation(
                chrom=variant_call.chromosome_2,
                pos=variant_call.position_2,
                region=GenomicRegionTypes.INTERGENIC
            )
            position_2_annotations.append(variant_annotation)
        else:
            for position_2_gene in position_2_genes:
                gene = Gene(
                    id=position_2_gene.gene_id,
                    source=AnnotationSources.ENSEMBL,
                    name=position_2_gene.gene_name,
                    chromosome=chromosome_2,
                    start=position_2_gene.start,
                    end=position_2_gene.end,
                    strand=position_2_gene.strand,
                    type=position_2_gene.biotype
                )
                is_exonic = False
                exon_ids = self.ensembl.exon_ids_of_gene_id(position_2_gene.gene_id)
                for exon_id in exon_ids:
                    exon = self.ensembl.exon_by_id(exon_id=exon_id)
                    if exon.start <= variant_call.position_2 <= exon.end:
                        is_exonic = True
                if is_exonic:
                    region = GenomicRegionTypes.EXONIC
                else:
                    region = GenomicRegionTypes.INTRONIC
                variant_annotation = VariantAnnotation(
                    chrom=variant_call.chromosome_2,
                    pos=variant_call.position_2,
                    region=region,
                    gene=gene
                )
                position_2_annotations.append(variant_annotation)

        return position_1_annotations, position_2_annotations

    def annotate_variants(self, variants_list) -> VariantsList:
        """
        Annotates variants using Ensembl data.

        Parameters
        ----------
        variants_list   :   VariantsList object.

        Returns
        -------
        variants_list   :   VariantsList object.
        """
        # Check if member variables for pyensembl have been set
        if self.use_pyensembl:
            if self.release is None or self.species is None:
                logger.error('release and species must be set to query pyensembl.')
                exit(1)

        if self.use_pyensembl:
            for i in range(0, variants_list.size):
                for j in range(0, variants_list.variants[i].size):
                    position_1_annotations, position_2_annotations = self.annotate_variant_call_using_pyensembl(
                        variants_list.variants[i].variant_calls[j]
                    )
                    for annotation in position_1_annotations:
                        variants_list.variants[i].variant_calls[j].position_1_annotations.append(annotation)
                    for annotation in position_2_annotations:
                        variants_list.variants[i].variant_calls[j].position_2_annotations.append(annotation)
        else:
            # todo: implement variant annotation using Ensembl TXT file
            a = 1
        return variants_list
