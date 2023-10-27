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
The purpose of this python3 script is to implement the Ensembl dataclass.
"""


import pyensembl
from dataclasses import dataclass
from typing import List, Tuple
from .annotation_db import AnnotationDb
from .constants import *
from .logging import get_logger
from .variant_annotation import VariantAnnotation
from .variant_call import VariantCall
from .variants_list import VariantsList


logger = get_logger(__name__)


@dataclass
class Ensembl(AnnotationDb):
    release: int
    species: str
    _ensembl = None

    @property
    def source(self):
        return AnnotationSources.ENSEMBL

    @property
    def ensembl(self):
        if self._ensembl is None:
            self._ensembl = pyensembl.EnsemblRelease(release=self.release, species=self.species)
        return self._ensembl

    def annotate_position_using_pyensembl(
            self,
            chromosome: str,
            position: int
    ) -> List[VariantAnnotation]:
        """
        Annotates a position using pyensembl and returns a list of
        VariantAnnotation objects.

        Parameters
        ----------
        chromosome              :   Chromosome.
        position                :   Position.

        Returns
        -------
        variant_annotations     :   List of VariantAnnotation objects.
        """
        variant_annotations = []
        chromosome = chromosome.replace('chr', '')
        genes = self.ensembl.genes_at_locus(contig=chromosome, position=position)
        if len(genes) == 0:
            variant_annotation = VariantAnnotation(
                region=GenomicRegionTypes.INTERGENIC,
                source=AnnotationSources.ENSEMBL,
                source_version=str(self.release),
                gene_id='',
                gene_stable_id='',
                gene_version='',
                gene_name='',
                gene_type='',
                gene_strand='',
                species=self.species
            )
            variant_annotations.append(variant_annotation)
        else:
            for gene in genes:
                region = GenomicRegionTypes.INTRONIC
                exon_ids = self.ensembl.exon_ids_of_gene_id(gene.gene_id)
                for exon_id in exon_ids:
                    exon = self.ensembl.exon_by_id(exon_id=exon_id)
                    if exon.start <= position <= exon.end:
                        region = GenomicRegionTypes.EXONIC
                variant_annotation = VariantAnnotation(
                    region=region,
                    source=AnnotationSources.ENSEMBL,
                    source_version=str(self.release),
                    gene_id=gene.gene_id,
                    gene_stable_id=gene.gene_id,
                    gene_version='',
                    gene_name=gene.gene_name,
                    gene_type=gene.biotype,
                    gene_strand=gene.strand,
                    species=self.species
                )
                variant_annotations.append(variant_annotation)
        return variant_annotations

    def annotate_variant_call_using_pyensembl(
            self,
            variant_call: VariantCall
    ) -> Tuple[List[VariantAnnotation], List[VariantAnnotation]]:
        """
        Annotates a VariantCall object and returns two lists of VariantAnnotation objects.

        Parameters
        ----------
        variant_call            :   VariantCall object.

        Returns
        -------
        position_1_annotations  :   List of VariantAnnotation objects for position 1.
        position_2_annotations  :   List of VariantAnnotation objects for position 2.
        """
        position_1_annotations = self.annotate_position_using_pyensembl(
            chromosome=variant_call.chromosome_1,
            position=variant_call.position_1
        )
        position_2_annotations = self.annotate_position_using_pyensembl(
            chromosome=variant_call.chromosome_2,
            position=variant_call.position_2
        )
        return position_1_annotations, position_2_annotations

    def annotate_variants_list(self, variants_list) -> VariantsList:
        """
        Annotates a VariantsList object.

        Parameters
        ----------
        variants_list   :   VariantsList object.

        Returns
        -------
        variants_list   :   VariantsList object.
        """
        for i in range(0, len(variants_list.variants)):
            for j in range(0, len(variants_list.variants[i].variant_calls)):
                position_1_annotations, position_2_annotations = self.annotate_variant_call_using_pyensembl(
                    variants_list.variants[i].variant_calls[j]
                )
                for annotation in position_1_annotations:
                    variants_list.variants[i].variant_calls[j].position_1_annotations.append(annotation)
                for annotation in position_2_annotations:
                    variants_list.variants[i].variant_calls[j].position_2_annotations.append(annotation)
        return variants_list

