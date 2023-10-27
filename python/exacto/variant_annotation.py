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
The purpose of this python3 script is to implement VariantAnnotation dataclass.
"""


import pandas as pd
from dataclasses import dataclass
from .constants import GenomicRegionTypes, AnnotationSources, Strands


@dataclass(frozen=True)
class VariantAnnotation:
    region: str
    source: str
    source_version: str
    gene_id: str
    gene_stable_id: str
    gene_version: str
    gene_name: str
    gene_type: str
    gene_strand: str
    species: str

    def __post_init__(self):
        if self.region not in GenomicRegionTypes.ALL:
            raise Exception(
                'region must be one of the following: %s' %
                (', '.join(GenomicRegionTypes.ALL))
            )
        if self.source not in AnnotationSources.ALL:
            raise Exception(
                'source must be one of the following: %s' %
                (', '.join(AnnotationSources.ALL))
            )
        if self.gene_strand is not None and \
            self.gene_strand not in [Strands.POSITIVE, Strands.NEGATIVE, '']:
            raise Exception(
                'gene_strand must be one of the following: %s' %
                (', '.join([Strands.POSITIVE, Strands.NEGATIVE, '']))
            )

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())

    def to_dataframe_row(self):
        return {
            'region': [self.region],
            'source': [self.source],
            'source_version': [self.source_version],
            'gene_id': [self.gene_id],
            'gene_stable_id': [self.gene_stable_id],
            'gene_version': [self.gene_version],
            'gene_name': [self.gene_name],
            'gene_type': [self.gene_type],
            'gene_strand': [self.gene_strand],
            'species': [self.species]
        }

    def to_dict(self):
        return {
            'region': self.region,
            'source': self.source,
            'source_version': self.source_version,
            'gene_id': self.gene_id,
            'gene_stable_id': self.gene_stable_id,
            'gene_version': self.gene_version,
            'gene_name': self.gene_name,
            'gene_type': self.gene_type,
            'gene_strand': self.gene_strand,
            'species': self.species
        }
