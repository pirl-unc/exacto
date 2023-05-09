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
    gene_id: str = None
    gene_stable_id: str = None
    gene_version: str = None
    gene_name: str = None
    gene_type: str = None
    gene_source: str = None
    gene_strand: str = None

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
            self.gene_strand not in [Strands.POSITIVE, Strands.NEGATIVE]:
            raise Exception(
                'gene_strand must be one of the following: %s' %
                (', '.join([Strands.POSITIVE, Strands.NEGATIVE]))
            )

    def to_dict(self):
        data = {
            'region': [self.region],
            'source': [self.source],
            'gene_id': ['' if self.gene_id is None else self.gene_id],
            'gene_stable_id': ['' if self.gene_stable_id is None else self.gene_stable_id],
            'gene_version': ['' if self.gene_version is None else self.gene_version],
            'gene_name': ['' if self.gene_name is None else self.gene_name],
            'gene_type': ['' if self.gene_type is None else self.gene_type],
            'gene_source': ['' if self.gene_source is None else self.gene_source],
            'gene_strand': ['' if self.gene_strand is None else self.gene_strand]
        }
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())

