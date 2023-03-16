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
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .gene import Gene
from .exon import Exon


@dataclass
class VariantAnnotation:
    chrom: str = None
    pos: int = None
    region: str = None
    gene: Gene = None

    @property
    def exons(self) -> List[Exon]:
        exons = []
        for transcript in self.gene.transcripts:
            for exon in transcript.exons:
                if exon.start <= self.pos <= exon.end:
                    exons.append(exon)
        return exons

    @property
    def exon_ids(self) -> List[str]:
        exon_ids = []
        for transcript in self.gene.transcripts:
            for exon in transcript.exons:
                if exon.start <= self.pos <= exon.end:
                    exon_ids.append(exon.id)
        return exon_ids

    def to_dataframe(self) -> pd.DataFrame:
        df = pd.DataFrame({
            'chrom': [self.chrom],
            'pos': [self.pos],
            'region': [self.region],
        })
        if self.gene is not None:
            return pd.concat([df, self.gene.to_dataframe()], axis=1)
        else:
            return df

