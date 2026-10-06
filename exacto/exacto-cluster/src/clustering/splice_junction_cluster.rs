// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use std::collections::HashSet;


#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SpliceJunctionCluster {
    pub(crate) id: usize,

    /// Splice-junction chain in genomic order.
    /// Empty splice_junctions means this is a single-exon (junction-less) gene cluster
    pub(crate) splice_junctions: Vec<SpliceJunction>,

    /// Reads in this cluster. Each read is committed to exactly one cluster
    /// (its single best backbone chain, or single best single-exon gene).
    pub(crate) read_ids: HashSet<ReadID>,

    /// Reference gene and transcript IDs
    pub(crate) reference_gene_transcript_ids: Vec<(ReferenceGeneID, ReferenceTranscriptID)>
}

impl SpliceJunctionCluster {
    pub(crate) fn new(
        id: usize,
        splice_junctions: Vec<SpliceJunction>,
        read_ids: HashSet<ReadID>,
        reference_gene_transcript_ids: Vec<(ReferenceGeneID, ReferenceTranscriptID)>
    ) -> Self {
        Self { 
            id,
            splice_junctions,
            read_ids,
            reference_gene_transcript_ids
        }
    }
}
