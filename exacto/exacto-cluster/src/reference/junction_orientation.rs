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


use bimap::BiMap;
use exacto_caller::prelude::SpliceJunction;
use exacto_core::prelude::{FastaMap, Strand};


pub(crate) fn identify_splice_motif_strand(
    intron: (u16, u32, u32),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap
) -> Option<Strand> {
    let (chromosome_id, start, end): (u16, u32, u32) = intron;
    let chromosome: &Box<str> = chromosome_names_map.get_by_right(&chromosome_id)?;
    if end < start.saturating_add(3) {
        return None;
    }
    let donor: String = fasta_map.try_get_sequence(chromosome, start as usize, start as usize + 1)?.to_uppercase();
    let acceptor: String = fasta_map.try_get_sequence(chromosome, end as usize - 1, end as usize)?.to_uppercase();
    match (donor.as_str(), acceptor.as_str()) {
        ("GT", "AG") | ("GC", "AG") | ("AT", "AC") => Some(Strand::Forward),
        ("CT", "AC") | ("CT", "GC") | ("GT", "AT") => Some(Strand::Reverse),
        _ => None
    }
}


pub(crate) fn is_aligned_against_transcript(
    splice_junctions: &[SpliceJunction],
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap
) -> bool {
    let (mut num_along, mut num_against): (usize, usize) = (0, 0);
    for junction in splice_junctions.iter().filter(|junction| junction.chromosome_1 == junction.chromosome_2) {
        match identify_splice_motif_strand(junction.intron_span_key(), chromosome_names_map, fasta_map) {
            Some(strand) if strand == junction.strand_1 => num_along += 1,
            Some(_) => num_against += 1,
            None => {}
        }
    }
    num_against > num_along
}


pub(crate) fn mirror_splice_junction(junction: &SpliceJunction) -> SpliceJunction {
    let flip = |strand: &Strand| -> Strand {
        match strand {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            other => other.clone()
        }
    };
    SpliceJunction::new(
        junction.chromosome_2,
        junction.chromosome_1,
        junction.position_2,
        junction.position_1,
        flip(&junction.strand_2),
        flip(&junction.strand_1)
    )
}
