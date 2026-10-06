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


use exacto_core::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::prelude::{RNAReadCharacterizationSummary, SplicingEventReadCounts};


pub struct RNAJunctionReadCounts {
    /// Vec<(chromosome ID, first intron base, last intron base)>
    introns: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)>,

    /// Vec<(chromosome ID, last intron base, first intron base)>
    introns_by_end: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)>,

    num_spliced: Vec<ReadSupport>,
    num_retained: Vec<ReadSupport>
}

impl RNAJunctionReadCounts {
    pub fn new(summaries: &[RNAReadCharacterizationSummary], num_threads: usize) -> Self {
        let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        // Step 1. Count the reads splicing each intron.
        // HashMap<(ReferenceChromosomeID, ReferencePosition, ReferencePosition), number of reads>
        let junction_counts: HashMap<(ReferenceChromosomeID, ReferencePosition, ReferencePosition), ReadSupport> = thread_pool.install(|| {
            summaries
                .par_iter()
                .fold(
                    HashMap::new,
                    |mut counts: HashMap<(ReferenceChromosomeID, ReferencePosition, ReferencePosition), ReadSupport>, summary| {
                        let spliced: BTreeSet<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = summary
                            .splice_junctions
                            .iter()
                            .filter(|junction| junction.chromosome_1 == junction.chromosome_2)
                            .map(|junction| junction.intron_span_key())
                            .collect();
                        for intron in spliced {
                            *counts.entry(intron).or_insert(0) += 1;
                        }
                        counts
                    }
                )
                .reduce(HashMap::new, |counts_1, counts_2| {
                    // Merge the smaller map into the larger one.
                    let (mut larger, smaller) = if counts_1.len() >= counts_2.len() { (counts_1, counts_2) } else { (counts_2, counts_1) };
                    for (intron, count) in smaller {
                        *larger.entry(intron).or_insert(0) += count;
                    }
                    larger
                })
        });
        let mut introns: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = junction_counts.keys().copied().collect();
        introns.sort_unstable();
        let num_spliced: Vec<ReadSupport> = introns.iter().map(|intron| junction_counts[intron]).collect();

        // Step 2. Count the reads retaining each intron.
        let mut introns_by_end: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = introns
            .iter()
            .map(|&(chromosome, intron_start, intron_end)| (chromosome, intron_end, intron_start))
            .collect();
        introns_by_end.sort_unstable();
        let num_retained: Vec<ReadSupport> = thread_pool.install(|| {
            summaries
                .par_iter()
                .fold(
                    || vec![0; introns.len()],
                    |mut counts: Vec<ReadSupport>, summary| {
                        for index in summary.identify_retained_introns(&introns, &introns_by_end) {
                            counts[index] += 1;
                        }
                        counts
                    }
                )
                .reduce(
                    || vec![0; introns.len()],
                    |mut counts_1, counts_2| {
                        counts_1.iter_mut().zip(counts_2).for_each(|(count_1, count_2)| *count_1 += count_2);
                        counts_1
                    }
                )
        });

        Self {
            introns: introns,
            introns_by_end: introns_by_end,
            num_spliced: num_spliced,
            num_retained: num_retained
        }
    }

    pub fn get_depths(&self) -> HashSet<ReadDepth> {
        self.num_spliced
            .iter()
            .zip(self.num_retained.iter())
            .map(|(num_spliced, num_retained)| num_spliced + num_retained)
            .collect()
    }
}

impl SplicingEventReadCounts for RNAJunctionReadCounts {
    type Input = RNAReadCharacterizationSummary;

    fn get_read_counts(&self, summary: &RNAReadCharacterizationSummary) -> Vec<(ReadSupport, ReadDepth)> {
        summary
            .identify_retained_introns(&self.introns, &self.introns_by_end)
            .into_iter()
            .map(|index| (self.num_retained[index], self.num_spliced[index] + self.num_retained[index]))
            .collect()
    }
}


#[cfg(test)]
#[path = "../../tests/filtering/rna/junction_read_counts.rs"]
mod tests;