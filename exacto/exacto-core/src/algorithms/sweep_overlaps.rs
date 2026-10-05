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


use std::collections::BTreeMap;


pub fn sweep_overlaps(segments: &Vec<(u32, u32)>) -> Vec<(usize, usize)> {
    // Step 1. Get Vec<(start, end, original_index)>
    let mut segs: Vec<(u32, u32, usize)> = segments
        .iter()
        .enumerate()
        .map(|(idx, &(s, e))| (s, e, idx))
        .collect();

    // Step 2. Sort by start ascending, tie by end ascending
    segs.sort_by_key(|&(s, e, _)| (s, e));

    // Step 3. Active set keyed by end -> indices of segments with that end
    let mut active: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

    // Step 4. Get pairs of overlapping segments
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for (s, e, idx) in segs {
        // Remove segments with end < current start
        let expired_ends: Vec<u32> = active
            .range(..s) // keys < s
            .map(|(end, _)| *end)
            .collect();

        for end in expired_ends {
            active.remove(&end);
        }

        for (_end, ids) in active.iter() {
            for &other in ids.iter() {
                let (i, j) = if other < idx { (other, idx) } else { (idx, other) };
                pairs.push((i, j));
            }
        }

        // Add current segment to active
        active.entry(e).or_default().push(idx);
    }

    pairs
}
