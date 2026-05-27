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


pub fn join_box_strs(items: &Vec<Box<str>>) -> Box<str> {
    items.iter().map(|s| s.as_ref()).collect::<Vec<&str>>()
        .join(",").into_boxed_str()
}

pub fn join_u32_ids(ids: &[u32]) -> String {
    ids.iter().map(|id| id.to_string()).collect::<Vec<String>>().join(",")
}

/// Collapse a sorted slice of indices into a compact interval string.
/// Example: `[3, 4, 5, 12, 20, 21]` -> `"3-5,12,20-21"`.
/// Input must be sorted ascending; duplicates are tolerated.
pub fn format_intervals(indices: &[u32]) -> String {
    if indices.is_empty() {
        return String::new();
    }
    let mut parts: Vec<String> = Vec::new();
    let mut run_start: u32 = indices[0];
    let mut run_end: u32 = indices[0];
    for &i in indices.iter().skip(1) {
        if i == run_end || i == run_end + 1 {
            run_end = i;
        } else {
            parts.push(if run_start == run_end {
                run_start.to_string()
            } else {
                format!("{}-{}", run_start, run_end)
            });
            run_start = i;
            run_end = i;
        }
    }
    parts.push(if run_start == run_end {
        run_start.to_string()
    } else {
        format!("{}-{}", run_start, run_end)
    });
    parts.join(",")
}