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


use polars::prelude::*;

use crate::io::records::StitchedTranscriptRecord;


pub fn stitched_transcript_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = StitchedTranscriptRecord>
{
    let mut read_name: Vec<String> = Vec::new();
    let mut stitched_sequence: Vec<String> = Vec::new();
    let mut read_length: Vec<u32> = Vec::new();
    let mut stitched_length: Vec<u32> = Vec::new();
    let mut retained_read_start: Vec<u32> = Vec::new();
    let mut retained_read_end: Vec<u32> = Vec::new();
    let mut is_five_prime_stitched: Vec<bool> = Vec::new();
    let mut five_prime_reference_gene_id: Vec<String> = Vec::new();
    let mut five_prime_reference_transcript_id: Vec<String> = Vec::new();
    let mut five_prime_stitch_start: Vec<u32> = Vec::new();
    let mut five_prime_stitch_end: Vec<u32> = Vec::new();
    let mut five_prime_stitch_length: Vec<u32> = Vec::new();
    let mut is_three_prime_stitched: Vec<bool> = Vec::new();
    let mut three_prime_reference_gene_id: Vec<String> = Vec::new();
    let mut three_prime_reference_transcript_id: Vec<String> = Vec::new();
    let mut three_prime_stitch_start: Vec<u32> = Vec::new();
    let mut three_prime_stitch_end: Vec<u32> = Vec::new();
    let mut three_prime_stitch_length: Vec<u32> = Vec::new();
    let mut unmodeled_reason: Vec<String> = Vec::new();

    for r in records {
        read_name.push(r.read_name.into());
        stitched_sequence.push(r.stitched_sequence.into());
        read_length.push(r.read_length);
        stitched_length.push(r.stitched_length);
        retained_read_start.push(r.retained_read_start);
        retained_read_end.push(r.retained_read_end);
        is_five_prime_stitched.push(r.is_five_prime_stitched);
        five_prime_reference_gene_id.push(r.five_prime_reference_gene_id.into());
        five_prime_reference_transcript_id.push(r.five_prime_reference_transcript_id.into());
        five_prime_stitch_start.push(r.five_prime_stitch_start);
        five_prime_stitch_end.push(r.five_prime_stitch_end);
        five_prime_stitch_length.push(r.five_prime_stitch_length);
        is_three_prime_stitched.push(r.is_three_prime_stitched);
        three_prime_reference_gene_id.push(r.three_prime_reference_gene_id.into());
        three_prime_reference_transcript_id.push(r.three_prime_reference_transcript_id.into());
        three_prime_stitch_start.push(r.three_prime_stitch_start);
        three_prime_stitch_end.push(r.three_prime_stitch_end);
        three_prime_stitch_length.push(r.three_prime_stitch_length);
        unmodeled_reason.push(r.unmodeled_reason.into());
    }

    DataFrame::new(vec![
        Column::from(Series::new("read_name".into(), read_name)),
        Column::from(Series::new("stitched_sequence".into(), stitched_sequence)),
        Column::from(Series::new("read_length".into(), read_length)),
        Column::from(Series::new("stitched_length".into(), stitched_length)),
        Column::from(Series::new("retained_read_start".into(), retained_read_start)),
        Column::from(Series::new("retained_read_end".into(), retained_read_end)),
        Column::from(Series::new("is_five_prime_stitched".into(), is_five_prime_stitched)),
        Column::from(Series::new("five_prime_reference_gene_id".into(), five_prime_reference_gene_id)),
        Column::from(Series::new("five_prime_reference_transcript_id".into(), five_prime_reference_transcript_id)),
        Column::from(Series::new("five_prime_stitch_start".into(), five_prime_stitch_start)),
        Column::from(Series::new("five_prime_stitch_end".into(), five_prime_stitch_end)),
        Column::from(Series::new("five_prime_stitch_length".into(), five_prime_stitch_length)),
        Column::from(Series::new("is_three_prime_stitched".into(), is_three_prime_stitched)),
        Column::from(Series::new("three_prime_reference_gene_id".into(), three_prime_reference_gene_id)),
        Column::from(Series::new("three_prime_reference_transcript_id".into(), three_prime_reference_transcript_id)),
        Column::from(Series::new("three_prime_stitch_start".into(), three_prime_stitch_start)),
        Column::from(Series::new("three_prime_stitch_end".into(), three_prime_stitch_end)),
        Column::from(Series::new("three_prime_stitch_length".into(), three_prime_stitch_length)),
        Column::from(Series::new("unmodeled_reason".into(), unmodeled_reason))
    ]).unwrap()
}
