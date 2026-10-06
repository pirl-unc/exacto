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

use crate::prelude::*;
use crate::calling::alignment_model_record_calling::*;


pub struct DNAVariantRecordCaller {
    min_mapping_quality: MappingQuality,
    min_base_quality: BaseQuality,
    min_terminal_soft_clip_ins_len: u32
}

impl DNAVariantRecordCaller {
    pub fn new(
        min_mapping_quality: MappingQuality,
        min_base_quality: BaseQuality,
        min_terminal_soft_clip_ins_len: u32
    ) -> Self {
        DNAVariantRecordCaller { min_mapping_quality, min_base_quality, min_terminal_soft_clip_ins_len }
    }
}

impl VariantRecordCaller for DNAVariantRecordCaller {
    type Input = AlignmentModel;

    fn call(&self, model: &AlignmentModel) -> Vec<VariantRecord> {
        let read_id: ReadID = model.get_read_id();
        let records: Vec<AlignmentModelRecord> = identify_alignment_model_records(model, None);

        // Step 1. Mismatches and insertions.
        let mut variant_records: Vec<VariantRecord> = identify_base_variant_records(
            read_id,
            records.iter(),
            self.min_mapping_quality,
            self.min_base_quality
        );

        // Step 2. Terminal soft clips, each a single-ended breakpoint at its anchor.
        variant_records.extend(identify_terminal_soft_clip_variant_records(
            model,
            self.min_mapping_quality,
            self.min_base_quality,
            self.min_terminal_soft_clip_ins_len,
            |chromosome_id, position, strand, operation, sequence| GraphOperation::new(
                chromosome_id,
                position,
                strand.clone(),
                operation,
                chromosome_id,
                position,
                strand,
                GraphOperationType::Noop,
                sequence,
                VariantType::Breakpoint
            )
        ));

        // Step 3. Deletions, breakpoints and translocations.
        variant_records.extend(identify_event_variant_records(
            read_id,
            &records,
            self.min_mapping_quality,
            self.min_base_quality
        ));

        sort_variant_records(&mut variant_records);
        variant_records
    }
}
