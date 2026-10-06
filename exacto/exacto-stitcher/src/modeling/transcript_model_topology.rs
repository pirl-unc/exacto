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

use crate::prelude::TranscriptModelTopology;


pub(crate) fn classify_topology(model: &TranscriptModel) -> TranscriptModelTopology {
    let mut has_fusion: bool = false;
    for (read_position_1, read_position_2, context) in get_annotated_events(model) {
        match context {
            AlignmentModelEventContext::BackSplicing => {
                return TranscriptModelTopology::BackSplicing
            },
            AlignmentModelEventContext::FusionGene => {
                has_fusion = true
            },
            _ => {}
        }
    }
    if has_fusion {
        TranscriptModelTopology::Fusion
    } else {
        TranscriptModelTopology::Linear
    }
}


fn get_annotated_events(
    model: &TranscriptModel
) -> impl Iterator<Item = (u32, u32, &AlignmentModelEventContext)> {
    let annotation: &TranscriptModelAnnotation = model.get_annotation();
    model.get_alignment_model().get_events().keys()
        .filter_map(move |&(read_position_1, read_position_2)| {
            annotation.get_event(read_position_1, read_position_2)?
                .get_context()
                .as_ref()
                .map(|context| (read_position_1, read_position_2, context))
        })
}
