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


use exacto_caller::prelude::{
    AlignmentModelBaseContext,
    AlignmentModelBaseKind,
    AlignmentModelEventKind,
    AlignmentModelEventContext,
    GraphOperationView
};
use serde::{Serialize, Deserialize};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct AssembledTranscriptAlignmentRecord {
    pub index: u32,
    pub read_start: u32,
    pub read_end: u32,
    pub item_type: TranscriptAlignmentItemType,
    pub graph_operation_view: GraphOperationView,
    pub annotation: TranscriptAlignmentAnnotation
}

impl AssembledTranscriptAlignmentRecord {
    pub fn new(
        index: u32,
        read_start: u32,
        read_end: u32,
        item_type: TranscriptAlignmentItemType,
        graph_operation_view: GraphOperationView,
        annotation: TranscriptAlignmentAnnotation
    ) -> Self {
        Self {
            index: index,
            read_start: read_start,
            read_end: read_end,
            item_type: item_type,
            graph_operation_view: graph_operation_view,
            annotation: annotation
        }
    }
}


impl Clone for AssembledTranscriptAlignmentRecord {
    fn clone(&self) -> Self {
        AssembledTranscriptAlignmentRecord {
            index: self.index,
            read_start: self.read_start,
            read_end: self.read_end,
            item_type: self.item_type.clone(),
            graph_operation_view: self.graph_operation_view.clone(),
            annotation: self.annotation.clone()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum TranscriptAlignmentItemType {
    Base {
        kind: AlignmentModelBaseKind,
        context: AlignmentModelBaseContext
    },
    Event {
        kind: AlignmentModelEventKind,
        context: Option<AlignmentModelEventContext>
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptAlignmentAnnotation {
    pub position_1_annotation: Annotation,
    pub position_2_annotation: Annotation
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Annotation {
    pub gene_id: Option<Box<str>>,
    pub transcript_id: Option<Box<str>>,
    pub exon_id: Option<Box<str>>
}