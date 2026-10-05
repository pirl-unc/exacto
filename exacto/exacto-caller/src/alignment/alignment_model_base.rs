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
use serde::{Serialize, Deserialize};

use crate::prelude::*;


#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AlignmentModelBase {
    read_position: ReadPosition,
    nucleotide: Nucleotide,
    base_quality: BaseQuality,
    kind: AlignmentModelBaseKind,
    placement: AlignmentModelBasePlacement,
    reference_nucleotide: Option<Nucleotide>,
    deleted_reference_bases: Vec<Nucleotide>
}

impl AlignmentModelBase {
    pub fn set_read_position(&mut self, position: ReadPosition) {
        self.read_position = position;
    }

    pub fn replace_nucleotide(&mut self, nucleotide: Nucleotide, quality: BaseQuality) {
        if matches!(self.kind, AlignmentModelBaseKind::Match | AlignmentModelBaseKind::Mismatch) {
            let reference = self.reference_nucleotide.get_or_insert_with(|| self.nucleotide.clone());
            self.kind = if *reference == nucleotide { AlignmentModelBaseKind::Match } else { AlignmentModelBaseKind::Mismatch };
        }
        self.nucleotide = nucleotide;
        self.base_quality = quality;
    }

    pub fn new(
        read_position: ReadPosition,
        nucleotide: Nucleotide,
        base_quality: BaseQuality
    ) -> Self {
        Self {
            read_position,
            nucleotide,
            base_quality,
            kind: AlignmentModelBaseKind::Unaligned,
            placement: AlignmentModelBasePlacement::Unplaced,
            reference_nucleotide: None,
            deleted_reference_bases: Vec::new()
        }
    }

    pub fn get_reference_nucleotide(&self) -> Option<&Nucleotide> {
        self.reference_nucleotide.as_ref()
    }

    pub fn get_deleted_reference_bases(&self) -> &[Nucleotide] {
        &self.deleted_reference_bases
    }

    pub fn get_deletion_read_position(&self) -> Option<ReadPosition> {
        if self.deleted_reference_bases.is_empty() {
            return None;
        }
        match self.placement.get_coordinate()?.2 {
            Strand::Forward => Some(self.read_position + 1),
            _ => Some(self.read_position)
        }
    }

    pub fn get_base_quality(&self) -> BaseQuality {
        self.base_quality
    }

    pub fn get_kind(&self) -> &AlignmentModelBaseKind {
        &self.kind
    }

    pub fn get_nucleotide(&self) -> &Nucleotide {
        &self.nucleotide
    }

    pub fn get_read_position(&self) -> ReadPosition {
        self.read_position
    }

    pub fn get_placement(&self) -> &AlignmentModelBasePlacement {
        &self.placement
    }

    pub fn is_aligned(&self) -> bool {
        self.placement.is_placed()
            && matches!(
                self.kind,
                AlignmentModelBaseKind::Match
                    | AlignmentModelBaseKind::Mismatch
                    | AlignmentModelBaseKind::Insertion
            )
    }

    pub fn set_deleted_reference_bases(&mut self, bases: Vec<Nucleotide>) {
        self.deleted_reference_bases = bases;
    }

    pub fn set_kind(&mut self, kind: AlignmentModelBaseKind) {
        self.kind = kind;
    }

    pub fn set_reference_nucleotide(&mut self, nucleotide: Option<Nucleotide>) {
        self.reference_nucleotide = nucleotide;
    }

    pub fn place(
        &mut self,
        chromosome_id: ReferenceChromosomeID,
        position: ReferencePosition,
        strand: Strand,
        mapping_quality: MappingQuality
    ) {
        self.placement = AlignmentModelBasePlacement::Placed {
            chromosome_id,
            position,
            strand,
            mapping_quality
        };
    }
}


#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum AlignmentModelBasePlacement {
    Unplaced,
    Placed {
        chromosome_id: ReferenceChromosomeID,
        position: ReferencePosition,
        strand: Strand,
        mapping_quality: MappingQuality
    }
}

impl AlignmentModelBasePlacement {
    pub fn get_coordinate(&self) -> Option<(ReferenceChromosomeID, ReferencePosition, &Strand)> {
        match self {
            AlignmentModelBasePlacement::Unplaced => {
                None
            },
            AlignmentModelBasePlacement::Placed {
                chromosome_id,
                position,
                strand,
                mapping_quality,
                ..
            } => {
                Some((*chromosome_id, *position, strand))
            }
        }
    }

    pub fn get_mapping_quality(&self) -> Option<MappingQuality> {
        match self {
            AlignmentModelBasePlacement::Unplaced => {
                None
            },
            AlignmentModelBasePlacement::Placed {
                mapping_quality,
                ..
            } => {
                Some(*mapping_quality)
            }
        }
    }

    pub fn is_placed(&self) -> bool {
        matches!(self, AlignmentModelBasePlacement::Placed { .. })
    }
}
