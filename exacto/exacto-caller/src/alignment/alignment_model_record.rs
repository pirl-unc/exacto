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
use exacto_core::prelude::*;
use serde::{Deserialize, Serialize};
use std::fmt;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct AlignmentModelRecord {
    start: ReadPosition,                                 // read position
    end: ReadPosition,                                   // read position
    sequence: Box<str>,
    base_quality_scores: Vec<BaseQuality>,
    record_type: AlignmentModelRecordType,
    kind: AlignmentModelKind,
    context: Option<AlignmentModelContext>,
    chromosome_1: ReferenceChromosomeID,
    position_1: ReferencePosition,
    operation_1: GraphOperationType,
    strand_1: Strand,
    mapping_quality_1: MappingQuality,
    chromosome_2: ReferenceChromosomeID,
    position_2: ReferencePosition,
    operation_2: GraphOperationType,
    strand_2: Strand,
    mapping_quality_2: MappingQuality,
    gene_id_1: Option<ReferenceGeneID>,
    transcript_id_1: Option<ReferenceTranscriptID>,
    exon_id_1: Option<ReferenceExonID>,
    gene_id_2: Option<ReferenceGeneID>,
    transcript_id_2: Option<ReferenceTranscriptID>,
    exon_id_2: Option<ReferenceExonID>,
    skipped: Option<Vec<Vec<ReferenceBase>>>
}

impl AlignmentModelRecord {
    pub fn new(
        start: ReadPosition,
        end: ReadPosition,
        sequence: &str,
        base_quality_scores: Vec<BaseQuality>,
        record_type: AlignmentModelRecordType,
        kind: AlignmentModelKind,
        context: Option<AlignmentModelContext>,
        chromosome_1: ReferenceChromosomeID,
        position_1: ReferencePosition,
        operation_1: GraphOperationType,
        strand_1: Strand,
        mapping_quality_1: MappingQuality,
        chromosome_2: ReferenceChromosomeID,
        position_2: ReferencePosition,
        operation_2: GraphOperationType,
        strand_2: Strand,
        mapping_quality_2: MappingQuality,
        gene_id_1: Option<ReferenceGeneID>,
        transcript_id_1: Option<ReferenceTranscriptID>,
        exon_id_1: Option<ReferenceExonID>,
        gene_id_2: Option<ReferenceGeneID>,
        transcript_id_2: Option<ReferenceTranscriptID>,
        exon_id_2: Option<ReferenceExonID>,
        skipped: Option<Vec<Vec<ReferenceBase>>>
    ) -> Self {
        assert_eq!(
            base_quality_scores.len(), sequence.len(), 
            "Base quality scores vector length must be equal to sequence length"
        );
        if record_type == AlignmentModelRecordType::Base {
            assert_eq!(
                strand_1, strand_2, 
                "Strands must be equal for base records"
            );
        }
        Self {
            start: start,
            end: end,
            sequence: sequence.into(),
            base_quality_scores: base_quality_scores,
            record_type: record_type,
            kind: kind,
            context: context,
            chromosome_1: chromosome_1,
            position_1: position_1,
            operation_1: operation_1,
            strand_1: strand_1,
            mapping_quality_1: mapping_quality_1,
            chromosome_2: chromosome_2,
            position_2: position_2,
            operation_2: operation_2,
            strand_2: strand_2,
            mapping_quality_2: mapping_quality_2,
            gene_id_1: gene_id_1,
            transcript_id_1: transcript_id_1,
            exon_id_1: exon_id_1,
            gene_id_2: gene_id_2,
            transcript_id_2: transcript_id_2,
            exon_id_2: exon_id_2,
            skipped: skipped
        }
    }
    
    pub fn get_start(&self) -> ReadPosition {
        self.start
    }

    pub fn get_end(&self) -> ReadPosition {
        self.end
    }

    pub fn get_sequence(&self) -> &Box<str> {
        &self.sequence
    }

    pub fn get_standardized_sequence(&self) -> Box<str> {
        if self.strand_1 == Strand::Forward {
            self.sequence.to_string().to_uppercase().into_boxed_str()
        } else {
            reverse_complement(&self.sequence).to_string().to_uppercase().into_boxed_str()
        }
    }

    pub fn get_base_quality_scores(&self) -> &Vec<BaseQuality> {
        &self.base_quality_scores
    }

    pub fn get_record_type(&self) -> &AlignmentModelRecordType {
        &self.record_type
    }

    pub fn get_kind(&self) -> &AlignmentModelKind {
        &self.kind
    }

    pub fn get_context(&self) -> &Option<AlignmentModelContext> {
        &self.context
    }

    pub fn get_chromosome_1(&self) -> ReferenceChromosomeID {
        self.chromosome_1
    }

    pub fn get_position_1(&self) -> ReferencePosition {
        self.position_1
    }

    pub fn get_operation_1(&self) -> &GraphOperationType {
        &self.operation_1
    }

    pub fn get_strand_1(&self) -> &Strand {
        &self.strand_1
    }

    pub fn get_mapping_quality_1(&self) -> MappingQuality {
        self.mapping_quality_1
    }

    pub fn get_chromosome_2(&self) -> ReferenceChromosomeID {
        self.chromosome_2
    }

    pub fn get_position_2(&self) -> ReferencePosition {
        self.position_2
    }

    pub fn get_operation_2(&self) -> &GraphOperationType {
        &self.operation_2
    }

    pub fn get_strand_2(&self) -> &Strand {
        &self.strand_2
    }

    pub fn get_mapping_quality_2(&self) -> MappingQuality {
        self.mapping_quality_2
    }

    pub fn get_gene_id_1(&self) -> &Option<ReferenceGeneID> {
        &self.gene_id_1
    }

    pub fn get_transcript_id_1(&self) -> &Option<ReferenceTranscriptID> {
        &self.transcript_id_1
    }

    pub fn get_exon_id_1(&self) -> &Option<ReferenceExonID> {
        &self.exon_id_1
    }

    pub fn get_gene_id_2(&self) -> &Option<ReferenceGeneID> {
        &self.gene_id_2
    }

    pub fn get_transcript_id_2(&self) -> &Option<ReferenceTranscriptID> {
        &self.transcript_id_2
    }

    pub fn get_exon_id_2(&self) -> &Option<ReferenceExonID> {
        &self.exon_id_2
    }

    pub fn get_skipped(&self) -> &Option<Vec<Vec<ReferenceBase>>> {
        &self.skipped
    }

    pub fn get_skipped_string(&self, chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>) -> String {
        if self.skipped.is_none() {
            "".to_string()
        } else {
            let mut skipped_reference_bases: Vec<Box<str>> = Vec::new();
            for cluster in self.skipped.as_ref().unwrap().iter() {
                let first_skipped_base: &ReferenceBase = cluster.first().unwrap();
                let last_skipped_base: &ReferenceBase = cluster.last().unwrap();
                let gene_id_1: ReferenceGeneID = first_skipped_base.reference_gene_id.as_deref().unwrap_or("").into();
                let transcript_id_1: ReferenceTranscriptID = first_skipped_base.reference_transcript_id.as_deref().unwrap_or("").into();
                let exon_id_1: ReferenceExonID = first_skipped_base.reference_exon_id.as_deref().unwrap_or("").into();
                let gene_id_2: ReferenceGeneID = last_skipped_base.reference_gene_id.as_deref().unwrap_or("").into();
                let transcript_id_2: ReferenceTranscriptID = last_skipped_base.reference_transcript_id.as_deref().unwrap_or("").into();
                let exon_id_2: ReferenceExonID = last_skipped_base.reference_exon_id.as_deref().unwrap_or("").into();
                let encoding: Box<str> = format!(
                    "{}:{}:S:{}|{}:{}:S:{}|{}:{}:{}|{}:{}:{}",
                    chromosome_names_map.get_by_right(&first_skipped_base.reference_chromosome_id).unwrap(),
                    first_skipped_base.reference_position,
                    first_skipped_base.reference_strand.as_str(),
                    chromosome_names_map.get_by_right(&last_skipped_base.reference_chromosome_id).unwrap(),
                    last_skipped_base.reference_position,
                    last_skipped_base.reference_strand.as_str(),
                    gene_id_1,
                    transcript_id_1,
                    exon_id_1,
                    gene_id_2,
                    transcript_id_2,
                    exon_id_2
                ).into();
                skipped_reference_bases.push(encoding);
            }
            skipped_reference_bases.join(LIST_SEPARATOR)
        }
    }
}

impl fmt::Display for AlignmentModelRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "AlignmentStructureRecord:")?;
        writeln!(f, "  start-end: {}-{}", self.start, self.end)?;
        writeln!(f, "  sequence: {}", self.sequence)?;
        writeln!(f, "  type: {:?}", self.record_type)?;
        writeln!(f, "  kind: {:?}", self.kind)?;
        writeln!(f, "  chromosome_1: {} position_1: {} strand_1: {:?} mq1: {}",
                 self.chromosome_1, self.position_1, self.strand_1, self.mapping_quality_1)?;
        writeln!(f, "  chromosome_2: {} position_2: {} strand_2: {:?} mq2: {}",
                 self.chromosome_2, self.position_2, self.strand_2, self.mapping_quality_2)?;
        if let Some(skipped) = &self.skipped {
            writeln!(f, "  skipped: {} groups", skipped.len())?;
        } else {
            writeln!(f, "  skipped: None")?;
        }
        Ok(())
    }
}

impl Clone for AlignmentModelRecord {
    fn clone(&self) -> Self {
        AlignmentModelRecord {
            start: self.start,
            end: self.end,
            sequence: self.sequence.clone(),
            base_quality_scores: self.base_quality_scores.clone(),
            record_type: self.record_type.clone(),
            kind: self.kind.clone(),
            context: self.context.clone(),
            chromosome_1: self.chromosome_1,
            position_1: self.position_1,
            operation_1: self.operation_1.clone(),
            strand_1: self.strand_1.clone(),
            mapping_quality_1: self.mapping_quality_1,
            chromosome_2: self.chromosome_2,
            position_2: self.position_2,
            operation_2: self.operation_2.clone(),
            strand_2: self.strand_2.clone(),
            mapping_quality_2: self.mapping_quality_2,
            gene_id_1: self.gene_id_1.clone(),
            transcript_id_1: self.transcript_id_1.clone(),
            exon_id_1: self.exon_id_1.clone(),
            gene_id_2: self.gene_id_2.clone(),
            transcript_id_2: self.transcript_id_2.clone(),
            exon_id_2: self.exon_id_2.clone(),
            skipped: self.skipped.clone()
        }
    }
}
