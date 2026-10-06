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


/// Phred base quality of one read base.
pub type BaseQuality = u8;

/// ID of a read cluster.
pub type ClusterID = usize;

/// 1-based exon number within a transcript.
pub type ExonNumber = u16;

/// Mapping quality (MAPQ) of an alignment record.
pub type MappingQuality = u16;

/// Number of reads covering a position.
pub type ReadDepth = u32;

pub type ReadID = usize;

/// Read name (BAM QNAME).
pub type ReadName = Box<str>;

/// 0-based offset of a base in the read, in read orientation.
pub type ReadPosition = u32;

/// Number of reads supporting a variant. Same width as `ReadDepth`, since the two are compared and divided.
pub type ReadSupport = u32;

pub type ReferenceChromosomeID = u16;
pub type ReferenceChromosomeName = Box<str>;
pub type ReferenceChromosomeLength = u32;
pub type ReferenceGeneID = Box<str>;
pub type ReferenceGeneName = Box<str>;
pub type ReferenceTranscriptID = Box<str>;
pub type ReferenceExonID = Box<str>;

/// 1-based position on a contig of the reference genome.
pub type ReferencePosition = u32;

/// 0-based offset into the spliced reference transcript, 5' to 3'.
pub type ReferenceTranscriptPosition = u32;

/// 1-based splice junction number within a transcript.
pub type SpliceJunctionNumber = u16;

pub type VariantID = usize;
