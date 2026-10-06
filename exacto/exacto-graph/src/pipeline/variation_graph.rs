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
use exacto_core::prelude::*;
use exacto_core::log_info;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::str::FromStr;

use crate::prelude::*;


pub fn read_fasta_lengths(fasta_file: &str) -> Result<Vec<(Box<str>, u32)>, GraphError> {
    let file_error = |reason: String| GraphError::File { file: fasta_file.into(), reason: reason.into() };
    if !fasta_index_exists(fasta_file) {
        create_fai_file(fasta_file).map_err(|error| file_error(format!("no index could be made: {}", error)))?;
    }
    let fai_file: String = format!("{}.fai", fasta_file);
    let text: String = fs::read_to_string(&fai_file)
        .map_err(|error| file_error(format!("{} could not be read: {}", fai_file, error)))?;
    let mut lengths: Vec<(Box<str>, u32)> = Vec::new();
    for line in text.lines().filter(|line| !line.is_empty()) {
        let mut fields = line.split('\t');
        let name: &str = fields.next().unwrap_or("");
        let length: u32 = fields
            .next()
            .and_then(|field| field.parse::<u32>().ok())
            .ok_or_else(|| file_error(format!("{} has a line without a length: {}", fai_file, line)))?;
        lengths.push((name.into(), length));
    }
    Ok(lengths)
}

struct OperationColumns {
    chromosome_1: StringChunked,
    position_1: Int64Chunked,
    operation_1: StringChunked,
    strand_1: StringChunked,
    chromosome_2: StringChunked,
    position_2: Int64Chunked,
    operation_2: StringChunked,
    strand_2: StringChunked,
    sequence: StringChunked,
    num_cycles: Option<StringChunked>
}

fn read_column(df: &DataFrame, column: &str, data_type: &DataType) -> Result<Series, GraphError> {
    let column_error = |reason: String| GraphError::Column { column: column.into(), reason: reason.into() };
    df.column(column)
        .map_err(|_| column_error("not in the table".to_string()))?
        .as_materialized_series()
        .cast(data_type)
        .map_err(|error| column_error(format!("could not be read as {}: {}", data_type, error)))
}

fn read_string_column(df: &DataFrame, column: &str) -> Result<StringChunked, GraphError> {
    Ok(read_column(df, column, &DataType::String)?.str().unwrap().clone())
}

fn read_integer_column(df: &DataFrame, column: &str) -> Result<Int64Chunked, GraphError> {
    Ok(read_column(df, column, &DataType::Int64)?.i64().unwrap().clone())
}

impl OperationColumns {
    fn read(df: &DataFrame) -> Result<Self, GraphError> {
        Ok(OperationColumns {
            chromosome_1: read_string_column(df, "chromosome_1")?,
            position_1: read_integer_column(df, "position_1")?,
            operation_1: read_string_column(df, "operation_1")?,
            strand_1: read_string_column(df, "strand_1")?,
            chromosome_2: read_string_column(df, "chromosome_2")?,
            position_2: read_integer_column(df, "position_2")?,
            operation_2: read_string_column(df, "operation_2")?,
            strand_2: read_string_column(df, "strand_2")?,
            sequence: read_string_column(df, "sequence")?,
            num_cycles: if df.get_column_names().iter().any(|name| name.as_str() == "num_cycles") {
                Some(read_string_column(df, "num_cycles")?)
            } else {
                None
            }
        })
    }

    /// Reads row `i` into a GraphOperationView, or says what is wrong with it. `row` names the row
    /// in the error. Every position must be a base of its contig. A same-chromosome row that walks
    /// back (the U side below the D side) is a duplication and needs a number of copies: without
    /// num_cycles it takes `default_num_cycles`, and is an error when there is none.
    fn read_row(
        &self,
        i: usize,
        row: &str,
        chromosome_lengths_map: &HashMap<Box<str>, u32>,
        default_num_cycles: Option<u16>
    ) -> Result<GraphOperationView, GraphError> {
        let row_error = |reason: String| GraphError::Row { row: row.into(), reason: reason.into() };
        let read_side = |chromosomes: &StringChunked, positions: &Int64Chunked, operations: &StringChunked, strands: &StringChunked, side: u8|
            -> Result<(Box<str>, u32, GraphOperationType, Strand), GraphError> {
            let chromosome: &str = chromosomes
                .get(i)
                .ok_or_else(|| row_error(format!("chromosome_{} is missing", side)))?;
            let length: u32 = *chromosome_lengths_map
                .get(chromosome)
                .ok_or_else(|| row_error(format!("chromosome_{} {} is not in the FASTA file", side, chromosome)))?;
            let position: i64 = positions
                .get(i)
                .ok_or_else(|| row_error(format!("position_{} is missing", side)))?;
            if position < 1 || position > length as i64 {
                return Err(row_error(format!(
                    "position_{} {} is not a base of {} (1-{})", side, position, chromosome, length
                )));
            }
            let operation: &str = operations.get(i).unwrap_or("");
            let operation_type: GraphOperationType = GraphOperationType::from_str(operation)
                .map_err(|_| row_error(format!("operation_{} '{}' is not D, U or I", side, operation)))?;
            let strand: &str = strands.get(i).unwrap_or("");
            let strand: Strand = Strand::from_str(strand)
                .map_err(|_| row_error(format!("strand_{} '{}' is not +, - or *", side, strand)))?;
            Ok((chromosome.into(), position as u32, operation_type, strand))
        };
        let (chromosome_1, position_1, operation_type_1, strand_1) =
            read_side(&self.chromosome_1, &self.position_1, &self.operation_1, &self.strand_1, 1)?;
        let (chromosome_2, position_2, operation_type_2, strand_2) =
            read_side(&self.chromosome_2, &self.position_2, &self.operation_2, &self.strand_2, 2)?;

        let is_include_1: bool = operation_type_1 == GraphOperationType::Include;
        let is_include_2: bool = operation_type_2 == GraphOperationType::Include;
        if is_include_1 != is_include_2 {
            return Err(row_error("operation I must be on both sides".to_string()));
        }
        if is_include_1 && (chromosome_1 != chromosome_2 || position_1 > position_2) {
            return Err(row_error("an I row must span position_1 to position_2 on one chromosome".to_string()));
        }

        let is_cycle: bool = chromosome_1 == chromosome_2 && position_1 != position_2 && {
            let (lower, upper) = if position_1 < position_2 {
                (&operation_type_1, &operation_type_2)
            } else {
                (&operation_type_2, &operation_type_1)
            };
            *lower == GraphOperationType::Upstream && *upper == GraphOperationType::Downstream
        };
        let num_cycles: Option<u16> = match self.num_cycles.as_ref().and_then(|column| column.get(i)) {
            Some(value) if !value.trim().is_empty() => {
                let value: &str = value.trim();
                let num_cycles: u16 = value
                    .parse::<u16>()
                    .ok()
                    .or_else(|| {
                        // An integer column with gaps comes back from pandas as float (2.0)
                        value
                            .parse::<f64>()
                            .ok()
                            .filter(|v| v.fract() == 0.0 && *v >= 1.0 && *v <= u16::MAX as f64)
                            .map(|v| v as u16)
                    })
                    .filter(|num_cycles| *num_cycles >= 1)
                    .ok_or_else(|| row_error(format!("num_cycles '{}' is not a whole number of copies", value)))?;
                if !is_cycle {
                    return Err(row_error(format!("num_cycles {} is given for a row that is not a duplication", num_cycles)));
                }
                Some(num_cycles)
            },
            _ if is_cycle => Some(default_num_cycles.ok_or_else(|| {
                row_error("the row walks back (a duplication or back-splice) and needs num_cycles".to_string())
            })?),
            _ => None
        };

        Ok(GraphOperationView::new(
            &chromosome_1,
            position_1,
            operation_type_1,
            strand_1,
            &chromosome_2,
            position_2,
            operation_type_2,
            strand_2,
            self.sequence.get(i).unwrap_or(""),
            num_cycles
        ))
    }
}

/// The variants of one graph: the chromosomes they join, in FASTA order, and the variants by ID.
struct VariantGroup {
    chromosomes: Vec<Box<str>>,
    graph_operation_views: BTreeMap<usize, GraphOperationView>
}

/// Reads and checks every row of a genome variants table, then groups the variants by the
/// chromosomes they join. Only chromosomes with a variant are in a group. Groups are in FASTA order
/// of their first chromosome.
fn read_genome_variants(
    df_variants: &DataFrame,
    fasta_lengths: &Vec<(Box<str>, u32)>
) -> Result<Vec<VariantGroup>, GraphError> {
    let chromosome_lengths_map: HashMap<Box<str>, u32> = fasta_lengths.iter().cloned().collect();
    let chromosome_indices_map: HashMap<&str, u32> = fasta_lengths
        .iter()
        .enumerate()
        .map(|(i, (chromosome, _))| (chromosome.as_ref(), i as u32))
        .collect();

    // Step 1. Read every row
    log_info!("Started loading variants.");
    let col_variant_id: Int64Chunked = read_integer_column(df_variants, "variant_id")?;
    let columns: OperationColumns = OperationColumns::read(df_variants)?;
    let mut graph_operation_views: BTreeMap<usize, GraphOperationView> = BTreeMap::new();
    for i in 0..df_variants.height() {
        let variant_id: usize = match col_variant_id.get(i) {
            Some(variant_id) if variant_id >= 0 => variant_id as usize,
            _ => return Err(GraphError::Row {
                row: format!("row {}", i + 1).into(),
                reason: "variant_id is missing or negative".into()
            })
        };
        // exacto-caller writes a tandem duplication as a U/D row without a copy number: 2 copies
        let row: String = format!("variant_id {}", variant_id);
        let gov: GraphOperationView = columns.read_row(i, &row, &chromosome_lengths_map, Some(2))?;
        if graph_operation_views.insert(variant_id, gov).is_some() {
            return Err(GraphError::Row { row: row.into(), reason: "the variant_id is repeated".into() });
        }
    }
    log_info!("Finished loading {} variants.", graph_operation_views.len());

    // Step 2. Join the chromosomes each variant connects
    let mut uf: UnionFind = UnionFind::new();
    for gov in graph_operation_views.values() {
        uf.union(
            chromosome_indices_map[gov.get_chromosome_1()],
            chromosome_indices_map[gov.get_chromosome_2()]
        );
    }
    let mut clusters: Vec<BTreeSet<u32>> = uf.get_clusters();
    clusters.sort_by_key(|cluster| *cluster.first().unwrap());
    let mut group_indices_map: HashMap<u32, usize> = HashMap::new();
    let mut groups: Vec<VariantGroup> = Vec::new();
    for (group_index, cluster) in clusters.iter().enumerate() {
        for &chromosome_index in cluster.iter() {
            group_indices_map.insert(chromosome_index, group_index);
        }
        groups.push(VariantGroup {
            chromosomes: cluster.iter().map(|&i| fasta_lengths[i as usize].0.clone()).collect(),
            graph_operation_views: BTreeMap::new()
        });
    }

    // Step 3. Put each variant in the group of its chromosomes
    for (variant_id, gov) in graph_operation_views {
        let group_index: usize = group_indices_map[&chromosome_indices_map[gov.get_chromosome_1()]];
        groups[group_index].graph_operation_views.insert(variant_id, gov);
    }

    Ok(groups)
}

/// Builds the graph of one group. The chromosomes are split at every variant position first, so
/// no reference node is split while variants are added.
fn build_variant_group_graph(
    fasta_file: &str,
    group: &VariantGroup,
    chromosome_lengths_map: &HashMap<Box<str>, u32>,
    graph_type: &VarGraphTypes
) -> VarGraph {
    let reference_chromosomes_string: String = group.chromosomes
        .iter()
        .map(|s| s.as_ref())
        .collect::<Vec<&str>>()
        .join(LIST_SEPARATOR);

    // Collect variant breakpoints (positions where we want 1-bp reference nodes)
    let mut breakpoint_positions: HashMap<Box<str>, BTreeSet<u32>> = HashMap::new();
    for gov in group.graph_operation_views.values() {
        breakpoint_positions
            .entry(gov.get_chromosome_1().into())
            .or_insert_with(BTreeSet::new)
            .insert(gov.get_position_1());
        breakpoint_positions
            .entry(gov.get_chromosome_2().into())
            .or_insert_with(BTreeSet::new)
            .insert(gov.get_position_2());
    }

    // Build pre-split VarGraphReferenceNode lists for each chromosome
    let mut reference_nodes: Vec<Vec<VarGraphReferenceNode>> = Vec::new();
    for chromosome in group.chromosomes.iter() {
        let chr: &str = chromosome.as_ref();
        let chr_len: u32 = chromosome_lengths_map[chromosome];

        // Full sequence for this chromosome
        let seq: Box<str> = get_fasta_sequence(chr, 1, chr_len, fasta_file);

        // Breakpoints for this chromosome, sorted and deduped
        let breakpoints: Vec<u32> = breakpoint_positions
            .get(chr)
            .map(|set| set.iter().copied().collect())
            .unwrap_or_else(Vec::new);

        let mut chromosome_reference_nodes: Vec<VarGraphReferenceNode> = Vec::new();
        let mut curr_pos: u32 = 1;
        let mut seq_offset: usize = 0;
        for bp in breakpoints {
            // Region [curr_pos .. bp-1], if any
            if curr_pos < bp {
                let size: u32 = bp - curr_pos;
                let end: u32 = bp - 1;
                let slice_end: usize = seq_offset + size as usize;
                let subseq: &str = &seq[seq_offset..slice_end];
                chromosome_reference_nodes.push(VarGraphReferenceNode::new(
                    chr,
                    curr_pos,
                    end,
                    subseq,
                ));
                seq_offset = slice_end;
                curr_pos = bp;
            }

            // Single base [bp .. bp]
            if curr_pos == bp {
                let slice_end: usize = seq_offset + 1;
                let subseq: &str = &seq[seq_offset..slice_end];
                chromosome_reference_nodes.push(VarGraphReferenceNode::new(
                    chr,
                    bp,
                    bp,
                    subseq,
                ));
                seq_offset = slice_end;
                curr_pos = bp + 1;
            }
        }

        // Tail region [curr_pos .. chr_len], if any
        if curr_pos <= chr_len {
            let size: u32 = chr_len - curr_pos + 1;
            let slice_end: usize = seq_offset + size as usize;
            let subseq: &str = &seq[seq_offset..slice_end];
            chromosome_reference_nodes.push(VarGraphReferenceNode::new(
                chr,
                curr_pos,
                chr_len,
                subseq
            ));
        }
        reference_nodes.push(chromosome_reference_nodes);
    }

    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(reference_nodes);

    // Add variants, in variant ID order
    let mut count: usize = 0;
    let total: usize = group.graph_operation_views.len();
    for (variant_id, graph_op_view) in group.graph_operation_views.iter() {
        count += 1;
        if count % 10_000 == 0 {
            log_info!("VarGraph for {reference_chromosomes_string} processed {count}/{total}.");
        }
        let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(
            *variant_id,
            graph_op_view.clone()
        );
        vargraph.add_variant_node(variant_node);
    }

    if *graph_type == VarGraphTypes::Individual {
        // Only stitch adjacent variants if we're building an individual (personalized) graph
        vargraph.stitch_adjacent_variants();
    }

    vargraph
}

/// Builds one graph for each group of chromosomes that variants join, in FASTA order. Chromosomes
/// without a variant get no graph.
pub fn build_genome_variation_graph(
    fasta_file: &str,
    df_variants: &DataFrame,
    graph_type: VarGraphTypes,
    num_threads: usize
) -> Result<Vec<VarGraph>, GraphError> {
    let fasta_lengths: Vec<(Box<str>, u32)> = read_fasta_lengths(fasta_file)?;
    let chromosome_lengths_map: HashMap<Box<str>, u32> = fasta_lengths.iter().cloned().collect();
    let groups: Vec<VariantGroup> = read_genome_variants(df_variants, &fasta_lengths)?;

    log_info!("Started adding variants to the variation graphs.");
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let vargraphs: Vec<VarGraph> = thread_pool.install(|| {
        groups
            .par_iter()
            .map(|group| build_variant_group_graph(fasta_file, group, &chromosome_lengths_map, &graph_type))
            .collect()
    });
    log_info!("Finished adding variants to the variation graphs.");

    Ok(vargraphs)
}

/// The sequences build-genome-var-graph writes, in a fixed order, each handed to `emit` with
/// whether it carries variants.
///
/// First the variant sequences, group by group in FASTA order: paths that between them carry every
/// variant allele (`VarGraph::find_genome_paths`). Then, when `include_reference_sequences` is set,
/// the reference sequence of every contig in the FASTA file, so each variant's reference allele is
/// written too. `num_threads` groups are built at a time, and each graph is dropped once its
/// sequences are made.
pub fn find_genome_variation_graph_sequences<F>(
    fasta_file: &str,
    df_variants: &DataFrame,
    graph_type: VarGraphTypes,
    include_reference_sequences: bool,
    num_threads: usize,
    mut emit: F
) -> Result<(), GraphError>
where
    F: FnMut(Box<str>, bool) -> Result<(), GraphError>
{
    let fasta_lengths: Vec<(Box<str>, u32)> = read_fasta_lengths(fasta_file)?;
    let chromosome_lengths_map: HashMap<Box<str>, u32> = fasta_lengths.iter().cloned().collect();
    let groups: Vec<VariantGroup> = read_genome_variants(df_variants, &fasta_lengths)?;

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    for batch in groups.chunks(num_threads.max(1)) {
        let sequences_list: Vec<Vec<Box<str>>> = thread_pool.install(|| {
            batch
                .par_iter()
                .map(|group| {
                    let vargraph: VarGraph = build_variant_group_graph(
                        fasta_file,
                        group,
                        &chromosome_lengths_map,
                        &graph_type
                    );
                    vargraph
                        .find_genome_paths(&vargraph.get_variant_node_ids().into_iter().collect(), &Default::default())
                        .iter()
                        .map(|path| path.get_sequence())
                        .collect()
                })
                .collect()
        });
        for sequence in sequences_list.into_iter().flatten() {
            emit(sequence, true)?;
        }
    }

    if include_reference_sequences {
        for (chromosome, length) in fasta_lengths.iter() {
            emit(get_fasta_sequence(chromosome, 1, *length, fasta_file), false)?;
        }
    }

    Ok(())
}

/// The rows of one transcript model by index, or the first of its rows that could not be read.
struct TranscriptModelRows {
    name: Box<str>,
    rows: Result<Vec<(usize, GraphOperationView)>, GraphError>
}

/// Reads and checks every row of a transcript model alignments table (the
/// `_exacto_assembled_transcript_model_alignments.tsv` of call-rna-transcript-vars), grouped by
/// `assembled_transcript_name` in the order the names first appear. A missing column is an error
/// for the table; a row that cannot be read is an error for its model only. A row that walks back
/// needs num_cycles: the number of laps a read makes is not a copy number.
fn read_transcript_models(
    df_transcript_structures: &DataFrame,
    chromosome_lengths_map: &HashMap<Box<str>, u32>
) -> Result<Vec<TranscriptModelRows>, GraphError> {
    let col_name: StringChunked = read_string_column(df_transcript_structures, "assembled_transcript_name")?;
    let col_index: Int64Chunked = read_integer_column(df_transcript_structures, "index")?;
    let columns: OperationColumns = OperationColumns::read(df_transcript_structures)?;
    let mut models: Vec<TranscriptModelRows> = Vec::new();
    let mut model_indices_map: HashMap<Box<str>, usize> = HashMap::new();
    for i in 0..df_transcript_structures.height() {
        let (name, index): (&str, usize) = match (col_name.get(i), col_index.get(i)) {
            (Some(name), Some(index)) if index >= 0 => (name, index as usize),
            _ => return Err(GraphError::Row {
                row: format!("row {}", i + 1).into(),
                reason: "assembled_transcript_name or index is missing, or index is negative".into()
            })
        };
        let row: String = format!("transcript {} index {}", name, index);
        let model_index: usize = *model_indices_map.entry(name.into()).or_insert_with(|| {
            models.push(TranscriptModelRows { name: name.into(), rows: Ok(Vec::new()) });
            models.len() - 1
        });
        if let Ok(rows) = models[model_index].rows.as_mut() {
            match columns.read_row(i, &row, chromosome_lengths_map, None) {
                Ok(gov) => rows.push((index, gov)),
                Err(error) => models[model_index].rows = Err(error)
            }
        }
    }
    Ok(models)
}

/// Builds the graph of one transcript model. Its rows are taken in index order; a repeated index
/// is an error for this model only.
fn build_transcript_model_graph(
    fasta_file: &str,
    model: &TranscriptModelRows,
    chromosome_lengths_map: &HashMap<Box<str>, u32>,
    graph_type: &VarGraphTypes
) -> Result<VarGraph, GraphError> {
    let mut graph_op_views: Vec<(usize, GraphOperationView)> = match &model.rows {
        Ok(rows) => rows.clone(),
        Err(error) => return Err(GraphError::Transcript { name: model.name.clone(), reason: error.to_string().into() })
    };
    graph_op_views.sort_by_key(|(index, _)| *index);
    if let Some(pair) = graph_op_views.windows(2).find(|pair| pair[0].0 == pair[1].0) {
        return Err(GraphError::Transcript {
            name: model.name.clone(),
            reason: format!("index {} is repeated", pair[0].0).into()
        });
    }

    // Build per-chromosome (start, end) spans from graph operations
    let mut chromosome_spans: BTreeMap<Box<str>, (u32, u32)> = BTreeMap::new();
    for (_, gov) in graph_op_views.iter() {
        for (chr, pos) in [
            (gov.get_chromosome_1(), gov.get_position_1()),
            (gov.get_chromosome_2(), gov.get_position_2()),
        ] {
            let entry = chromosome_spans
                .entry(chr.into())
                .or_insert((pos, pos));
            entry.0 = entry.0.min(pos);
            entry.1 = entry.1.max(pos);
        }
    }

    // Load only the needed region per chromosome
    let mut reference_nodes: Vec<Vec<VarGraphReferenceNode>> = Vec::new();
    for (chromosome, (min_pos, max_pos)) in chromosome_spans.iter() {
        let chromosome_length: u32 = chromosome_lengths_map[chromosome];
        let start: u32 = (*min_pos).max(1);
        let end: u32 = (*max_pos).min(chromosome_length);
        let sequence: Box<str> = get_fasta_sequence(chromosome, start, end, fasta_file);
        let reference_node = VarGraphReferenceNode::new(chromosome, start, end, &*sequence);
        reference_nodes.push(vec![reference_node]);
    }

    // Create a VarGraph from the reference nodes
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(reference_nodes);

    // Add variants
    for (variant_id, graph_op_view) in graph_op_views.iter() {
        let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(
            *variant_id,
            graph_op_view.clone()
        );
        vargraph.add_variant_node(variant_node);
    }

    // Stitch adjacent variants
    if *graph_type == VarGraphTypes::Individual {
        // Only stitch adjacent variants if we're building an individual (personalized) graph
        vargraph.stitch_adjacent_variants();
    }

    Ok(vargraph)
}

/// Builds one graph per transcript model, in the order the models appear in the table. A problem
/// with the table (a column, a row, a chromosome) is an error; a model that cannot be built is
/// logged and left out.
pub fn build_transcriptome_variation_graph(
    fasta_file: &str,
    df_transcript_structures: &DataFrame,
    graph_type: VarGraphTypes,
    num_threads: usize
) -> Result<Vec<(Box<str>, VarGraph)>, GraphError> {
    let chromosome_lengths_map: HashMap<Box<str>, u32> = read_fasta_lengths(fasta_file)?.into_iter().collect();
    let models: Vec<TranscriptModelRows> = read_transcript_models(df_transcript_structures, &chromosome_lengths_map)?;
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let results: Vec<Result<VarGraph, GraphError>> = thread_pool.install(|| {
        models
            .par_iter()
            .map(|model| build_transcript_model_graph(fasta_file, model, &chromosome_lengths_map, &graph_type))
            .collect()
    });
    let mut vargraphs: Vec<(Box<str>, VarGraph)> = Vec::new();
    for (model, result) in models.into_iter().zip(results) {
        match result {
            Ok(vargraph) => vargraphs.push((model.name, vargraph)),
            Err(error) => log::warn!("Skipped: {}", error)
        }
    }
    Ok(vargraphs)
}

/// The sequences build-transcriptome-var-graph writes: one per transcript model, named by its
/// `assembled_transcript_name`, in the order the models appear in the table. Models are built
/// `batch_size` at a time (0 = all at once) and dropped once their sequence is made. A model that
/// does not make exactly one path is logged and skipped; the others are still written.
pub fn find_transcriptome_variation_graph_sequences<F>(
    fasta_file: &str,
    df_transcript_structures: &DataFrame,
    graph_type: VarGraphTypes,
    num_threads: usize,
    batch_size: usize,
    mut emit: F
) -> Result<(), GraphError>
where
    F: FnMut(Box<str>, Box<str>) -> Result<(), GraphError>
{
    let chromosome_lengths_map: HashMap<Box<str>, u32> = read_fasta_lengths(fasta_file)?.into_iter().collect();
    let models: Vec<TranscriptModelRows> = read_transcript_models(df_transcript_structures, &chromosome_lengths_map)?;
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let batch_size: usize = if batch_size == 0 { models.len().max(1) } else { batch_size };
    let mut num_skipped: usize = 0;
    for batch in models.chunks(batch_size) {
        let results: Vec<Result<Box<str>, GraphError>> = thread_pool.install(|| {
            batch
                .par_iter()
                .map(|model| {
                    let vargraph: VarGraph = build_transcript_model_graph(
                        fasta_file,
                        model,
                        &chromosome_lengths_map,
                        &graph_type
                    )?;
                    let path: VarGraphPath = vargraph.find_transcript_path(
                        &model.name,
                        &vargraph.get_variant_node_ids().into_iter().collect()
                    )?;
                    Ok(path.get_sequence())
                })
                .collect()
        });
        for (model, result) in batch.iter().zip(results) {
            match result {
                Ok(sequence) => emit(model.name.clone(), sequence)?,
                Err(error) => {
                    log::warn!("Skipped: {}", error);
                    num_skipped += 1;
                }
            }
        }
    }
    if num_skipped > 0 {
        log::warn!("Skipped {} of {} transcript models.", num_skipped, models.len());
    }
    Ok(())
}


#[cfg(test)]
#[path = "../tests/pipeline/variation_graph.rs"]
mod tests;